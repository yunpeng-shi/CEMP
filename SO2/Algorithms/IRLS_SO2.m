%% Author: Yunpeng Shi
%%------------------------------------------------
%% Iteratively Reweighted Least Squares for Angular Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% thetaij: vector that stores the given relative rotations corresponding to Ind

%% Output:
%% theta_est: Estimated group element


function theta_est = IRLS_SO2(Ind,thetaij,niter)
              
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    m = length(Ind_i);
    n=max(Ind,[],'all');       
    AdjMat = sparse(Ind_i,Ind_j,1,n,n); % Adjacency matrix
    AdjMat = full(AdjMat + AdjMat');
    thetaijMat1 = sparse(Ind_i,Ind_j,thetaij,n,n);
    thetaijMat1 = thetaijMat1 - thetaijMat1';
    aijMat = exp(1i*thetaijMat1).*AdjMat;
    %%% Spectral 
    [V,~] = eigs(aijMat,1,'la');
    theta_sp = mod(real(-1i*log(V./abs(V)))+2*pi, 2*pi);

    theta_irls = theta_sp;
    thetaij_irls = zeros(1,m);
    irls_stop = 1;
    iter = 1;
    while iter <= niter && irls_stop >1e-3
        for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        thetaij_irls(k)=theta_irls(i)-theta_irls(j);
        end
        thetaij_irls = mod(thetaij_irls+2*pi, 2*pi);
        resid_vec = mod(thetaij_irls-thetaij+2*pi, 2*pi)/pi;
        resid_vec = min(resid_vec,2-resid_vec) + 1e-4;
        weight_vec = 1./resid_vec;
        Weights = sparse(Ind_i, Ind_j, weight_vec, n, n);
        Weights = Weights + Weights';

        Weights = diag(1./sum(Weights,2))*Weights;
        aijW = aijMat.*Weights;
        [V,~] = eigs(aijW,1,'la');
        theta_irls_new = mod(real(-1i*log(V./abs(V)))+2*pi, 2*pi);
        diff_vec = mod(theta_irls_new-theta_irls+2*pi, 2*pi)/pi;
        irls_stop = mean(min(diff_vec,2-diff_vec));

        theta_irls =  theta_irls_new;
        iter = iter +1;

    end

    theta_est = theta_irls;
                    
end
