%% Author: Yunpeng Shi
%% 
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind


%% Output:
%% R_est: Estimated rotations (3x3xn)



function R_est = IRLS_SO3(Ind,RijMat,niter)

 % building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    m=size(Ind_i,1);
    d=3;
    mat_size = ones(1,n)*d;
    cum_ind = [0,cumsum(mat_size)];
    Rij_blk = zeros(n*d);
    for k = 1:m
       i = Ind_i(k); j=Ind_j(k);
       Rij_blk((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1))= RijMat(:,:,k);    
    end

    Rij_blk = Rij_blk+Rij_blk';

    %%% Spectral 

    [V,~] = eigs(Rij_blk,d, 'la');

    V(:,1) = V(:,1)*sign(det(V(1:d,:))); % ensure det = 1
    R_est = zeros(d,d,n);
    for i=1:n
       Ri = V((cum_ind(i)+1):cum_ind(i+1), :); 
       [Ur,~,Vr] = svd(Ri);
       S0 = diag([ones(1,d-1),det(Ur*Vr')]);
       R_est(:,:,i) = Ur*S0*Vr';
    end

    R_irls = R_est;
    Rij_irls = zeros(d,d,m);
    Rij_irls_new = zeros(d,d,m);
    irls_stop = 1;
    iter = 1;
    while iter <= niter && irls_stop >1e-3
        for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        Rij_irls(:,:,k)=R_irls(:,:,i)*(R_irls(:,:,j))';
        end
        
        MSE_res = zeros(d,d,m);
        for j = 1:d
            MSE_res = MSE_res + bsxfun(@times,Rij_irls(:,j,:),RijMat(:,j,:));
        end
    
    
        res_trace = zeros(1,m);
    
        for j = 1:d
            res_trace = res_trace + reshape(MSE_res(j,j,:),[1,m]);
        end
    
        resid_vec = sqrt(abs((1-res_trace/d)*0.5)) + 1e-4;
        weight_vec = 1./resid_vec;
        Weights = sparse(Ind_i, Ind_j, weight_vec, n, n);
        Weights = Weights + Weights';
        
        Weights = diag(1./sum(Weights,2))*Weights;
        Weights = kron(Weights, ones(d));    
        RijW = Rij_blk.*Weights;
        
        [V,~] = eigs(RijW,d,'la');
        V(:,1) = V(:,1)*sign(det(V(1:d,:))); % ensure det = 1
        R_est = zeros(d,d,n);
        for i=1:n
           Ri = V((cum_ind(i)+1):cum_ind(i+1), :); 
           [Ur,~,Vr] = svd(Ri);
           S0 = diag([ones(1,d-1),det(Ur*Vr')]);
           R_est(:,:,i) = Ur*S0*Vr';
    
        end
        R_irls_new = R_est;
        for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        Rij_irls_new(:,:,k)=R_irls_new(:,:,i)*(R_irls_new(:,:,j))';
        end
        diff_R = zeros(d,d,m);
        for j = 1:d
          diff_R = diff_R + bsxfun(@times,Rij_irls(:,j,:),Rij_irls_new(:,j,:));
        end
    
    
        diff_trace = zeros(1,m);
    
        for j = 1:d
            diff_trace = diff_trace + reshape(diff_R(j,j,:),[1,m]);
        end
    
    
        irls_stop = sqrt(mean(abs(1-diff_trace/d)*0.5));
        
        R_irls =  R_irls_new;
        iter = iter+1;
       
    end



end