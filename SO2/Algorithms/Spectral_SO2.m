%% Author: Yunpeng Shi
%%------------------------------------------------
%% Eigenvector method for Angular Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% thetaij: vector that stores the given relative rotations corresponding to Ind

%% Output:
%% theta_est: Estimated group element


function theta_est = Spectral_SO2(Ind,thetaij)
              
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    
    n=max(Ind,[],'all');       
    AdjMat = sparse(Ind_i,Ind_j,1,n,n); % Adjacency matrix
    AdjMat = full(AdjMat + AdjMat');
    thetaijMat1 = sparse(Ind_i,Ind_j,thetaij,n,n);
    thetaijMat1 = thetaijMat1 - thetaijMat1';
    aijMat = exp(1i*thetaijMat1).*AdjMat;
    %%% Spectral 
    [V,~] = eigs(aijMat,1,'la');
    theta_est = mod(real(-1i*log(V./abs(V)))+2*pi, 2*pi);
                    
end
