%% Author: Yunpeng Shi
%%------------------------------------------------
%% Eigenvector method for Z2 synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% zij: vector that stores the given relative rotations corresponding to Ind


%% Output:
%% z_est: Estimated group element


function z_est = Spectral_Z2(Ind,zij)
              
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');       
    zijMat = sparse(Ind_i,Ind_j,zij,n,n);
    zijMat = zijMat + zijMat';
    [V,~] = eigs(zijMat,1,'la');
    z_est = sign(V);
                    
end
