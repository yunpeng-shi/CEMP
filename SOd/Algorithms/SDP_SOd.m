%% Author: Yunpeng Shi
%%------------------------------------------------
%% SDP relaxation for SO(d) synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind


%% Output:
%% R_est: Estimated rotations (3x3xn)



function R_est = SDP_SOd(Ind,RijMat)

    d = size(RijMat,1);
 % building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    m=size(Ind_i,1);
    mat_size = ones(1,n)*d;
    cum_ind = [0,cumsum(mat_size)];
    Rij_blk = zeros(n*d);
    for k = 1:m
       i = Ind_i(k); j=Ind_j(k);
       Rij_blk((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1))= RijMat(:,:,k);    
    end

    Rij_blk = Rij_blk+Rij_blk';

    %%% Spectral 

    cvx_begin sdp quiet
        variable X_sdp(n*d,n*d)  semidefinite
        diag(X_sdp) == 1
        %X_sdp == semidefinite(n)
        maximize (trace(Rij_blk*X_sdp))
    cvx_end
    [V,~] = eigs(X_sdp,d,'la');
    V(:,1) = V(:,1)*sign(det(V(1:d,:))); % ensure det = 1
    R_est = zeros(d,d,n);
    for i=1:n
       Ri = V((cum_ind(i)+1):cum_ind(i+1), :); 
       [Ur,~,Vr] = svd(Ri);
       S0 = diag([ones(1,d-1),det(Ur*Vr')]);
       R_est(:,:,i) = Ur*S0*Vr';
   
    end

end
