%% Author: Yunpeng Shi
%%------------------------------------------------
%% Cycle-Edge Message Passing for Angular Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% thetaij: vector that stores the given relative rotations corresponding to Ind

%% Output:
%% theta_est: Estimated group element


function theta_est = SDP_SO2(Ind,thetaij)
              
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');       
    thetaijMat1 = sparse(Ind_i,Ind_j,thetaij,n,n);
    thetaijMat1 = thetaijMat1 - thetaijMat1';
    aijMat = exp(1i*thetaijMat1).*AdjMat;

    cvx_begin sdp quiet
        variable X_sdp(n,n) hermitian complex semidefinite
        diag(X_sdp) == 1
        %X_sdp == semidefinite(n)
        maximize (trace(aijMat*X_sdp))
    cvx_end
    [V,~] = eigs(X_sdp,1,'la');
    theta_est = mod(real(-1i*log(V./abs(V)))+2*pi, 2*pi);
                    
end

              
