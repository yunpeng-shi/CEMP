%% Author: Yunpeng Shi

%%------------------------------------------------
%% Cycle-Edge Message Passing for Z2 Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% zij: vector that stores the given relative rotations corresponding to Ind
%% parameters.beta_init: initial reweighting parameter beta for CEMP
%% parameters.beta_max: the maximal reweighting parameter beta for CEMP
%% parameters.rate: beta is updated by beta = beta*rate until it hits beta_max



%% Output:
%% z_est: Estimated group element
%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing", Foundations of Computational Mathematics, 2021.



function z_est = CEMP_GCW_Z2(Ind,zij,parameters)

    %CEMP parameters
    beta_init = parameters.beta_init;
    beta_max = parameters.beta_max;
    rate = parameters.rate;
            
    % building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    AdjMat = sparse(Ind_i,Ind_j,1,n,n); % Adjacency matrix
    AdjMat = full(AdjMat + AdjMat');
        
    zijMat = sparse(Ind_i,Ind_j,zij,n,n);
    zijMat = zijMat + zijMat';
    Weights = AdjMat;
    iter = 1;
    beta = beta_init;
    while beta <= beta_max
        S = zijMat.*Weights;
        S2 = S^2;
        Weights2 = Weights^2;
        Weights2_inv = zeros(n,n);
        Weights2_inv(Weights2>0) = 1./Weights2(Weights2>0);
        SMat = (0.5-0.5*(S2.*Weights2_inv).*zijMat).*AdjMat;
        Weights = exp(-beta.*SMat).*AdjMat;
        fprintf('Reweighting Iteration %d Completed!\n', iter)   
        iter=iter+1;
        beta = beta*rate;
    end
    
    disp('CEMP DONE');
    
    Weights = exp(-beta/rate.*SMat).*AdjMat;
    row_sum = sum(Weights,2);
    Weights(row_sum==0,:)=1; % avoid zero rowsums
    row_sum(row_sum==0)=n;
    Weights = diag(1./row_sum)*Weights;
    zijW = zijMat.*Weights;
    [V,~] = eigs(zijW,1,'la');
    z_est = sign(V);
                    
end
