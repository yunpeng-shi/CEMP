%% Author: Yunpeng Shi
%%------------------------------------------------
%% Cycle-Edge Message Passing for Z2 Synchronization
%% This is a equivalent implementation of CEMP using matrix multiplications. 
%% Its speed is usually faster than its original version "CEMP_Z2.m"
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% zij: vector that stores the given relative rotations corresponding to Ind
%% parameters.beta_init: initial reweighting parameter beta for CEMP
%% parameters.beta_max: the maximal reweighting parameter beta for CEMP
%% parameters.rate: beta is updated by beta = beta*rate until it hits beta_max


%% Output:
%% SVec: Estimated corruption levels of all edges

%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing", Foundations of Computational Mathematics, 2021.




function[SVec] = CEMP_fast_Z2(Ind,zij,parameters)
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
    
    SVec = SMat(sub2ind([n,n], Ind_i, Ind_j));

        

end
