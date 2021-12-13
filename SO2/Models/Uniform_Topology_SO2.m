%% Author: Yunpeng Shi
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% generation of the synthetic data
%%------------------------------------------------
%% Input Parameters: 
%% n: the number of the graph nodes
%% p: the probability of connecting a pair of vertices. G([n],E) is Erdos-Renyi graph G(n,p).
%% q: the probability of corrupting an edge
%% sigma: the noise level (>0)
%% model (optional): choose 'uniform' or 'self-consistent'. The default choice is 'uniform'.

%% Output:
%% model_out.AdjMat: n by n adjacency matrix of the generated graph
%% model_out.Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). edge_num is the number of edges.
%% model_out.thetaij: vector that stores the given relative rotation angles
%% model_out.thetaij_orig: vector that stores the ground truth relative rotations
%% model_out.theta_orig: vector that stores the ground truth absolute rotations
%% model_out.ErrVec: the true corruption level of each edge
%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing", Foundations of Computational Mathematics, 2021.


function[model_out]=Uniform_Topology_SO2(n,p,q,crpt_type)
    if ~exist('crpt_type','var')
        crpt_type = 'uniform';
    end
    
    G = rand(n,n) < p;
    G = tril(G,-1);
    % generate adjacency matrix
    [Ind_j, Ind_i] = find(G==1);
    Ind = [Ind_i, Ind_j];
    m = length(Ind_i);
    AdjMat = sparse(Ind_i,Ind_j,1,n,n);
    AdjMat = AdjMat + AdjMat';


    %generate rotation matrices
    theta_orig = 2*pi*rand(1,n);
    
    
    thetaij_orig = zeros(1,m);
    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        thetaij_orig(k)=theta_orig(i)-theta_orig(j);
    end
    
    thetaij_orig = mod(thetaij_orig+2*pi, 2*pi);
    thetaij = thetaij_orig;
    corrIndLog = rand(1,m)<=q;
    % indices of corrupted edges
    corrInd=find(corrIndLog);
    
    thetaij = thetaij + sigma * randn(1,m);            
    thetaij = mod(thetaij+2*pi, 2*pi);
    
    theta_corr = 2*pi*rand(1,n);
    
    if strcmp(crpt_type,'uniform')  
        thetaij(corrInd) = 2*pi*rand(1,length(corrInd));
    else 
        for k = corrInd
            i=Ind_i(k); j=Ind_j(k); 
            thetaij(k) = theta_corr(i)-theta_corr(j) + sigma * randn(1);
        end  
    end
    
    thetaij = mod(thetaij+2*pi, 2*pi);

    ErrVec = mod(thetaij-thetaij_orig+2*pi, 2*pi)/pi;
    ErrVec = min(ErrVec,2-ErrVec);

    
    
    model_out.AdjMat = AdjMat;
    model_out.Ind = Ind;
    model_out.thteaij = thetaij;
    model_out.thteaij_orig = thetaij_orig;
    model_out.theta_orig = theta_orig;
    model_out.ErrVec = ErrVec;
    
    
end
