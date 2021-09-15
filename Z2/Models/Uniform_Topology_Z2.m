%% Author: Yunpeng Shi
%% generation of the synthetic data with uniform graph topology
%%------------------------------------------------
%% Input Parameters: 
%% n: the number of the graph nodes
%% p: the probability of connecting a pair of vertices. G([n],E) is Erdos-Renyi graph G(n,p).
%% q: the probability of corrupting an edge
%% crpt_type (optional): choose 'uniform' or 'adv'. The default choice is 'uniform'.

%% Output:
%% model_out.AdjMat: n by n adjacency matrix of the generated graph
%% model_out.Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). edge_num is the number of edges.
%% model_out.zij: vector that stores the given group ratios
%% model_out.zij_orig: vector that stores the ground truth group ratios
%% model_out.z_orig: vector that stores the ground truth absolute group elements
%% model_out.ErrVec: the true corruption level of each edge
%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing", Foundations of Computational Mathematics, 2021.


function[model_out]=Uniform_Topology_Z2(n,p,q,crpt_type)
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
    z_orig = 2*(randi([0,1],1,n)-0.5);
    zij_orig = zeros(1,m);
    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        zij_orig(k)=z_orig(i)*z_orig(j);
    end

    zij = zij_orig;
    crptIndLog = rand(1,m)<=q;
    % indices of corrupted edges
    crptInd=find(crptIndLog);
    z_crpt = 2*(randi([0,1],1,n)-0.5);
  
    if strcmp(crpt_type,'uniform')
        zij(crptInd) = 2*(randi([0,1],1,length(crptInd))-0.5);
    else
        for k = crptInd
        i=Ind_i(k); j=Ind_j(k); 
        zij(k) = -z_crpt(i)*z_crpt(j);
        end  
    end
        
    ErrVec = (1-zij.*zij_orig)*0.5;
    
    
    model_out.AdjMat = AdjMat;
    model_out.Ind = Ind;
    model_out.zij = zij;
    model_out.zij_orig = zij_orig;
    model_out.z_orig = z_orig;
    model_out.ErrVec = ErrVec;
    
    
end
