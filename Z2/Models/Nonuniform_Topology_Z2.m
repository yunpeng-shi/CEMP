%% Author: Yunpeng Shi

%%------------------------------------------------
%% generation of the synthetic data with nonuniform topology
%%------------------------------------------------
%% Input Parameters: 
%% n: the number of the graph nodes
%% p: the probability of connecting a pair of vertices. G([n],E) is Erdos-Renyi graph G(n,p).
%% p_node_crpt: the probability of corrupting a node
%% p_node_edge: the probability of corrupting an edge
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

function[model_out]=Nonuniform_Topology_Z2(n,p, p_node_crpt,p_edge_crpt, crpt_type)
    if ~exist('crpt_type','var')
        crpt_type = 'uniform';
    end
    
    G = rand(n,n) < p;
    G = tril(G,-1);
    % generate adjacency matrix
    [Ind_j, Ind_i] = find(G==1);
    Ind = [Ind_i,Ind_j];
    Ind_full = [Ind_j, Ind_i;Ind_i, Ind_j];
    m = length(Ind_i);

    AdjMat = sparse(Ind_i,Ind_j,1,n,n);
    AdjMat = AdjMat + AdjMat';
    %generate rotation matrices
    z_orig = 2*(randi([0,1],1,n)-0.5);
    zij_orig = zeros(1,m);
    IndMat = zeros(n,n);
 
    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        zij_orig(k)=z_orig(i)*z_orig(j);
        IndMat(i,j)=k;
        IndMat(j,i)=-k;
    end

    zij = zij_orig;

    z_crpt = 2*(randi([0,1],1,n)-0.5);
    

    node_crpt = randperm(n);
    n_node_crpt = floor(n*p_node_crpt);
    node_crpt = node_crpt(1:n_node_crpt);

    for i = node_crpt
        neighbor_cand = Ind_full(Ind_full(:,1)==i,2);
        neighbor_cand = reshape(neighbor_cand, 1, length(neighbor_cand));
        neighbor_crpt = randperm(length(neighbor_cand));
        n_neighbor = floor(p_edge_crpt * length(neighbor_cand));
        neighbor_crpt =    neighbor_crpt(1:n_neighbor);
        neighbor_crpt = neighbor_cand(neighbor_crpt);

        for j = neighbor_crpt 

            k = abs(IndMat(i,j));
          
            if strcmp(crpt_type,'adv')
                zij(k) = -z_crpt(i)*z_crpt(j);
            else
                zij(k) = 2*(randi([0,1])-0.5);
            end
            
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
