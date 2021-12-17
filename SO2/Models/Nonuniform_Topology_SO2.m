%% Author: Yunpeng Shi

%%------------------------------------------------
%% generation of the synthetic data with nonuniform topology
%%------------------------------------------------
%% Input Parameters: 
%% n: the number of the graph nodes
%% p: the probability of connecting a pair of vertices. G([n],E) is Erdos-Renyi graph G(n,p).
%% p_node_crpt: the probability of corrupting a node
%% p_node_edge: the probability of corrupting an edge
%% sigma_in: the noise level for inliers
%% sigma_out: the noise level for outliers
%% crpt_type (optional): choose 'uniform' or 'self-consistent', or 'adv'. The default choice is 'uniform'.

%% Output:
%% model_out.AdjMat: n by n adjacency matrix of the generated graph
%% model_out.Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). edge_num is the number of edges.
%% model_out.thetaij: vector that stores the given relative angles
%% model_out.thetaij_orig: vector that stores the ground truth relative angles
%% model_out.theta_orig: vector that stores the ground truth absolute angle
%% model_out.ErrVec: the true corruption level of each edge
%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing", Foundations of Computational Mathematics, 2021.

function[model_out]=Nonuniform_Topology_SO2(n,p, p_node_crpt, p_edge_crpt, sigma_in, sigma_out, crpt_type)
    if ~exist('crpt_type','var')
        crpt_type = 'uniform';
    end
    
    G = rand(n,n) < p;
    G = tril(G,-1);
    AdjMat = G + G'; 
    % generate adjacency matrix
    [Ind_j, Ind_i] = find(G==1);
    Ind = [Ind_i, Ind_j];
    Ind_full = [Ind_j, Ind_i;Ind_i, Ind_j];
    m = length(Ind_i);
    theta_orig = 2*pi*rand(1,n);
    
    IndMat = zeros(n,n);
    thetaij_orig = zeros(1,m);
    for k = 1:m
        i=Ind_i(k); j=Ind_j(k); 
        thetaij_orig(k)=theta_orig(i)-theta_orig(j);
        IndMat(i,j)=k;
        IndMat(j,i)=-k;
    end
    
    thetaij_orig = mod(thetaij_orig+2*pi, 2*pi);
    thetaij = thetaij_orig;
    
    theta_crpt = 2*pi*rand(1,n);
    node_crpt = randperm(n);
    n_node_crpt = floor(n*p_node_crpt);
    node_crpt = node_crpt(1:n_node_crpt);
    crptInd = false(1,m);
    theta_0 = 2*pi*rand(1,m);
    
    
    for i = node_crpt
        neighbor_cand = Ind_full(Ind_full(:,1)==i,2);
        neighbor_cand = reshape(neighbor_cand, 1, length(neighbor_cand));
        neighbor_crpt = randperm(length(neighbor_cand));
        n_neighbor = floor(p_edge_crpt * length(neighbor_cand));
        neighbor_crpt =    neighbor_crpt(1:n_neighbor);
        neighbor_crpt = neighbor_cand(neighbor_crpt);

        for j = neighbor_crpt 

            k = IndMat(i,j);
            crptInd(abs(k))=1; 
            

            if strcmp(crpt_type,'uniform')
                if k>0
                    thetaij(k)= theta_0(k);
                else
                    thetaij(-k)= -theta_0(-k);
                end

            end

            if strcmp(crpt_type,'self-consistent')
                if k>0
                    thetaij(k)=theta_crpt(i) - theta_crpt(j);
                else
                    thetaij(-k)=theta_crpt(j)-theta_crpt(i);
                end
            end

            if strcmp(crpt_type,'adv')

                if k>0
                    thetaij(k)=theta_crpt(i) - theta_orig(j);
                else
                    thetaij(-k)=theta_orig(j)-theta_crpt(i);
                end
                
            end 

        end

    end


    noiseInd = ~crptInd;
    % indices of corrupted edges
    thetaij(noiseInd)= ...
    thetaij(noiseInd)+sigma_in*randn(1,sum(noiseInd)); % add noise

    thetaij(crptInd)= ...
    thetaij(crptInd)+sigma_out*randn(1,sum(crptInd)); % add noise
   

    thetaij = mod(thetaij+2*pi, 2*pi);
    
    ErrVec = mod(thetaij-thetaij_orig+2*pi, 2*pi)/pi;
    ErrVec = min(ErrVec,2-ErrVec);
    
    
    model_out.AdjMat = AdjMat;
    model_out.Ind = Ind;
    model_out.thetaij = thetaij;
    model_out.thetaij_orig = thetaij_orig;
    model_out.theta_orig = theta_orig;
    model_out.ErrVec = ErrVec;
    
end