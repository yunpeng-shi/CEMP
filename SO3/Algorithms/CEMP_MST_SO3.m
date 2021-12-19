%% Author: Yunpeng Shi
%%------------------------------------------------
%% Cycle-Edge Message Passing for Rotation Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind
%% parameters.beta_init: initial reweighting parameter beta for CEMP
%% parameters.beta_max: the maximal reweighting parameter beta for CEMP
%% parameters.rate: beta is updated by beta = beta*rate until it hits beta_max


%% 
%% Output:
%% R_est: Estimated absolute rotations (3x3xn)

%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing" arXiv preprint, 2019




function R_est = CEMP_MST_SO3(Ind,RijMat,parameters)
    %CEMP parameters
    beta_init = parameters.beta_init;
    beta_max = parameters.beta_max;
    rate = parameters.rate;
            
    % building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    n=max(Ind,[],'all');
    m=size(Ind_i,1);
    AdjMat = sparse(Ind_i,Ind_j,1,n,n); % Adjacency matrix
    AdjMat = full(AdjMat + AdjMat');
    % start CEMP iterations as initialization   
    CoDeg = (AdjMat*AdjMat).*AdjMat;
    CoDeg((CoDeg==0)&(AdjMat>0))=-1;
    CoDeg_low = tril(CoDeg,-1);
    CoDeg_vec = CoDeg_low(:);
    CoDeg_vec(CoDeg_vec==0)=[];
    CoDeg_pos_ind = find(CoDeg_vec>0);
    CoDeg_vec_pos = CoDeg_vec(CoDeg_pos_ind);
    CoDeg_zero_ind = find(CoDeg_vec<0);
    cum_ind = [0;cumsum(CoDeg_vec_pos)];
    m_pos = length(CoDeg_pos_ind);
    m_cycle = cum_ind(end);
    
    Ind_ij = zeros(1,m_cycle);
    Ind_jk = zeros(1,m_cycle);
    Ind_ki = zeros(1,m_cycle);
    
    RijMat4d = zeros(3,3,n,n);
    
    % construct edge index matrix (for 2d-to-1d index conversion)
    for l = 1:m
        i=Ind_i(l);j=Ind_j(l);
        RijMat4d(:,:,i,j)=RijMat(:,:,l);
        RijMat4d(:,:,j,i)=(RijMat(:,:,l))';
        IndMat(i,j)=l;
        IndMat(j,i)=-l;
    end
   
    Rjk0Mat = zeros(3,3,m_cycle);
    Rki0Mat = zeros(3,3,m_cycle);
    
    for l = 1:m_pos
        IJ = CoDeg_pos_ind(l);
        i=Ind_i(IJ); j=Ind_j(IJ);
        CoInd_ij= find(AdjMat(:,i).*AdjMat(:,j));
        Ind_ij((cum_ind(l)+1):cum_ind(l+1)) =  IJ;
        Ind_jk((cum_ind(l)+1):cum_ind(l+1)) =  abs(IndMat(j,CoInd_ij));
        Ind_ki((cum_ind(l)+1):cum_ind(l+1)) =  abs(IndMat(CoInd_ij,i));
        Rjk0Mat(:,:,(cum_ind(l)+1):cum_ind(l+1)) =  RijMat4d(:,:,j,CoInd_ij);
        Rki0Mat(:,:,(cum_ind(l)+1):cum_ind(l+1)) =  RijMat4d(:,:,CoInd_ij,i);
    end

    Rij0Mat = RijMat(:,:,Ind_ij);

    
    
    disp('compute R cycle')

    R_cycle0 = zeros(3,3,m_cycle);
    R_cycle = zeros(3,3,m_cycle);
    for j = 1:3
      R_cycle0 = R_cycle0 + bsxfun(@times,Rij0Mat(:,j,:),Rjk0Mat(j,:,:));
    end
    
    for j = 1:3
      R_cycle = R_cycle + bsxfun(@times,R_cycle0(:,j,:),Rki0Mat(j,:,:));
    end

    R_trace = reshape(R_cycle(1,1,:)+R_cycle(2,2,:)+R_cycle(3,3,:),[1,m_cycle]);
    S0_long = abs(acos((R_trace-1)./2))/pi;
    S0_vec = ones(1,m);
    
    Weight_vec = ones(1,m_cycle);
    S0_weight = S0_long.*Weight_vec;
    
    for l=1:m_pos
        IJ = CoDeg_pos_ind(l);
        S0_vec(IJ) = sum(S0_weight((cum_ind(l)+1):cum_ind(l+1)))/sum(Weight_vec((cum_ind(l)+1):cum_ind(l+1)));
    end
    
        disp('Initialization completed!')
    
   
        disp('Reweighting Procedure Started ...')
    
    iter = 0;
    
    SVec = S0_vec;
    beta = beta_init;
    while beta <= beta_max
        Sjk = SVec(Ind_jk);
        Ski = SVec(Ind_ki);
        S_sum = Ski+Sjk;
        
        Weight_vec = exp(-beta*S_sum);
        S0_weight = S0_long.*Weight_vec;
    
        for l=1:m_pos
            IJ = CoDeg_pos_ind(l);
            SVec(IJ) = sum(S0_weight((cum_ind(l)+1):cum_ind(l+1)))/sum(Weight_vec((cum_ind(l)+1):cum_ind(l+1)));
        end
        
        fprintf('Reweighting Iteration %d Completed!\n',iter)   
        iter = iter+1;
        % parameter controling the decay rate of reweighting function
        beta = beta*rate;
    end

    disp('Building minimum spanning tree ...')
    SMatij = sparse(Ind_j,Ind_i,SVec+1,n,n);
    [MST,~]=graphminspantree(SMatij);
    AdjTree = logical(MST+MST');
    
    % compute Ri by multiplying Rij along the spanning tree
    rootnodes = 1;
    added=zeros(1,n);
    R_est = zeros(3,3,n);
    R_est(:,:,rootnodes)=eye(3);
    added(rootnodes)=1;
    newroots = [];
            
    while sum(added)<n
        for node_root = rootnodes
            leaves = find((AdjTree(node_root,:).*(1-added))==1);
            newroots = [newroots, leaves];
            for node_leaf=leaves
                edge_leaf = IndMat(node_leaf,node_root);
                if edge_leaf>0
                    R_est(:,:,node_leaf)=RijMat(:,:,abs(edge_leaf))*R_est(:,:,node_root);
                else
                    R_est(:,:,node_leaf)=(RijMat(:,:,abs(edge_leaf)))'*R_est(:,:,node_root);
                end
                added(node_leaf)=1;
            end
        end
        rootnodes = newroots;
    end


end
