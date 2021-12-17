%% Author: Yunpeng Shi
%%------------------------------------------------
%% Cycle-Edge Message Passing for Rotation Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind
%% CEMP_parameters.max_iter: the number of iterations of CEMP
%% CEMP_parameters.reweighting: the sequence of reweighting parameter beta_t
%% 
%% Output:
%% R_est: Estimated corruption levels of all edges

%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing" arXiv preprint, 2019




function R_est = CEMP_GCW_SO3(Ind,RijMat,parameters)
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
    % store pairwise directions in 3 by n by n tensor
    % construct edge index matrix (for 2d-to-1d index conversion)
    for l = 1:m
        i=Ind_i(l);j=Ind_j(l);
        RijMat4d(:,:,i,j)=RijMat(:,:,l);
        RijMat4d(:,:,j,i)=(RijMat(:,:,l))';
        IndMat(i,j)=l;
        IndMat(j,i)=l;
    end
   
    Rjk0Mat = zeros(3,3,m_cycle);
    Rki0Mat = zeros(3,3,m_cycle);
    
    for l = 1:m_pos
        IJ = CoDeg_pos_ind(l);
        i=Ind_i(IJ); j=Ind_j(IJ);
        CoInd_ij= find(AdjMat(:,i).*AdjMat(:,j));
        Ind_ij((cum_ind(l)+1):cum_ind(l+1)) =  IJ;
        Ind_jk((cum_ind(l)+1):cum_ind(l+1)) =  IndMat(j,CoInd_ij);
        Ind_ki((cum_ind(l)+1):cum_ind(l+1)) =  IndMat(CoInd_ij,i);
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
        % compute weight matrix (nsample by m)
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

    d=3;
    mat_size = ones(1,n)*d;
    cum_ind = [0,cumsum(mat_size)];
    Rij_blk = zeros(n*d);
    for k = 1:m
       i = Ind_i(k); j=Ind_j(k);
       Rij_blk((cum_ind(i)+1):cum_ind(i+1), (cum_ind(j)+1):cum_ind(j+1))= RijMat(:,:,k);    
    end
    
    Rij_blk = Rij_blk+Rij_blk';
    
    %%% Spectral 
    beta_T = beta/rate;
    SMat = sparse(Ind_i, Ind_j, SVec, n, n);
    SMat = SMat + SMat';
    Weights = exp(-beta_T.*SMat).*AdjMat;
    Weights = diag(1./sum(Weights,2))*Weights;
    Weights = kron(Weights, ones(d));    
    RijW = Rij_blk.*Weights;
    clear 'Rij_blk';
    
    
    [V,~] = eigs(RijW,d,'la');
    V(:,1) = V(:,1)*sign(det(V(1:d,:))); % ensure det = 1
    R_est = zeros(d,d,n);
    for i=1:n
       Ri = V((cum_ind(i)+1):cum_ind(i+1), :); 
       [Ur,~,Vr] = svd(Ri);
       S0 = diag([ones(1,d-1),det(Ur*Vr')]);
       R_est(:,:,i) = Ur*S0*Vr';
       
    end


end
