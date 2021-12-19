%% Author: Yunpeng Shi
%%------------------------------------------------
%% Cycle-Edge Message Passing for Angular Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j) that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% thetaij: vector that stores the given relative rotations corresponding to Ind
%% parameters.beta_init: initial reweighting parameter beta for CEMP
%% parameters.beta_max: the maximal reweighting parameter beta for CEMP
%% parameters.rate: beta is updated by beta = beta*rate until it hits beta_max


%% Output:
%% theta_est: Estimated absolute group elements (angles)

%% Reference
%% [1] Gilad Lerman and Yunpeng Shi. "Robust Group Synchronization via Cycle-Edge Message Passing", Foundations of Computational Mathematics, 2021.




function[theta_est] = CEMP_GCW_SO2(Ind,thetaij,parameters)
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
    % Matrix of codegree:
    % CoDeg(i,j) = 0 if i and j are not connected, otherwise,
    % CoDeg(i,j) = # of vertices that are connected to both i and j
    CoDeg = (AdjMat*AdjMat).*AdjMat;
    CoDeg((CoDeg==0)&(AdjMat>0))=-1;
    CoDeg_low = tril(CoDeg,-1);
    CoDeg_vec = CoDeg_low(:);
    CoDeg_vec(CoDeg_vec==0)=[];
    CoDeg_pos_ind = find(CoDeg_vec>0);
    CoDeg_vec_pos = CoDeg_vec(CoDeg_pos_ind);
    cum_ind = [0;cumsum(CoDeg_vec_pos)];
    m_pos = length(CoDeg_pos_ind);
    m_cycle = cum_ind(end);

    Ind_ij = zeros(1,m_cycle);
    Ind_jk = zeros(1,m_cycle);
    Ind_ki = zeros(1,m_cycle);


  
    IndMat = sparse(Ind_i,Ind_j,(1:m),n,n);
    IndMat = IndMat+IndMat';

    thetaijMat = sparse(Ind_i,Ind_j,thetaij,n,n);
    thetaji = mod(-thetaij, 2*pi);
    thetaijMat = thetaijMat + sparse(Ind_j, Ind_i, thetaji,n,n);
    thetaijMat1 = sparse(Ind_i,Ind_j,thetaij,n,n);
    thetaijMat1 = thetaijMat1 - thetaijMat1';

    
    jkvec = zeros(1,m_cycle);
    kivec = zeros(1,m_cycle);

    for l = 1:m_pos
        IJ = CoDeg_pos_ind(l);
        i=Ind_i(IJ); j=Ind_j(IJ);
        CoInd_ij= find(AdjMat(:,i).*AdjMat(:,j));
        Ind_ij((cum_ind(l)+1):cum_ind(l+1)) =  IJ;
        Ind_jk((cum_ind(l)+1):cum_ind(l+1)) =  IndMat(j,CoInd_ij);
        Ind_ki((cum_ind(l)+1):cum_ind(l+1)) =  IndMat(CoInd_ij,i);
        jkvec((cum_ind(l)+1):cum_ind(l+1)) =  thetaijMat(j,CoInd_ij);
        kivec((cum_ind(l)+1):cum_ind(l+1)) =  thetaijMat(CoInd_ij,i);
    end
   
    
    ijvec = thetaij(Ind_ij);
    


    disp('compute cycle-inconsistencies')
    
    theta_cycle = mod(ijvec + jkvec +kivec +6*pi, 2*pi)/pi;
    S0_long = min(theta_cycle,2-theta_cycle);
    S0_vec = ones(1,m);

    


    Weight_vec = ones(1,m_cycle);
    S0_weight = S0_long.*Weight_vec;

    for l=1:m_pos
        IJ = CoDeg_pos_ind(l);
        S0_vec(IJ) = sum(S0_weight((cum_ind(l)+1):cum_ind(l+1)))/sum(Weight_vec((cum_ind(l)+1):cum_ind(l+1)));
    end

    disp('Initialization completed!')

    
    disp('Reweighting Procedure Started ...')

    

   SVec = S0_vec;
   iter=1;
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

        beta = beta*rate;
        fprintf('Reweighting Iteration %d Completed!\n',iter)   
        iter = iter+1;

    end

    disp('CEMP stage Completed!')

    SMat = sparse(Ind_i, Ind_j, SVec, n, n);
    SMat = SMat + SMat';
    aijMat = exp(1i*thetaijMat1).*AdjMat;

    Weights = exp(-beta/rate.*SMat).*AdjMat;
    Weights = diag(1./sum(Weights,2))*Weights;
    aijW = aijMat.*Weights;
    [V,~] = eigs(aijW,1,'la');
    theta_est = mod(real(-1i*log(V./abs(V)))+2*pi, 2*pi);        

end
