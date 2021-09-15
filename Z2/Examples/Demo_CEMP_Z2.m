n=200; p=1; q=0.6; crpt_type='uniform';

model_out = Uniform_Topology_Z2(n,p,q,crpt_type);

Ind = model_out.Ind; % matrix of edge indices (m by 2)
zij = model_out.zij; % given corrupted and noisy group ratios
ErrVec = model_out.ErrVec; % ground truth corruption levels

% set CEMP defult parameters
parameters.beta_init = 1;
parameters.beta_max = 40;
parameters.rate = 1.2;

% run CEMP
SVec = CEMP_fast_Z2(Ind,zij,parameters);

%visualize sij^* and sij,t, ideally one should see two points at (0,0) and
%(1,1).
plot(ErrVec,SVec,'b.');
title('Scatter Plot of s_{ij}^* v.s. s_{ij,T}')
xlabel('s_{ij}^*') 
ylabel('s_{ij,T}') 
