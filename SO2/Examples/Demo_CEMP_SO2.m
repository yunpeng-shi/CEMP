n=200; p=1; q=0.7; sigma=0; crpt_type='uniform';

model_out = Uniform_Topology_SO2(n,p,q,sigma,crpt_type);

Ind = model_out.Ind; % matrix of edge indices (m by 2)
thetaij = model_out.thetaij; % given corrupted and noisy group ratios
ErrVec = model_out.ErrVec; % ground truth corruption levels

% set CEMP defult parameters
parameters.beta_init = 1;
parameters.beta_max = 40;
parameters.rate = 1.2;

% run CEMP
SVec = CEMP_SO2(Ind,thetaij,parameters);

%visualize sij^* and sij,t, ideally one should see the line "y=x"
plot(ErrVec,SVec,'b.');
title('Scatter Plot of s_{ij}^* v.s. s_{ij,T}')
xlabel('s_{ij}^*') 
ylabel('s_{ij,T}') 
