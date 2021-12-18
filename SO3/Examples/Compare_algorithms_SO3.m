rng(1);
% parameters with uniform topology
% n=100; p=0.5; q=0.6; sigma=0; crpt_type='uniform';

% use the following commented line to generate data with uniform topology
% model_out = Uniform_Topology_SO3(n,p,q,sigma,crpt_type);


% The following code is for nonuniform topology. This is a more malicious
% scenario where corrupted edges have cluster behavior (so local coruption
% level can be extremely high)

% parameters with nonuniform topology
n=100; p=0.5; p_node_crpt=0.5; p_edge_crpt=0.75; sigma_in=0; sigma_out=2; crpt_type='self-consistent';
model_out = Nonuniform_Topology_SO3(n,p, p_node_crpt,p_edge_crpt, sigma_in, sigma_out, crpt_type);

Ind = model_out.Ind; % matrix of edge indices (m by 2)
RijMat = model_out.RijMat; % given corrupted and noisy relative rotations
ErrVec = model_out.ErrVec; % ground truth corruption levels
R_orig = model_out.R_orig; % ground truth absolute rotations
Rij_orig = model_out.Rij_orig; % ground truth relative rotations


% set CEMP defult parameters
parameters.beta_init = 1;
parameters.beta_max = 40;
parameters.rate = 1.2;


% run Spectral
R_SP = Spectral_SO3(Ind, RijMat);

% run SDP  
R_SDP = SDP_SO3(Ind, RijMat); %commented this out if cvx is not installed

% run IRLS
niter = 100;
R_IRLS = IRLS_SO3(Ind, RijMat, niter);

% run CEMP+MST
R_CEMP_MST = CEMP_MST_SO3(Ind, RijMat, parameters);

% run CEMP+GCW
R_CEMP_GCW = CEMP_GCW_SO3(Ind, RijMat, parameters);

[mean_error_SP] = evaluate_error_SO3(R_SP, Rij_orig, Ind);
[mean_error_SDP] = evaluate_error_SO3(R_SDP, Rij_orig, Ind);
[mean_error_IRLS] = evaluate_error_SO3(R_IRLS, Rij_orig, Ind);
[mean_error_CEMP_MST] = evaluate_error_SO3(R_CEMP_MST, Rij_orig, Ind);
[mean_error_CEMP_GCW] = evaluate_error_SO3(R_CEMP_GCW, Rij_orig, Ind);


% Report estimation error
sz = [5 2];
varTypes = {'string','double'};
varNames = {'Algorithms','MeanError'};
Results = table('Size',sz,'VariableTypes',varTypes, 'VariableNames',varNames);
Results(1,:)={'Spectral', mean_error_SP};
Results(2,:)={'SDP', mean_error_SDP};
Results(3,:)={'IRLS', mean_error_IRLS};
Results(4,:)={'CEMP+MST', mean_error_CEMP_MST};
Results(5,:)={'CEMP+GCW', mean_error_CEMP_GCW};



Results
