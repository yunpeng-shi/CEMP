rng(1)
% parameters with uniform topology
% n =200; p=0.5; q=0.8; sigma=0; crpt_type='uniform';
% model_out = Uniform_Topology_SO2(n,p,q,sigma,crpt_type);

% parameters with nonuniform topology
n =200; p=0.5; p_node_crpt=0.8; p_edge_crpt=0.75; sigma_in=0; sigma_out=2;crpt_type='self-consistent';


% The following code is for nonuniform topology. This is a more malicious
% scenario where corrupted edges have cluster behavior (so local coruption
% level can be extremely high)
model_out = Nonuniform_Topology_SO2(n,p,p_node_crpt,p_edge_crpt,sigma_in,sigma_out,crpt_type);

Ind = model_out.Ind; % matrix of edge indices (m by 2)
thetaij = model_out.thetaij; % given corrupted and noisy group ratios
ErrVec = model_out.ErrVec; % ground truth corruption levels
theta_orig = model_out.theta_orig; % ground truth absolute group element
thetaij_orig = model_out.thetaij_orig; % ground truth group ratios
% set CEMP defult parameters
parameters.beta_init = 1;
parameters.beta_max = 40;
parameters.rate = 1.2;


% run Spectral
theta_SP = Spectral_SO2(Ind, thetaij);

% run SDP  
theta_SDP = SDP_SO2(Ind, thetaij); %commented this out if cvx is not installed

% run IRLS
niter = 100;
theta_IRLS = IRLS_SO2(Ind, thetaij, niter);

% run CEMP+MST
theta_CEMP_MST = CEMP_MST_SO2(Ind, thetaij, parameters);

% run CEMP+GCW
theta_CEMP_GCW = CEMP_GCW_SO2(Ind, thetaij, parameters);

[mean_error_SP] = evaluate_error_SO2(theta_SP, thetaij_orig, Ind);
[mean_error_SDP] = evaluate_error_SO2(theta_SDP, thetaij_orig, Ind);
[mean_error_IRLS] = evaluate_error_SO2(theta_IRLS, thetaij_orig, Ind);
[mean_error_CEMP_MST] = evaluate_error_SO2(theta_CEMP_MST, thetaij_orig, Ind);
[mean_error_CEMP_GCW] = evaluate_error_SO2(theta_CEMP_GCW, thetaij_orig, Ind);


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

