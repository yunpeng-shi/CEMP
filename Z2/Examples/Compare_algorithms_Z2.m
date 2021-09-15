% parameters with uniform topology
%n =200; p=0.5; q=0.8; sigma=0; crpt_type='uniform';
% parameters with nonuniform topology
n =200; p=0.5; p_node_crpt=0.9; p_edge_crpt=0.75; crpt_type='adv';


% The following code is for nonuniform topology. This is a more malicious
% scenario where corrupted edges have cluster behavior (so local coruption
% level can be extremely high)
model_out = Nonuniform_Topology_Z2(n,p, p_node_crpt,p_edge_crpt, crpt_type);

Ind = model_out.Ind; % matrix of edge indices (m by 2)
zij = model_out.zij; % given corrupted and noisy group ratios
ErrVec = model_out.ErrVec; % ground truth corruption levels
z_orig = model_out.z_orig; % ground truth absolute group element
zij_orig = model_out.zij_orig; % ground truth group ratios
% set CEMP defult parameters
parameters.beta_init = 1;
parameters.beta_max = 40;
parameters.rate = 1.2;


% run Spectral
z_SP = Spectral_Z2(Ind, zij);

% run SDP
z_SDP = SDP_Z2(Ind, zij);

% run CEMP+GCW
z_CEMP_GCW = CEMP_GCW_Z2(Ind, zij, parameters);

[mean_error_SP] = evaluate_error_Z2(z_SP, zij_orig, Ind);
[mean_error_SDP] = evaluate_error_Z2(z_SDP, zij_orig, Ind);
[mean_error_CEMP_GCW] = evaluate_error_Z2(z_CEMP_GCW, zij_orig, Ind);


% Report estimation error
sz = [3 2];
varTypes = {'string','double'};
varNames = {'Algorithms','MeanError'};
Results = table('Size',sz,'VariableTypes',varTypes, 'VariableNames',varNames);
Results(1,:)={'Spectral', mean_error_SP};
Results(2,:)={'SDP', mean_error_SDP};
Results(3,:)={'CEMP+GCW', mean_error_CEMP_GCW};



Results

