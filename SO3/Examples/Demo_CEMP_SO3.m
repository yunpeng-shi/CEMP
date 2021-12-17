
model_out = Uniform_Topology_SO3(200,0.5,0.3,0,'uniform');

Ind = model_out.Ind; % matrix of edge indices (m by 2)
RijMat = model_out.RijMat; % given corrupted and noisy relative rotations
ErrVec = model_out.ErrVec; % ground truth corruption levels
R_orig = model_out.R_orig; % ground truth rotations

% set CEMP defult parameters
parameters.beta_init = 1;
parameters.beta_max = 40;
parameters.rate = 1.2;

% run CEMP
SVec = CEMP_SO3(Ind,RijMat,parameters);

%visualize sij^* and sij,t, ideally one should see a straight line "y=x"
plot(ErrVec,SVec,'b.');
title('Scatter Plot of s_{ij}^* v.s. s_{ij,T}')
xlabel('s_{ij}^*') 
ylabel('s_{ij,T}') 
