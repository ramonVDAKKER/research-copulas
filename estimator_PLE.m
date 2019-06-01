function hat_theta=estimator_PLE(Data,model)
% hat_theta=estimator_PLE(Data,model)
% 
% This function computes the Pseudo Likelihood Estimator
% INPUTS:
% Data: n x p matrix of observations
% model: structure describing the model
% OUTPUT:
% hat_theta: k x 1 vector containing estimate of theta

%%
pseudoU=data2rank(Data); % compute marginal ranks (so estimator is invariant with respect to margins
%%
hat_theta=estimator_IML(pseudoU,model); % compute PLE by computing parametric MLE using rank data

