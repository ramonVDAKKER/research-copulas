function [theta_upd_I,theta_upd_EI]=estimator_OSE(Data,theta_init,model)
% theta_upd=estimator_OSE(Data,theta_init,model)
% 
% This function computes the One Step Estimator
% INPUTS:
% Data:  n x p matrix of observations
% theta_init: k x 1 vector containing  initial estimate of \theta [ should
% be an element of parameter space ]
% model: structure describing model
% OUTPUT:
% theta_upd_I: k x 1 vector containing one step estimate (using Fisher
% information matrix)
% theta_upd_EI: k x 1 vector containing one step estimate (using outer-poduct of scores 
% as estimate of Fisher % information matrix)

%%
Data=data2rank(Data);  % transform data to marginal ranks
[n,p]=size(Data); % number of observations
k=length(theta_init); % dimension parameter
%% calculate correlation matrix, inverse and its derivatives at initial estimate of theta
[R,S,dotR,dotS]=Corrmatrix(theta_init,model);
%% compute efficient information matrix
[~,info_matrix,~]=information_matrices(theta_init,model);
%% compute efficient scores and  compute estimate of efficient information matrix using outer product of efficient scores
eff_score_matrix=zeros(k,n);
%cum_eff_score=zeros(k,1);
X=norminv(Data'); % calculate Gaussianized data(p x n matrix)
%A=zeros(p,p,k);% 3-d matrix containing matrices yielding quadratic forms
for i=1:n
    X_i=X(:,i);
    for m=1:k
        A=EfficientScore_QuadraticForm(R,S,dotR(:,:,m),dotS(:,:,m));
        eff_score_matrix(m,i)=X_i'*A*X_i;
    end
end
cum_eff_score=sum(eff_score_matrix,2);
estimate_info_matrix=cov(eff_score_matrix');
estimate_info_matrix=estimate_info_matrix+(mean(eff_score_matrix'))'*mean(eff_score_matrix') ;

%% compute updates
theta_upd_I=theta_init+inv(info_matrix)*cum_eff_score/n;
theta_upd_EI=theta_init+inv(estimate_info_matrix)*cum_eff_score/n;





