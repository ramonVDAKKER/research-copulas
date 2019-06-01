function Output=Gaussian_MC(model,theta,n,reps)
% Output=Gaussian_MC(model,theta,n,reps) 
% model is structure describing the model
%   model.type:  availabe are 'unrestricted', 'exchangeable', 'circular',
%   'AR1', 'MA', 'Toeplitz', 'factor'
% theta determines DGP
% n is number of observations
% reps is number of replications

%%
model_checker(model,theta) % some checks to validate input
%% determine correlation matrix DGP
[R_0,~,~,~]=Corrmatrix(theta,model);
%%
for sim=1:reps  % replication loop
    U=copularnd('Gaussian',R_0,n);  % generate n iid draws from Gaussian copula
    %% compute estimators
    % infeasible parametric maximum likelihood
    theta_IMLE(:,sim)=estimator_IML(U,model);
    % pseudo MLE (data is transformed to rank data within function, so we can just U here)
    theta_PMLE(:,sim)=estimator_PLE(U,model);
    % One step estimator (data is transformed to rank data within function, so we can just U here)
    % PLE is used as initial estimate
    [theta_upd_I,theta_upd_EI]=estimator_OSE(U,theta_PMLE(:,sim),model);
    theta_OSE_I(:,sim)=theta_upd_I;
    theta_OSE_EI(:,sim)=theta_upd_EI;
end
%% organize estimation results in structure Output
Output.theta_IMLE=theta_IMLE;
Output.theta_PMLE=theta_PMLE;
Output.theta_OSE_I=theta_OSE_I;
Output.theta_OSE_EI=theta_OSE_EI;
end
