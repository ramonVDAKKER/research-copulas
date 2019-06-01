function hat_theta=estimator_IML(Data,model)
% hat_theta=estimator_IML(Data,model) 
%
% This function calculates the (infeasible) parametric MLE (for U[0,1]
% margins)
% INPUTS:
% Data:   n x p matrix with observations from p-dimensional Gaussian copula
% model:  structure describing the specific model 
% OUTPUT:
% hat_theta: k x 1 vector containing estimate 
% NOTE:
% the MATLAB routine "fmincon" is used to maximize the likelihood; in
% the application of this algorithm  the likelihood is set equal to minus infinity if  the supplied "correlation
% matrix" is not positive definite.


%%
k=model.k;  % dimension parameter
p=model.p;  % dimension copula
%% Settings for fmincon
if strcmp(model.type,'exchangeable')
    lb=-ones(k,1)/(p-1); % for exchangeable model we know that \theta>-1/(p-1) is the exact requirement for positive definite correlation matrix an we thus pass this as information to fmincon
else
    lb=-ones(k,1); % for all other models the parameters should  be >=-1 (necessary, not sufficient [non positive definite matrices can still pop up])
end
options = optimset('GradObj','on','Algorithm','interior-point' );
x0=zeros(k,1); % for all currently implemented models x0 is a point of the parameter space and can thus be used as starting value of the algorithm
ub=ones(k,1); % for all implemented models the parameters should be at most +1 
A=[]; b=[]; Aeq=[]; beq=[]; nonlcon=[]; % no other linear (in)equalities
%% calculate estimator by minimizing minus log likelihood
hat_theta=fmincon(@(theta)aux_mLL(Data,theta,model),x0,A,b,Aeq,beq,lb,ub,nonlcon,options); 
    
function [f,g]=aux_mLL(Data,theta,model)
% [f,g]=aux_mLL(Data,theta,model)
%
% auxiliary function that calculates minus the log-likelihood (f) and its
% gradient (g [ k x 1 vector ] )

%%
n=size(Data,1);
k=model.k;
[R,S,dotR,dotS]=Corrmatrix(theta,model);  % R is correlation-matrix calculated at \theta
%% calculate f (minus log likelihood)
if  min(eig(R))<=0 % if R is not positive definite likelihood is - infinity
    f=Inf;   
else
    f=-sum(log(copulapdf('Gaussian',Data,R)));
end
%% calculate gradient g [note g=-\ell_\theta ]
g=zeros(k,1);
X=(norminv(Data'));
for m=1:k
    dotRm=dotR(:,:,m);
    dotSm=dotS(:,:,m);
    g(m,1)=.5*(n*trace(dotRm*S)+   (X(:))'*(kron(eye(n),dotSm))    *(X(:)) ) ;
end

end

end

