function  A=EfficientScore_QuadraticForm(R,S,dotR,dotS)
% A=EfficientScore_QuadraticForm(R,S,dotR,dotS)
% 
% This auxiliary function is used to compute the matrix yielding the quadratic form for the
% efficient score
% INPUTS:
% R: p x p correlation matrix
% S: p x p inverse correlation matrix
% dotR: p x p matrix corresponding to derivative of R w.r.t. \theta_m
% dotS: p x p matrix corresponding to derivative of S w.r.t. \theta_m
% OUTPUT:
% p x p matrix A; interpretation x'*A*x is the m-th component of the efficient score evaluated at
% x=\Phi^{-1}(u)

%%
p=size(R,1); % dimension copula
%% determine g
g=-inv(eye(p)+(R.*S))*( (dotR.*S)*ones(p,1) );
%% determine matrix quadratic form  efficient score
A= -.5*dotS + 0.5*(diag(g)*S+S*diag(g));
    