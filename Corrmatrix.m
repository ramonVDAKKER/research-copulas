function [R,S,dotR,dotS]=Corrmatrix(theta,model) 
% [R,S,dotS]=Corrmatrix(theta,model)
%
% INPUTS:
% theta: k x 1 vector of parameters
% model: structure describing model  [see Gaussian_MC.m for description of
% available choices]
% OUTPUT:
% R: p x p correlation matrix evaluated at \theta
% S: inverse of R
% dotS: p x p x k matrix of derivatives with respect to \theta_m
% NOTE:
% This function does not check that the input theta indeed yields a
% correlation matrix, i.e. a positive definite matrix

%%
p=model.p;  % dimension copula
k=model.k;  % dimension parameter
type=model.type; % model type
%% Determine correlation matrix
R=zeros(p,p);
dotR=zeros(p,p,k);
if strcmp(type,'unrestricted')
    counter=1;
    for i=1:p-1
        for j=i+1:p
            R(i,j)=theta(counter,1);
            dotR(i,j,counter)=1; dotR(j,i,counter)=1; 
            counter=counter+1;
        end
    end
     R=R+R'+eye(p);
elseif strcmp(type,'exchangeable')
      R=theta*ones(p,p);
      R(1:p+1:p^2)=ones(p,1);
      dotR=ones(p,p)-eye(p);
elseif strcmp(type,'circular')
    if abs(p-4)>0
        err('The circular model is only available for p=4');
    end
    R=toeplitz([1 theta theta^2 theta]);
    dotR=toeplitz([0 1 2*theta 1]);
elseif strcmp(type,'AR1')
    aux_v=zeros(1,p);
    aux_v2=zeros(1,p);
    aux_v(1,1)=1;
    for j=2:p
      aux_v(1,j)=theta^(j-1);
      aux_v2(1,j)=(j-1)*theta^(j-2);
    end
    R=toeplitz(aux_v);
    dotR=toeplitz(aux_v2);
elseif strcmp(type,'MA')
   theta_ext=[1 ;theta; zeros(p-k-1,1)];
   R=toeplitz(theta_ext);
   for m=1:k
       aux_v=zeros(1,p);
       aux_v(1,m+1)=1;
       dotR(:,:,m)=toeplitz(aux_v); 
   end
elseif strcmp(type,'Toeplitz')
    R=toeplitz([1 theta']);
    for m=1:k  
       aux_v=zeros(1,p);
       aux_v(1,m+1)=1;
       dotR(:,:,m)=toeplitz(aux_v);
    end
elseif strcmp(type,'factor')
    aux=[1 ;theta];
    R=aux*aux';  R(1:p+1:p^2)=ones(p,1);
    for m=1:k
        aux_v=zeros(p,1);
        aux_v(m+1,1)=1;
       dotR(:,:,m)=aux_v*aux'+aux*aux_v'; dotR(1:p+1:p^2)=zeros(p,1); 
    end
else
    error('This model is not yet available')
end
%%
S=inv(R);  % inverse of R
dotS=zeros(p,p,k);
for j=1:k
    dotS(:,:,j)=-S*dotR(:,:,j)*S;  % p x p x k matrix containing derivatives of S with respect to \theta_m
end


