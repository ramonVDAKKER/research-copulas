function model_checker(model,theta)
% This function yields an error in case the model structure in combination
% with theta does not lead to a valid correlation matrix
 

%% Check length theta
if size(theta,2)>1
    error('theta should be a column vector.')
end
%%
if abs(model.k-length(theta))>0
    error('The dimension of theta should match the dimension of the parameter space.')
end
%% Modelspecific check on supplied values of p and k
if strcmp(model.type,'unrestricted')
    if abs(0.5*model.p*(model.p-1)-model.k)>0
       error('For the unrestricted model k=p(p-1)/2.')
    end
elseif strcmp(model.type,'circular')
    if abs(model.p-4)>0
        error('The circular model is only available for p=4.')
    end
elseif strcmp(model.type,'MA')
    if model.p<1+model.k
        error('For the MA model p should be >=1+k.')
    end
elseif strcmp(model.type,'Toeplitz')
   if model.p< 1+model.k
       error('For the Toeplitz model p should be >=1+k.')
   end
elseif strcmp(model.type,'factor')
   if abs(model.p-( 1+model.k))>0
       error('For the factor model p should be  equal to 1+k.')
   end
end
%% Check that supplied value of theta yields positive definite correlation matrix
[R,~,~,~]=Corrmatrix(theta,model);
if min(eig(R))<=0
    error('The input value of theta does not lead to a positive definite matrix R.')
end
