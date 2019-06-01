function [I_P,I_SP,I_PLE]=information_matrices(theta,model)
% [I_P,I_SP,I_PLE]=information_matrices(theta,model)
% 
% This function calculates, for a given model and at the parameter value theta, the parametric Fisher information matrix (I_P), the
% semiparametric efficient information matrix (I_SP), and the asymptotic variance
% matrix of the PLE (I_PLE)  [ which all are k x k matrices, where k
% denotes the dimension of the parameter space ]


%%
[R,S,dotR,dotS]=Corrmatrix(theta,model); % determine correlation matrix, its inverse and relevant derivatives
k=model.k; % dimension parameter space

%%  (parametric) Fisher information matrix
I_P=zeros(k,k); % initialization
for m=1:k
    aux1=-.5*dotS(:,:,m);
    for ell=1:k
        aux2=-.5*dotS(:,:,ell);
        I_P(m,ell)=2*trace(aux1*R*aux2*R);
    end
end
%% semiparamatric efficient information matrix
I_SP=zeros(k,k); % initialization
for m=1:k
    aux1=EfficientScore_QuadraticForm(R,S,dotR(:,:,m),dotS(:,:,m));
    for ell=1:k
        aux2=EfficientScore_QuadraticForm(R,S,dotR(:,:,ell),dotS(:,:,ell));
        I_SP(m,ell)=2*trace(aux1*R*aux2*R);
    end
end
%% PLE asymptotic variance matrix
PLE_aux=zeros(k,k); % initialization
for m=1:k
    aux_m=.5*diag(diag(dotR(:,:,m)*S));
    for ell=1:k
        aux_ell=.5*diag(diag(dotR(:,:,ell)*S));
        PLE_aux(m,ell)=2*trace(aux_m*R*aux_ell*R);
    end
end
asyvar_PLE=inv(I_P)+inv(I_P)*PLE_aux*inv(I_P); % asymptotic variance PLE
I_PLE=inv(asyvar_PLE);
