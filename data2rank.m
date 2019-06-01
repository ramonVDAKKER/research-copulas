function  r=data2rank(Data)  
% r=data2rank(Data) 
%
% INPUT:
% Data: n x p matrix (n is number of observations) 
% OUTPUT:
% r: is n x p matrix whose columns contain ``marginal ranks'' scaled by
% 1/(n+1)


%%
[n,p]=size(Data);
r=zeros(n,p); 
for j=1:p
   r(:,j)=tiedrank(Data(:,j))/(n+1);  % compute marginal EDFs evaluated at ranks / (n+1)
end