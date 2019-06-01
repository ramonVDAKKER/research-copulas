%% AR(1) 
model.type='AR1';
model.k=1;
model.p=5;
reps=1000;
n=100;
theta=.5;
%%
O=Gaussian_MC(model,theta,n,reps);
  