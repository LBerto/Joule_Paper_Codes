function [] = var_init()

load('mesh&time_param.mat');
load('scaling_param.mat');

%Initialization of the variables
yvac_init = yvac*ones(Nt,Nx);  
ynx_init = yn0*ones(Nt,Nx); 
ypx_init = ypL*ones(Nt,Nx); 
yPHI_init = zeros(Nt,Nx);


save var_init.mat ynx_init ypx_init yPHI_init yvac_init

end
