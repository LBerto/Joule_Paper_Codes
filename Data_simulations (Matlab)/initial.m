function value = initial(x)
%MATLAB function M-file that specifies the initial condition
%for a PDE in time and one space dimension.

global yn_init1 yp_init1 yphi_init1 yvac_init1 

value = [yn_init1(x); yp_init1(x); yvac_init1(x); yphi_init1(x)];

