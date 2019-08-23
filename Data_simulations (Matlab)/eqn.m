
function [c,b,s] = eqn(x,t,y,DyDx)
%MATLAB function M-file that specifies
%a PDE in time and one space dimension.

global g01 g02 yG a X0 kD kD_0 yR yi yvac c0 L

%Generation term: if g01 = 1, illumination occurs from the ETL and/or if 
%g02 = 1, illumination occurs from the HTL.

yg = g01*yG*exp(-a*X0*x) + g02*yG*exp(-a*X0*(L/X0 - x));

%Bulk PDEs for electrons, holes, vacancies and the electrostatic potential.
c = [1; 1; 1; 0];
b = [(DyDx(1) - c0*y(1)*DyDx(4)); (DyDx(2) + c0*y(2)*DyDx(4))*kD; (DyDx(3) + c0*y(3)*DyDx(4))*kD_0; DyDx(4)];
s = [yg - yR*(y(1)*y(2) - yi^2); yg - yR*(y(1)*y(2) - yi^2); 0; y(2) - y(1) + (y(3) - yvac)];