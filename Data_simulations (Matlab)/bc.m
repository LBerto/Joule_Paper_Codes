function [pl,ql,pr,qr] = bc(xl,yl,xr,yr,t)
%BC1: MATLAB function M-file that specifies boundary conditions
%for a PDE in time and one space dimension.

global yV0 yV yn0 ypL ySp ySn 


pl = [yl(1)-yn0; -ySp*yl(2); 0; yl(4)];
ql = [0; 1; 1; 0];
pr = [ySn*yr(1); yr(2)-ypL; 0; yr(4)+yV0-yV];
qr = [1; 0; 1; 0];

