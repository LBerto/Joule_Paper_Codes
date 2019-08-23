function [current] = current(ynx,ypx,xpos)

global q yi yR ySn ySp X0 L Dn N0 a yG g01 g02 Nx


%Extrapolation of the electron and hole densities for integration.
un_x = pchip(xpos,ynx);
un_tx = @(x)ppval(un_x,x);
up_x = pchip(xpos,ypx);
up_tx = @(x)ppval(up_x,x);
 
%Calculation of the photogeneration and bulk recombination current 
yjbulk=@(x)g01*yG.*exp(-a*X0.*x) + g02*yG.*exp(-a*X0.*(L/X0-x))...
    -yR*(un_tx(x).* up_tx(x)-yi^2);

%Calculation of the total current.   
yjbulk = integral(yjbulk,0,L/X0) ;
yjSn = ySn*ynx(Nx);
yjSp = ySp*ypx(1);
current = yjbulk - yjSn - yjSp; %scaled current





