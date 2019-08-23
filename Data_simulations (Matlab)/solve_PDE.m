function [ynx1, ypx1, yvac1, yPHI1] = solve_PDE(xpos, T, m, yvac0, ynx0, ypx0, yPHI0)

global Nt yvac_init1 yn_init1 yp_init1 yphi_init1

yvac_x0  = pchip(xpos,yvac0(Nt,:));
yvac_init1 = @(x)ppval(yvac_x0 ,x);
yn_x0  = pchip(xpos,ynx0(Nt,:));
yn_init1 = @(x)ppval(yn_x0 ,x);
yp_x0  = pchip(xpos,ypx0(Nt,:));
yp_init1 = @(x)ppval(yp_x0 ,x);
yphi_x0  = pchip(xpos,yPHI0(Nt,:));
yphi_init1 = @(x)ppval(yphi_x0 ,x);


sol1 = pdepe(m,@eqn,@initial,@bc,xpos,T);
ynx1 = sol1(:,:,1); ypx1 = sol1(:,:,2); yvac1 = sol1(:,:,3); yPHI1 = sol1(:,:,4);

save varinit.mat yvac_init1 yn_init1 yp_init1 yphi_init1


end

