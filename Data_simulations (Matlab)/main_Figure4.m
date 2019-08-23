clc; close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load simulation parameters and initizalization;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global g01 g02 a  L Nt Nx
global X0 yG yR yn0 ypL yi yvac c0 kD kD_0 yV  ySp ySn yV0 
global yn_init1 yp_init1 yphi_init1 yvac_init1 


physical_param(); scaling_param(); mesh_param(); var_init()

load('physical_param.mat'); load('scaling_param.mat');
load('mesh&time_param.mat'); load('var_init.mat')

imax = 60; %Number of ion concentrations

g01 = 1; %Illumination from ETL 
g02 = 0; %No light from HTL

%Stored-variable initialization
PCE = 0; %Power conversion efficiency

%Initial guess on the max power point voltage (VMPP = Vapp).
Voc = VT*log(G/a*(1-exp(-a*L))/(L*B*Nc^2*exp(-Eg/VT))+1);
Vapp = Voc-0.1;
yV = Vapp/VT; %scaling
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of the max power point voltage (MPP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We solve the PDEs using the function "solve_PDE.m" and calculate the 
%guessed value of the current, voltage and PCE at max power point.
    
[ynx, ypx, yvacx, yPHI] = solve_PDE(xpos, T0, m, yvac_init, ynx_init, ypx_init, yPHI_init);

k = 1;
JMPPth(k) = current(ynx(Nt,:),ypx(Nt,:),xpos)*q*Dn*N0/X0*1e3;
VMPPth(k) = Vapp;
PCEth(k) = JMPPth(k).*VMPPth(k)
    
%We repeat the previous step but the new value of the guessed VMPP is
%Voc+0.01. In addition, we initialize the densities, current and 
%potential with the values calculated for the previous step (k=1).
    
Vapp = Vapp + 0.01;
yV = Vapp/VT;
yvac_init = yvacx; ynx_init = ynx; ypx_init = ypx; yPHI_init = yPHI;
    
[ynx, ypx, yvacx, yPHI] = solve_PDE(xpos, T0, m, yvac_init, ynx_init, ypx_init, yPHI_init);
 yvac_init = yvacx; ynx_init = ynx; ypx_init = ypx; yPHI_init = yPHI;
        
k = 2;
JMPPth(k) = current(ynx(Nt,:),ypx(Nt,:),xpos)*q*Dn*N0/X0*1e3;
VMPPth(k) = Vapp;
PCEth(k) = JMPPth(k).*VMPPth(k)

%We probe the maximum PCE value. If PCE(k=2) - PCE(k=1)>0, we keep
%increasing the guessed value of VMPP by +0.01 and repeat the previous 
%steps until the value of the PCE decreases. If PCE(k=2) - PCE(k=1)<0, 
%we keep decreasing the guessed value of VMPP by -0.01 and repeat the 
%previous steps until the value of the PCE decreases. 

if(PCEth(2)-PCEth(1) > 0)
    
    DeltaPCE = PCEth(2)-PCEth(1);
    
    while(DeltaPCE > 0)
            
            Vapp = Vapp + 0.01;
            yV = Vapp/VT;
            [ynx, ypx, yvacx, yPHI] = solve_PDE(xpos, T0, m, yvac_init, ynx_init, ypx_init, yPHI_init);
            yvac_init = yvacx; ynx_init = ynx; ypx_init = ypx; yPHI_init = yPHI;
            
            k = k+1;
            JMPPth(k) = current(ynx(Nt,:),ypx(Nt,:),xpos)*q*Dn*N0/X0*1e3;
            VMPPth(k) = Vapp;
            PCEth(k) = JMPPth(k).*VMPPth(k)
            DeltaPCE = PCEth(k) - PCEth(k-1);
    
    end
else
    
    DeltaPCE = PCEth(1)-PCEth(2);
    
    while(DeltaPCE > 0)
            
            Vapp = Vapp - 0.01
            yV = Vapp/VT;
            [ynx, ypx, yvacx, yPHI] = solve_PDE(xpos, T0, m, yvac_init, ynx_init, ypx_init, yPHI_init);
            yvac_init = yvacx; ynx_init = ynx; ypx_init = ypx; yPHI_init = yPHI;
            
            k = k+1;
            JMPPth(k) = current(ynx(Nt,:),ypx(Nt,:),xpos)*q*Dn*N0/X0*1e3;
            VMPPth(k) = Vapp;
            PCEth(k) = JMPPth(k).*VMPPth(k)
            DeltaPCE = PCEth(k) - PCEth(k-1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of the band diagram at MPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VMPP = VMPPth(k-1)
PCE = PCEth(k-1)
Vapp = VMPP;
yV = Vapp/VT;

[ynx, ypx, yvacx, yPHI] = solve_PDE(xpos, T0, m, yvac_init, ynx_init, ypx_init, yPHI_init);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot and save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos = pos*1E7;
Ec=-yPHI(Nt,:)*VT+Eg;
Ev=-yPHI(Nt,:)*VT;
EFp=-VT*log(ypx(Nt,:)*N0/Nc)+Ev;
EFn= VT*log(ynx(Nt,:)*N0/Nc)+Ec;

plot(pos, Ec,'k', pos, Ev, 'k', pos, EFn, '--', pos, EFp, '--')
axis([0 max(pos) 0 max(Ec)])

save BD_S_1E5_1E17.txt pos Ec Ev Efn Efp -ascii

