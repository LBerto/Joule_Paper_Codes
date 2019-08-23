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

imax = 37; %Number of bulk recombination rates
jmax = 46; %Number of surface recombination rates 

g01 = 1; %Illumination from ETL 
g02 = 0; %No light from HTL

%Stored-variable initialization
krec = 0; %for bulk recombination rate
Ssurf = 0; %for surace recombination velocity
PCE = 0; %Power conversion efficiency


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of the solar cell efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for i = 1:imax
    
    %We vary the band to band recombination (B) rate at the ETL (Sn) by 
    %step of 10 for every decade from 10^-10 to 10^-7 cm^3.s^-1. For each 
    %value of B we calculate an estimate of Voc and and estimate of the max
    %power point voltage (VMPP_guess = Voc - 0.1)
    
    B = (i-9*floor((i-1)/9))*10^(-10+floor((i-1)/9)); 
    Voc = VT*log(G/a*(1-exp(-a*L))/(L*B*Nc^2*exp(-Eg/VT))+1);
    Vapp = Voc-0.1;
    yV = Vapp/VT; %scaling
    
    for j = 1:jmax
    
    %We vary the surface recombination rate at the ETL (Sn) and HTL (Sp) by
    %step of 10 for every decade from 1 to 10^5 cm.s^-1. 
    
    S = (j-9*floor((j-1)/9))*10^(floor((j-1)/9));
    Sn = S; Sp = S;
    ySp = X0/Dn*Sp; ySn = X0/Dn*Sn; %Scaling

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

    %The VMPP guess for the next value of the bulk/surface recombination 
    %rate is the final value of VMPP for the current value of the
    %bulk/surface recombination.
    
    Vapp = VMPPth(k-1); 
    yV = Vapp/VT;
    
    %Storage of the recombination rates and PCE
    PCE(i,j) = PCEth(k-1)
    Ssurf(j) = S
    krec(i) = B

    
    save PCE_SvsB_1E17.txt PCE -ascii
    save S.txt Ssurf  -ascii
    save B.txt krec  -ascii

    end
end
