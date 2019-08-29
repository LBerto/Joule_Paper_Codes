function [] = physical_param()

%Device thickness
L = 8E-5;

%Physical constants
q=1.6*10^-19; %absolute value of the charge of the electron (in C)
VT=0.026; %Thermal voltage at room temperature (in V)
e0=8.85*10^-14; %Dielectric permittivity of the vacuum (in F.cm^-1)
er=24; %Dielectric constant of the material

%Band gap and injection barriers
Eg=1.66; %Semiconductor bandgap (in eV)
bn0=0.23; %Electron injection barrier in x=0 (in eV)
bpL=0.23; %Hole injection barrier in x=L (in eV)
Vbi=Eg-bn0-bpL; %Built-in potential within the device (in V)

%Densities
Nvac=1*10^17;%Concentration of ion vacancies 
Nc=3*10^18; %Density of states at the bottom of the conduction band 
ni=Nc*exp(-Eg/(2*VT)); %Intrinsic equilibrium carrier concentration 
n0s=Nc*exp(-bn0/VT); %Surface concentration of electrons in x=0
pLs=Nc*exp(-bpL/VT); %Surface concentration of holes in x=L

%Optoelectronic parameters
a = 2.5*10^4; %absorption coefficient (in cm^-1)
G = a*1.485*10^17; %Generation rate (in cm^-3.s^-1)
B = 5*10^-10; %Band to band recombination factor (in cm^3.s^-1)
Sn = 10^5; %Electron surface recombination velocity at HTL (in cm.s^-1)
Sp = 10^5; %Hole surface recombination velocity at ETL (in cm.s^-1)
Dn = 0.8; %Electron diffusion coefficient (in cm^2.s^-1)
Dp = Dn; %Hole diffusion coefficient
kD=Dp/Dn; 
Dion = 7.8e-9; %Vacancy diffusion coefficient
kD_0 = Dion/Dn;

save physical_param.mat q e0 er VT Eg bn0 bpL Vbi Nvac Nc ni n0s pLs  ...  
                        a G B Sn Sp Dn Dion kD_0 kD L

