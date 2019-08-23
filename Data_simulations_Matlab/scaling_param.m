function [] = scaling_param()

load('physical_param.mat');

%Normalization/scaling factors.
G0=1*10^30; %Generation rate normalization factor
V0 = VT; %Potential normalization factor
c0 = V0/VT; 
N0=sqrt(e0*er*V0*G0/(q*Dn)); % Density normalization factor
X0=sqrt(e0*er*V0/(q*N0)); %Length normalization factor
tau=X0^2/Dn;% Time normalization factor
j0 = q*Dn*N0/X0; %Current normalization factor

%Normalization/scaling of the generation/recombination rates and densities.
yG=G/G0; %Normalized generation rate 
yR=e0*er*V0*B/(q*Dn); %Normalized recombination rate 
ySp = X0/Dn*Sp;  %Normalize hole surface recombination velocity
ySn = X0/Dn*Sn; %Normalize electron surface recombination velocity
yi=ni/N0;  %Normalized intrinsic free carrier density
yvac=Nvac/N0; %Normalized vacancy density
yn0=n0s/N0; %Normalized free electron density at the left contact (x=0)
ypL=pLs/N0; %Normalized free hole density at the right contact (x=L)
yV0=Vbi/V0; %Normalized built-in potential


save scaling_param.mat V0 c0 N0 X0 tau j0 yG yR ySn ySp yi yn0 ypL yvac yV0 

end