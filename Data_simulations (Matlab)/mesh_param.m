function [] = mesh_param()

load('physical_param.mat');
load('scaling_param.mat');

%Space and time meshing constants
m = 0; % For pdepe solver
Nx=1000; %Number of mesh points for the space meshing
Nt=101; %Number of mesh points for the time meshing
tmax0=1; %Time scale for the simulation (in seconds)

%Space and time meshing
pos=linspace(0,L,Nx); %Space mesh
xpos=linspace(0,L/X0,Nx); %Normalized Space mesh
t0 = linspace(0,tmax0,Nt); %Time mesh
T0 = linspace(0,tmax0/tau,Nt); %Normalized Time mesh

save mesh&time_param.mat m Nt Nx xpos pos t0 T0 

end