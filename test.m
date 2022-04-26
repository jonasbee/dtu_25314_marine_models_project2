%% Advection diffusion model
clear all
close all
clc

"Running"

% hello git
%% 1) Parameters:
p.u= 0.96; %sinking speed [m/day]
p.D= 10;%diffusion rate [m^2/day] 
p.depth= 100; %depth [m]
p.n= 100; %amount of sections

p.m=0.24; %1/day
p.mu_max=0.96;%max speciffic growth rate [1/d]

p.I0=350; % incident light intensity [mumol photons m^-2 s^-1]
p.Kbg=0.2; %(0.045) %background turbidity [1/m]
p.k=15e-12; %(6e-10) %specific light attenuation of phytoplankton [m^2/cell]
p.Hi=30; %Half satuaration ligth growth [mumol photons/(m^2s)]

p.alpha=1*10^(-9); %nutrient content of pp [mmol N/cell]
p.N_b=50; %bottom con af nutrients [mmol/m^3]
p.Hn=0.04; %Half satuaration nutrients growth [mmol N/(m^3)]

p.tau=0.1; %1/d  remineralization
p.w=5; %sinking speed of detritus m/d

% New parameters
p.r=1/3
p.c=0.5
p.b=0.5
p.h=0.5
p.kl=0.5

p.DeltaT=1


%2) tThe grid
p.dz=p.depth/p.n; %width of seciton
p.z=0.5*p.dz:p.dz:(p.depth-0.5*p.dz); %The grid
z=0.5*p.dz:p.dz:(p.depth-0.5*p.dz); %The grid


%3) initial conditions Plankton
P0 = 2e9*exp(-(p.z-p.depth/4).^2/1000); %Gauss distribution

%3) initial condition for Nutrients
N0 = p.N_b*exp(-(p.z-p.depth/1.8).^2/500); %Gauss distribution

%3) initial condition for Ditritus
D0=zeros(1,p.n);

y=[P0,N0,D0];


%% RUN ODE
% 4) Set time step !!!

tt=3*365

t1=[0:tt];

%4) Run the model
[t,y]=ode45(@func_diff,t1,y,[],p);

%Splits the outpuy into PP, N and D:
Ps=y(:,1:p.n);
Ns=y(:,p.n+1:2*p.n);
Ds=y(:,2*p.n+1:end);

%% CALCULATE LIGHT

%Calculate light
I=func_light(z,Ps(end,:),p);

%% AGENT BASED MODEL
% Copepods

g0=p.b*Ps(t,z(t))/(1+p.b*p.h*z(t)
m0=p.kl*I(t+1,z(t))
w0=g0/m0
dwdz=
Z(1)=-30;

for i=1:tt
Z(i+1)=Z(i)+(-1+2*rand())*(2*p.r^(-1)*p.D*p.DeltaT)^(1/2)+p.c*dwdz*p.DeltaT;
if Z(i+1)>0;
    Z(i+1)=0;
end
if Z(i+1)<100;
    Z(i+1)=100;
end
end

%+p.c*dwdz*p.DeltaT

%g=p.b*Ps(t+1,z(t+1))/(1+p.b*p.h*z(t+1)

%m=p.kl*I(t+1,z(t+1))

%w=g/m 

%% 5) plot Plankton
figure() %plot
contourf(t,-z,Ps')
c=colorbar;
c.Label.String = 'Concentration of PP [cells/m^3]';
hold on
plot(t,Z','Color','#A2142F','Linewidth',1)
hold off
ylabel("Depth [m]")
xlabel("Time [days]")
title("Distribution of Phytoplankton")
grid on

%%

"done"
