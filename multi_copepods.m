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

% Parameters for project 2
p.r=1/3; %Randomization parameter [unitless]
p.c=10; %coefficient expressing the magnitude of 
%impact the fitness gradient has on movement [unitless]
p.b=0.1e-8; %Encounter rate [1/day]
p.Cmax=1.1; %max comsumption [(cell/m^3)/day]
p.eps=0.9; %consumption efficiency [1/(cell/m^3)]
p.m0=0.2; %Baseline mortality [1/day]
p.kl=0.008; %predation pr light pr day [1/((mumol photons/(m^2s))*day)]


p.M=p.Cmax*0.1; %Matabolic rate [1/day]

p.DeltaT=1; %time step [day]


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

tt=3*365;

t1=[0:tt];

%4) Run the model
%[t,y]=ode45(@func_diff,t1,y,[],p);
% Seasonal
[t,y]=ode45(@func_diff_season,t1,y,[],p);



%Splits the outpuy into PP, N and D:
Ps=y(:,1:p.n);
Ns=y(:,p.n+1:2*p.n);
Ds=y(:,2*p.n+1:end);

%% CALCULATE LIGHT

%Calculate light
%I=func_light(z,Ps(end,:),p);
%Season
I=func_light_s(z,t,Ps(end,:),p);


%% Growth rate, Mortality and fitness
g1=p.b*Ps./(p.b*Ps+p.Cmax).*p.Cmax-p.M;
m1=p.kl.*I+p.m;
w1=g1./m1';

figure()
contourf(t,-z,w1','Color','#77AC30','Linewidth',2)
ylabel("depth")
xlabel("Fitness")
title("fitness")
grid on
%% Contour plots of m, g, w, and light
figure()
contourf(t,-z,m1)
c = colorbar;
c.Label.String = 'Mortality';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Mortality")
grid on

figure()
contourf(t,-z,g1')
c = colorbar;
c.Label.String = 'Growth';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Growth")
grid on

figure()
contourf(t,-z,w1')
c = colorbar;
c.Label.String = 'Fitness';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Fitness")
grid on

%5) plot light
figure() %plot
contourf(t,-z,I)
c = colorbar;
c.Label.String = 'Light intensity [mumol/(m^2s^1)]';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 0.5; % Make tick marks thicker.
title("Seasonal forcing on the light intensity")
grid on

%% AGENT BASED MODEL


% amount of agents (copepods)
Copepods=3


Z=ones(1,Copepods)
for j=1:length(Z)

Z(1,j)=(-100+100*rand())

for i=1:tt
%     
pos(i,j)=round(-Z(i,j));
Dpos(i,j)=pos(i,j)+1;

g0=p.b*Ps(i,pos(i,j))./(p.b*Ps(i,pos(i,j))+p.Cmax).*p.Cmax-p.M;
m0=p.kl.*I(pos(i,j))+p.m0;
w0=g0./m0';

% 
g0d=p.b*Ps(i,Dpos(i,j))./(p.b*Ps(i,Dpos(i,j))+p.Cmax).*p.Cmax-p.M;
m0d=p.kl.*I(Dpos(i,j))+p.m0;
w0d=g0d./m0d';

dwdz(i)=(w0d-w0)/p.DeltaT;

    
    
Z(i+1,j)=Z(i,j)+(-1+2*rand())*(2*p.r^(-1)*p.D*0.5*p.DeltaT)^(1/2)-p.c*dwdz(i)*p.DeltaT;
if Z(i+1,j)>-1;
    Z(i+1,j)=-1;
end
if Z(i+1,j)<-99;
    Z(i+1,j)=-99;
end

end
end

%5) plot Plankton with copepods movement
figure() %plot
contourf(t,-z,Ps')
c=colorbar;
c.Label.String = 'Concentration of PP [cells/m^3]';
hold on
plot(t,Z','Linewidth',0.7);
hold off
ylabel("Depth [m]")
xlabel("Time [days]")
title("Distribution of Phytoplankton")
grid on

%%
figure()
contourf(t,-z,w1')
c = colorbar;
c.Label.String = 'Fitness';
ylabel("Depth [m]")
hold on
plot(t,Z','Linewidth',1);
hold off
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Movement of 3 copepods")
grid on

"done"
