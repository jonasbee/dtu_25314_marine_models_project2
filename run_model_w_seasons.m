%% Advection diffusion model
clear all
close all
clc


%1) Parameters:
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

%2) the grid dz=z(depth)/n ->  z=0.5*dz
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

%%
%4) Ode45
t1=[0:1/24:3*360];
[t,y]=ode45(@func_diff_season,t1,y,[],p);

Ps=y(:,1:p.n);
Ns=y(:,p.n+1:2*p.n);
Ds=y(:,2*p.n+1:end);


 %%
I=func_light_s(z,t,Ps(end,:),p);
I_0=func_light_s(z,t,0*Ps(end,:),p);
% 
% %% light test
% I=func_light_s(z,t,Ps(end,:),p);
% 
% figure() %plot
% plot(t,I,'Linewidth',2)
% ylabel("Depth [m]")
% xlabel("time [days]")
% title("Ligth intensity")
% %legend("Light intensity with plankton","Light intensity without plankton")
% grid on


%%
%{
figure()
subplot(1,3,1)

plot(Ns(end,:),-z,'Linewidth',2)
ylabel("Depth [m]")
xlabel("concentration [mmol/m^3]")
title("End distribution of Nutrients")
grid on

subplot(1,3,2)

plot(Ps(end,:),-z,'Linewidth',2)
ylabel("Depth [m]")
xlabel("concentration [cell/m^3]")
title("End distribution of Plankton")
grid on

subplot(1,3,3)

plot(Ds(end,:),-z,'Linewidth',2)
ylabel("Depth [m]")
xlabel("concentration [mmol/m^3]")
title("End distribution of detritus")
grid on
%}
%%

%5) plot Plankton
figure() %plot
contourf(t,-z,Ps')
c=colorbar;
c.Label.String = 'Concentration of PP [cells/m^3]';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
%xticklabels({'Mar','Jun','Sep','Dec','Mar','Jun','Sep','Dec','Mar','Jun','Sep','Dec'})
xticklabels({'Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Distribution of PP with seasonal forcing")
grid on

%%
%5) plot Nutrients
figure() %plot
contourf(t,-z,Ns')
c = colorbar;
c.Label.String = 'Concentration of N [mmol N/m^3]';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Distribution of N with seasonal forcing")
grid on

%%
%5) plot Detritus
figure() %plot
contourf(t,-z,Ds')
c = colorbar;
c.Label.String = 'Concentration of D [mmol N/m^3]';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Distribution of D with seasonal forcing")
grid on
%%
%5) plot light
figure() %plot
contourf(t,-z,I)
c = colorbar;
c.Label.String = 'Light intensity [mumol/(m^2s^1)]';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Seasonal forcing on the light intensity")
grid on