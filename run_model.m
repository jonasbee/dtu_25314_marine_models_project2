%% Advection diffusion model
clear all
close all
clc

"Running"

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
%4) Run the model
t1=[0:3*365];
[t,y]=ode45(@func_diff,t1,y,[],p);

%Splits the outpuy into PP, N and D:
Ps=y(:,1:p.n);
Ns=y(:,p.n+1:2*p.n);
Ds=y(:,2*p.n+1:end);

%% CALCULATE LIGHT

%Calculate light with and with the effect of pp above
I=func_light(z,Ps(end,:),p);
I_0=func_light(z,0*Ps(end,:),p);

% %plot the two light conditions
% figure()
% hold on
% plot(I,-z,'g','Linewidth',2)
% plot(I_0,-z,'b','Linewidth',2)
% ylabel("Depth [m]")
% xlabel("I [mumol/(m^2s^1)]")
% hold off
% title("Ligth intensity")
% legend("With the effect of plankton above","Without the effect of plankton above")
% grid on


%% limiting factors
%Calculate the light limited growht and the nutrients limited growth
gL=I./(I+p.Hi);

Nend=Ns(end,:)';
gN=Nend./(Nend+p.Hn);

%Plot the limiting factors and highlight the minimum values (liebigs law)
figure() %plot
hold on
plot(gL,-z,'color','#EDB120','Linewidth',2.5)
plot(gN,-z,'color','#0072BD','Linewidth',2.5)
plot(gN(1:9),-z(1:9),'--','color','#77AC30','Linewidth',1.5)
plot(gL(9:100),-z(9:100),'--','color','#77AC30','Linewidth',1.5)
ylabel("Depth [m]")
xlabel("[unitless]")
hold off
title("Limiting factors")
legend("Light limited growth","Nutrient limited growth","Growth according to Liebig's law")
grid on

%% STEADY STATE
% Plot end distrubution of Nutrients, PP and Detritus as subplots:
figure()

subplot(1,3,1)

plot(Ns(end,:),-z,'Color','#0072BD','Linewidth',2)
ylabel("Depth [m]")
xlabel("Concentration N [mmol N/m^3]")
title("Steady State Solution of Nutrients")
grid on


subplot(1,3,2)

plot(Ps(end,:),-z,'Color','#77AC30','Linewidth',2)
ylabel("Depth [m]")
xlabel("Concentration PP [cell/m^3]")
title("Steady State Solution of Phytoplankton")
grid on

subplot(1,3,3)

plot(Ds(end,:),-z,'Color','#A2142F','Linewidth',2)
ylabel("Depth [m]")
xlabel("Concentration D [mmol N/m^3]")
title("Steady State Solution of Detritus")
grid on
%% Contur plots of N, PP and D

%5) plot Plankton
figure() %plot
contourf(t,-z,Ps')
c=colorbar;
c.Label.String = 'Concentration of PP [cells/m^3]';
ylabel("Depth [m]")
xlabel("Time [days]")
title("Distribution of Phytoplankton")
grid on


%%
%5) plot Nutrients
figure() %plot
contourf(t,-z,Ns')
c = colorbar;
c.Label.String = 'Concentration of N [mmol N/m^3]';
ylabel("Depth [m]")
xlabel("Time [days]")
title("Distribution of Nutrients")
grid on

%%
%5) plot Detritus
figure() %plot
contourf(t,-z,Ds')
c = colorbar;
c.Label.String = 'Concentration of D [mmol N/m^3]';
ylabel("Depth [m]")
xlabel("Time [days]")
title("Distribution of Detritus")
grid on

%% done :)
"Done :)"
