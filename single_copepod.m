%% Advection diffusion model
clear all
close all
clc

"Running"

%% 1) Parameters and initial condidtions
%Call function: "call_param"
p=call_param()

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


%% RUN ODE for NPD model
% 4) Set time step !!!

tt=3*365

t1=[0:tt];

%4) Run the model
% without seasons
[t,y]=ode45(@func_diff,t1,y,[],p);

% Seasonal
%[t,y]=ode45(@func_diff_season,t1,y,[],p);

%Splits the outpuy into PP, N and D:
Ps=y(:,1:p.n);
Ns=y(:,p.n+1:2*p.n);
Ds=y(:,2*p.n+1:end);

%% CALCULATE LIGHT

%Calculate light
I=func_light(z,Ps(end,:),p);
%Season
%I=func_light_s(z,t,Ps(end,:),p);


%% Growth rate, Mortality and fitness
g0=p.eps*p.b*Ps(end,:)./(p.b*Ps(end,:)+p.Cmax).*p.Cmax-p.M;
m0=p.kl.*I+p.m0;
w0=g0./m0';

% Plot for section 6.2
%comparing mortality, growth and fitness
figure()
plot(w0,-z,'Color','#0072BD','Linewidth',2)
hold on
plot(g0,-z,'Color','#77AC30','Linewidth',2)
plot(m0,-z,'Color','#A2142F','Linewidth',2)
hold off
ylabel("Depth")
xlabel("[1/day],[unitless]")
title("Fitness proxy based on growth rate and mortality")
legend("Fitness [unitless]","Growth rate [1/day]","Mortality [1/day]")
grid on

%% AGENT BASED MODEL
% Copepod

Z(1)=(-100+100*rand());
S(1)=1
R(1)=0
for i=1:tt %Loop for each time step
%     
pos(i)=round(-Z(i));
Dpos(i)=pos(i)+1;

% Fitness proxy
g0=p.b*Ps(i,pos(i))./(p.b*Ps(i,pos(i))+p.Cmax).*p.eps*p.Cmax-p.M;
m0=p.kl.*I(pos(i))+p.m0;
w0=g0./m0';

% 
g0d=p.b*Ps(i,Dpos(i))./(p.b*Ps(i,Dpos(i))+p.Cmax).*p.eps*p.Cmax-p.M;
m0d=p.kl.*I(Dpos(i))+p.m0;
w0d=g0d./m0d';

dwdz(i)=(w0d-w0)/p.DeltaT;

% Survival chance
S(i+1)=S(i)-m0*S(i)*p.DeltaT;
if S(i+1)<0.001
    S(i+1)=0;
end

% Reproduction output
R(i+1)=R(i)+S(i)*g0*p.DeltaT;

if R(i+1)<0.001
    R(i+1)=0;
end

% Movements
Z(i+1)=Z(i)+(-1+2*rand())*(2*p.r^(-1)*p.Dc*p.DeltaT)^(1/2)-p.c*dwdz(i)*p.DeltaT;
if Z(i+1)>-1;
    Z(i+1)=-1;
end
if Z(i+1)<-99;
    Z(i+1)=-99;
end

end


%5) Plot of movement on top of Plankton concentrations
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

%% Script done
"done"
