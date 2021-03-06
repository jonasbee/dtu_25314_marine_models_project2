%% Advection diffusion model
clear all
close all
clc

"Running"

%% 1) Parameters and initial conditions
%Call function: "call_param"
p=call_param();

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

% 4) Set time step
tt=3*365; %We run for 3 years

t1=[0:tt];

%4) Run the model
% Seasonal
[t,y]=ode45(@func_diff_season,t1,y,[],p);


%Splits the output into PP, N and D:
Ps=y(:,1:p.n);
Ns=y(:,p.n+1:2*p.n);
Ds=y(:,2*p.n+1:end);

%% CALCULATE LIGHT

%Calculate light
%Season
I=func_light_s(z,t,Ps(end,:),p);

%% Growth rate, Mortality and fitness
g1=p.eps*p.b*Ps./(p.b*Ps+p.Cmax).*p.Cmax-p.M;
m1=p.kl.*I+p.m;
w1=g1./m1';

%% Contour plots of m, g, w, and light
figure()
contourf(t,-z,m1)
c = colorbar;
c.Label.String = 'Mortality [1/day]';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Mortality")
grid on

figure()
contourf(t,-z,g1')
c = colorbar;
c.Label.String = 'Growth [1/day]';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Growth")
grid on

figure()
contourf(t,-z,w1')
c = colorbar;
c.Label.String = 'Fitness proxy';
ylabel("Depth [m]")
%xlabel("Time [days]")
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Fitness proxy")
grid on


%% AGENT BASED MODEL

% amount of agents (copepods)
% Choose a number :)
Copepods=1500


Z=ones(1,Copepods);
S=ones(1,Copepods);
R=zeros(1,Copepods);

for j=1:Copepods % Loop for each copepod

Z(1,j)=(-100+100*rand());

for i=1:tt % loop for each time step
%     
pos(i,j)=round(-Z(i,j));

% Making sure positions stays with frame
if pos(i,j)==0
    pos(i,j)=1;
end 

Dpos(i,j)=pos(i,j)+1;

if Dpos(i,j)==101
    Dpos(i,j)=100;
end

% Fitness proxy
g=p.b*Ps(i,pos(i,j))./(p.b*Ps(i,pos(i,j))+p.Cmax).*p.eps*p.Cmax-p.M;
m=p.kl.*I(pos(i,j))+p.m0;
w=g./m';

% 
gd=p.b*Ps(i,Dpos(i,j))./(p.b*Ps(i,Dpos(i,j))+p.Cmax).*p.eps*p.Cmax-p.M;
md=p.kl.*I(Dpos(i,j))+p.m0;
wd=gd./md';

dwdz(i)=(wd-w)/p.dz;

% Survival chance
S(i+1,j)=S(i,j)-m*S(i,j)*p.DeltaT;
if S(i+1,j)<0.001
    S(i+1,j)=0;
end

% Reproduction output
R(i+1,j)=R(i,j)+S(i,j)*g*p.DeltaT;

if R(i+1,j)<0.001
    R(i+1,j)=0;
end
    
% Movements    
Z(i+1,j)=Z(i,j)+(-1+2*rand())*(2*p.r^(-1)*p.Dc*p.DeltaT)^(1/2)-p.c*dwdz(i)*p.DeltaT;
if Z(i+1,j)>-1
    Z(i+1,j)=-1;
end
if Z(i+1,j)<-99
    Z(i+1,j)=-99;
end

end
end
%% Plot of movement on fitness contour
%5) plot Plankton with copepods movement
figure()
contourf(t,-z,w1')
c = colorbar;
c.Label.String = 'Fitness';
ylabel("Depth [m]")
hold on
plot(t,Z','Linewidth',1);
hold off
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Adaptive movement of 100 copepods")
grid on

%% Survival for all copepods
figure()
plot(t,S,'Linewidth',2)
ylabel("Survival change[%]")
xlim([0,100]) 
xlabel("Time [days]")
title("Survival")
grid on

%% Reproductive outout for all copepods
figure()
plot(t,R,'Linewidth',2)
ylabel("Reproduction output [mass eggs/mass copepod]")
xlim([0,100])
xlabel("Time [days]")
title("Fitness proxy = Mass specific reproduction output")
grid on
%% The average reproductive output 
"average fitness of copepods"
sum(R(end,:))/Copepods

% Script done
"done"
