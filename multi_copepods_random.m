%% Advection diffusion model
clear all
close all
clc

"Running"

% hello git
%% 1) Parameters 
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


%% RUN ODE
% 4) Set time step !!!

tt=3*365;

t1=0:tt;

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
contourf(t,-z,w1')
ylabel("depth")
xlabel("Fitness")
title("fitness")
grid on
%% Contour plots of m, g, w, and light
% figure()
% contourf(t,-z,m1)
% c = colorbar;
% c.Label.String = 'Mortality';
% ylabel("Depth [m]")
% %xlabel("Time [days]")
% xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
% xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
% ax = gca;
% ax.TickLength = [0.02,0]; % Make tick marks longer.
% ax.LineWidth = 1; % Make tick marks thicker.
% title("Mortality")
% grid on
% 
% figure()
% contourf(t,-z,g1')
% c = colorbar;
% c.Label.String = 'Growth';
% ylabel("Depth [m]")
% %xlabel("Time [days]")
% xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
% xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
% ax = gca;
% ax.TickLength = [0.02,0]; % Make tick marks longer.
% ax.LineWidth = 1; % Make tick marks thicker.
% title("Growth")
% grid on
% 
% figure()
% contourf(t,-z,w1')
% c = colorbar;
% c.Label.String = 'Fitness';
% ylabel("Depth [m]")
% %xlabel("Time [days]")
% xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
% xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
% ax = gca;
% ax.TickLength = [0.02,0]; % Make tick marks longer.
% ax.LineWidth = 1; % Make tick marks thicker.
% title("Fitness")
% grid on
% 
% %5) plot light
% figure() %plot
% contourf(t,-z,I)
% c = colorbar;
% c.Label.String = 'Light intensity [mumol/(m^2s^1)]';
% ylabel("Depth [m]")
% %xlabel("Time [days]")
% xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
% xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter'})
% ax = gca;
% ax.TickLength = [0.02,0]; % Make tick marks longer.
% ax.LineWidth = 0.5; % Make tick marks thicker.
% title("Seasonal forcing on the light intensity")
% grid on

%% AGENT BASED MODEL
% copepod's fitness based on random walk

% amount of agents (copepods)
Copepods=1500

Z=ones(1,Copepods);
S=ones(1,Copepods);
R=zeros(1,Copepods);

for j=1:Copepods

Z(1,j)=(-100+100*rand());

for i=1:tt
%     
pos(i,j)=round(-Z(i,j));

if pos(i,j)==0
    pos(i,j)=1;
end 

Dpos(i,j)=pos(i,j)+1;

if Dpos(i,j)==101
    Dpos(i,j)=100;
end

g=p.b*Ps(i,pos(i,j))./(p.b*Ps(i,pos(i,j))+p.Cmax).*p.Cmax-p.M;
m=p.kl.*I(pos(i,j))+p.m0;
w=g./m';

% 
gd=p.b*Ps(i,Dpos(i,j))./(p.b*Ps(i,Dpos(i,j))+p.Cmax).*p.Cmax-p.M;
md=p.kl.*I(Dpos(i,j))+p.m0;
wd=gd./md';

% dwdz(i)=(wd-w)/p.DeltaT;

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
Z(i+1,j)=Z(i,j)+(-1+2*rand())*(2*p.r^(-1)*p.D*0.5*p.DeltaT)^(1/2);%-p.c*dwdz(i)*p.DeltaT;
if Z(i+1,j)>-1;
    Z(i+1,j)=-1;
end
if Z(i+1,j)<-99;
    Z(i+1,j)=-99;
end



end
end
%%
%5) plot Plankton with copepods movement
figure() %plot
contourf(t,-z,Ps')
c=colorbar;
c.Label.String = 'Concentration of PP [cells/m^3]';
hold on
plot(t,Z','Linewidth',0.7);
hold off
xticks([0 365/4 365/2 365*3/4 365 365+365/4 365+365/2 365+365*3/4 365*2 2*365+365/4 2*365+365/2 2*365+365*3/4 365*3])
xticklabels({'Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring','Summer','Fall','Winter','Spring'})
ax = gca;
ax.TickLength = [0.02,0]; % Make tick marks longer.
ax.LineWidth = 1; % Make tick marks thicker.
title("Distribution of Phytoplankton and Copepod movement")
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

%%
figure()
plot(t,S,'Linewidth',2)
ylabel("Survival change[%]")
xlim([0,100]) 
xlabel("Time [days]")
title("Survival")
%legend("Fitness[unitless]","Growth rate[1/day]","Mortality[1/day]")
grid on

%%
figure()
plot(t,R,'Linewidth',2)
ylabel("Reproduction output [mass eggs/mass copepod]")
xlim([0,100])
xlabel("Time [days]")
title("Fitness proxy = Mass specific reproduction output")
%legend("Fitness[unitless]","Growth rate[1/day]","Mortality[1/day]")
grid on
%%
"average fitness of copepods"
sum(R(end,:))/Copepods

"done"
