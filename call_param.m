function p=call_param()

%% Parameters for project 1

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

%% Parameters for project 2

p.r=1/3; %Randomization parameter [unitless]
p.c=10; %coefficient expressing the magnitude of 
%impact the fitness gradient has on movement [unitless]
p.b=0.1e-8; %Encounter rate [1/day]
p.Cmax=1.1; %max comsumption [(cell/m^3)/day]
p.eps=0.9; %consumption efficiency [1/(cell/m^3)]
p.m0=0.1; %Baseline mortality [1/day]
p.kl=0.004; %predation pr light pr day [1/((mumol photons/(m^2s))*day)]
p.Dc=10 %diffusion rate og copepods [m^2/day] 


p.M=p.Cmax*0.1; %Matabolic rate [1/day]

p.DeltaT=1; %time step [day]


p = p;
end 