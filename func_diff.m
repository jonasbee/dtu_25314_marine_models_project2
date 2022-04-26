function dydt=flux_vec(t,y,p)
%splitting y into PP,N and D
P=y(1:p.n);
N=y(p.n+1:2*p.n);
D=y(2*p.n+1:end);

%% Plankton
% boundary condition #NO FLUX
Ja(1)=0; %surface
Jd(1)=0; %surface
Ja(p.n+1)=0; %bottom
Jd(p.n+1)=0; %bottom

in = 2:(p.n); %intorier
Ja(in) = p.u*P(in-1); %avection
Jd(in) = -p.D*(P(in)-P(in-1))/p.dz; %diffusion

J = Ja + Jd;

dPdt = -(J(2:(p.n+1))-J(1:p.n))/p.dz;


%% Nutrients
%No sinking only mixing. 
%boundary condition:
%No flux at the top, diffusion on the bottom
Jd_N(1)=0; %surface
Jd_N(p.n+1)=-p.D*(p.N_b-N(p.n))/p.dz; %bott1om

in = 2:(p.n); %intorier
Jd_N(in) = -p.D*(N(in)-N(in-1))/p.dz; %diffusion ONLY DIFF


N_mixing = -(Jd_N(2:(p.n+1))-Jd_N(1:p.n))/p.dz;

%% Detrius
% boundary condition #NO FLUX
Ja_D(1)=0; %surface
Jd_D(1)=0; %surface
Ja_D(p.n+1)=0; %bottom
Jd_D(p.n+1)=0; %bottom

in = 2:(p.n); %intorier
Ja_D(in) = p.w*D(in-1); %avection
Jd_D(in) = -p.D*(D(in)-D(in-1))/p.dz; %diffusion

J_D = Ja_D + Jd_D;

dDdt = -(J_D(2:(p.n+1))-J_D(1:p.n))/p.dz;

%%
% light influence on growth 

I=func_light(p.z,P',p);

% Growth due to nutrients

mu=p.mu_max*min(N./(p.Hn+N) , I./(p.Hi+I));

%% Change of plankton
%dPdt=growth-loss-sinking+mixing
P_growth=mu.*P;
P_loss=p.m.*P;
P_sink_mix=dPdt';

dPdt=P_growth-P_loss+P_sink_mix;

%% Change of nutrients
%dNdt=-uptake+recycling+mixing
N_up=p.alpha.*P_growth;
%N_re=p.alpha*p.eps.*P_loss;
N_detri=p.tau*D;
N_mix=N_mixing';


dNdt=-N_up+N_mix+N_detri;

%% change of detritus
D_in=p.alpha.*P_loss;
D_out=N_detri;
D_mix=dDdt';

dDdt=D_in-D_out+D_mix;

%% result
dydt=[dPdt',dNdt',dDdt']';

end 