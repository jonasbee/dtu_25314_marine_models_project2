function I=func_light_s(z,t,phi,p)

int=cumsum(phi.*p.dz,1)*p.k;

I=(p.I0*(1+cos(2*pi*t/365)))*exp(-p.Kbg*z-int);

% Make I a column vector:
I = I';
end 