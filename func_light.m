function I=func_light(z,phi,p)

%I=p.I0*exp(-p.Kbg*z);

int=cumsum(phi.*p.dz,1)*p.k;

I=p.I0*exp(-p.Kbg*z-int);

% Make I a column vector:
I = I';
end 