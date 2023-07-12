function moment= MomentFun(v,m,PHI)

global N0 V


Eta=v*N0/V;

Phi=PHI^2.*exp(-Eta*PHI);


nt=N0^2/V*Phi;
d=(6*v/pi).^(1/3);
moment=d.^m.*nt;

end





