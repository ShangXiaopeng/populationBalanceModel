function n= InitialFun(v)

global N0 V

% n0=N0^2/V*exp(-N0/V*v).*v*log(10);
n0=N0^2/V*exp(-N0/V*v);

n=n0;

end





