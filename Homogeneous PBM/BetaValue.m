function Beta = BetaValue(X,NS)

global beta0

Beta=zeros(NS,NS);

for ii=1:NS
    
    XI=X(ii);
    DI=(6*XI/pi)^(1/3);
    
    for k=1:NS
        
        XK=X(k);
        DK=(6*XK/pi)^(1/3);
        
%         Beta1=2*kB*T/3/miu*(DI+DK)^2/DI/DK;
%         Beta2=(DI+DK<lamd)*4/3*sqrt(3*pi*epsilon/10/niu)*(DI+DK)^3;

%         Beta2=4/3*sqrt(3*pi*epsilon/10/niu)*(DI+DK)^3;

        Beta(ii,k)= beta0;
%         Beta(ii,k)= 1;
        
    end
          
end
end



