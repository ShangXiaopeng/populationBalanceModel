function B=BValue(A,x,NS)

B=zeros(NS,NS);


for ii=1:NS-1
    
    XI=x(ii);
    XI1=x(ii+1);
    
    for k=ii:NS
        
        XK=x(k);
        v=XK/2;
%         B(ii,k)=v^A*2*(v>=XI && v<XI1);  
%         f=inline('v^A*2/XK');
        B(ii,k)=integral(@(v)v.^A*2/XK,XI,XI1);
        
    end

end

end