
function dN=NumCon(t,N,Eta0,Eta1,Beta,Chi,Gama,x,NS)


% global A


dN=zeros(NS,1);


for ii=1:NS
    
    sum1=0;
    sum2=0;
    sum3=0;
    
    for k=1:ii
        
        for j=k:ii
            
            sum1=sum1+(Eta0(ii,k,j)+Eta1(ii,k,j))*(1-0.5*(j==k))*Beta(k,j)*N(j)*N(k);
            
        end
   
    end
    
    for k=1:NS
        
        sum2=sum2+Beta(ii,k)*N(ii)*N(k);
        
    end
    
    for k=ii:NS
       
       sum3=sum3+Chi(ii,k)*Gama(k)*N(k);
       
    end
    
    dN(ii,1)=sum1-sum2+sum3-Gama(ii)*N(ii)*(ii~=1);
%     dN(ii,1)=sum1-sum2+sum3-Gama(ii)*N(ii);
    
end


disp(t);

end







