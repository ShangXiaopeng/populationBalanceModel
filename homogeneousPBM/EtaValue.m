function [Eta0,Eta1] =EtaValue(x,NS)


global A1 A2

Eta0=zeros(NS,NS,NS);
Eta1=zeros(NS,NS,NS);

for ii=1:NS
    
    XI=x(ii);
    
    if ii==1
        XI1=x(ii+1);
        
        v=x(1)+x(1);
%         Eta1(ii,1,1)= (x(ii+1)-v)/(x(ii+1)-x(ii));
        Eta1(ii,1,1)=(v>=x(ii) && v<x(ii+1))*(v^A1*XI1^A2-v^A2*XI1^A1)/(XI^A1*XI1^A2-XI^A2*XI1^A1);
       
    elseif ii==NS
        
        XI0=x(ii-1);
       
        for k=1:ii
            
            for j=k:ii
                
                v=x(j)+x(k);
                
                if v>=x(ii)
                    Eta1(ii,k,j)= v/x(ii);
                end
                
                if v>=x(ii-1) && v<=x(ii)
%                     Eta0(ii,k,j)= (v-x(ii-1))/(x(ii)-x(ii-1));
                    Eta0(ii,k,j)=(v^A1*XI0^A2-v^A2*XI0^A1)/(XI^A1*XI0^A2-XI^A2*XI0^A1);

                end
            end
        
        end
        
    else
        
        XI0=x(ii-1);
        XI1=x(ii+1);
        
        for k=1:ii
            
            for j=k:ii
                
                v=x(j)+x(k);
                
                if v>=x(ii) && v<=x(ii+1)
%                     Eta1(ii,k,j)= (x(ii+1)-v)/(x(ii+1)-x(ii));
                    Eta1(ii,k,j)=(v^A1*XI1^A2-v^A2*XI1^A1)/(XI^A1*XI1^A2-XI^A2*XI1^A1);
                end
          
                if v>=x(ii-1) && v<=x(ii)
%                     Eta0(ii,k,j)= (v-x(ii-1))/(x(ii)-x(ii-1));
                    Eta0(ii,k,j)=(v^A1*XI0^A2-v^A2*XI0^A1)/(XI^A1*XI0^A2-XI^A2*XI0^A1);
                end
            end
      
        end
       
       
       
   end
    
    
    
    
    
end



end