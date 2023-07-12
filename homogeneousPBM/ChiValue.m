function Chi = ChiValue(B1,B2,x,NS)

global A1 A2
Chi=zeros(NS,NS);

for ii=1:NS
            
    if ii==1
        
            XI=x(ii);
            XI1=x(ii+1);
            
            for k=ii:NS
                
                XK=x(k);
%                 Chi(ii,k)=integral(@ (v) ChiFun1(v,XI,XI1,XK),XI,XI1);%% 

%                 v=XK/2;
%                 Chi(ii,k)=(XI1-v)/(XI1-XI)*2*(v>=XI && v<XI1);

                Chi(ii,k)=(B1(ii,k)*XI1^A2-B2(ii,k)*XI1^A1)/(XI^A1*XI1^A2-XI^A2*XI1^A1);
          
            end
            
    elseif ii==NS
            
            XI=x(ii);
            XI0=x(ii-1);
            for k=ii:NS
                XK=x(k);
                
%                 Chi(ii,k)=integral(@ (v) ChiFun0(v,XI0,XI,XK),XI0,XI);%% xNS break up
%                 v=XK/2;
%                 Chi(ii,k)=(v-XI0)/(XI-XI0)*2*(v>=XI0 && v<XI);

                Chi(ii,k)=(B1(ii-1,k)*XI0^A2-B2(ii-1,k)*XI0^A1)/(XI^A1*XI0^A2-XI^A2*XI0^A1);
                
            end
                    
    else
            
            XI=x(ii);
            XI1=x(ii+1);
            XI0=x(ii-1);
            
            for k=ii:NS
               
               XK=x(k);
               
               v=XK/2;
               
               if k==ii
                   
%                    Chi(ii,k)=integral(@ (v) ChiFun0(v,XI0,XI,XK),XI0,XI);%% i=2:NS-1, k=i
%                     Chi(ii,k)=(v-XI0)/(XI-XI0)*2*(v>=XI0 && v<XI);

                    Chi(ii,k)=(B1(ii-1,k)*XI0^A2-B2(ii-1,k)*XI0^A1)/(XI^A1*XI0^A2-XI^A2*XI0^A1);

                   
               else
                   
%                    Chi(ii,k)=integral(@ (v) ChiFun1(v,XI,XI1,XK),XI,XI1)+integral(@ (v) ChiFun0(v,XI0,XI,XK),XI0,XI);%% i=2:NS-1, k>i

%                     Chi(ii,k)=(XI1-v)/(XI1-XI)*2*(v>=XI && v<XI1)+(v-XI0)/(XI-XI0)*2*(v>=XI0 && v<XI);

                    Chi(ii,k)=(B1(ii,k)*XI1^A2-B2(ii,k)*XI1^A1)/(XI^A1*XI1^A2-XI^A2*XI1^A1)+(B1(ii-1,k)*XI0^A2-B2(ii-1,k)*XI0^A1)/(XI^A1*XI0^A2-XI^A2*XI0^A1);
                  
               end
                
            end
    end            
                
        

end
end