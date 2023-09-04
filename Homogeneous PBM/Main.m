

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all
clc
tic
global beta0 N0 V s0 A1 A2 TAU

moment1=2;
moment2=3;

A1=moment1/3;
A2=moment2/3;

NS=50;
N0=1;
beta0=1e-3;
s0=1;
PHI_Inf=0.1*100;
TAU=10;
RT=TAU/N0/beta0;
tspan=[0,RT];

V=PHI_Inf^2*N0^2*beta0/2/s0;

v1=1e-4*V/N0;
v2=100*V/N0;
x1=log10(v1);
x2=log10(v2);

x=x1:(x2-x1)/NS:x2;
xm=(x(1:NS)+x(2:NS+1))/2;
x(1)=xm(1);
x(NS+1)=xm(NS);
v=10.^x;
vm=10.^xm;

dv=v(2:NS+1)-v(1:NS);

xx=vm;%%%%% unit  meter
eta=vm*N0/V;

for j=1:NS
    
    Q0(j)=integral(@InitialFun,v(j),v(j+1)); %%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

[Eta0,Eta1]=EtaValue(xx,NS);
B1=BValue(A1,xx,NS);
B2=BValue(A2,xx,NS);
Chi=ChiValue(B1,B2,xx,NS);

Beta=BetaValue(xx,NS);
Gama=GamaValue(xx,NS);

options=odeset('RelTol',1e-5,'MaxStep',10);

[t,Q]=ode45(@NumCon,tspan,Q0,options,Eta0,Eta1,Beta,Chi,Gama,xx,NS);

[m,n]=size(Q);
QQ=Q(m,:);
NT=sum(Q,2);
q=QQ./dv;

dpm=(vm*6/pi).^(1/3);
m4=sum(QQ.*dpm.^4);
m3=sum(QQ.*dpm.^3);
m2=sum(QQ.*dpm.^2);
m0=sum(QQ.*dpm.^0);

d43=m4./m3;
d32=m3./m2;
d30=m3./m0;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

NSS=1000;
X=x1:(x2-x1)/NSS:x2;
XM=(X(1:NSS)+X(2:NSS+1))/2;
VM=10.^XM;

Eta=N0*VM/V;

PHI=PHI_Inf*(1+PHI_Inf*tanh(PHI_Inf*TAU/2))/(PHI_Inf+tanh(PHI_Inf*TAU/2));

Nt=PHI*N0;

Phi0=exp(-Eta);

Phi=PHI^2*exp(-Eta.*PHI);

v1=v(1);
v2=v(NS+1);

M0=integral(@(v)MomentFun(v,0,PHI),v1,v2);
M2=integral(@(v)MomentFun(v,2,PHI),v1,v2);
M3=integral(@(v)MomentFun(v,3,PHI),v1,v2);
M4=integral(@(v)MomentFun(v,4,PHI),v1,v2);



% 
% for jj=1:NS
%     
% Nt(jj)=integral(@(v)MomentFun(v,0,PHI),v(jj),v(jj+1)); 
% 
% end
% 
% M0=sum(Nt.*dpm.^0)
% M2=sum(Nt.*dpm.^2)
% M3=sum(Nt.*dpm.^3)
% M4=sum(Nt.*dpm.^4)


D30=M3/M0;
D32=M3/M2

% Error1_d30=abs(d30-D30)/D30*100
Error2_d32=abs(d32-D32)/D32*100



semilogx(Eta,Phi0,'b--','LineWidth',1.5)
hold on
semilogx(Eta,Phi/PHI^2,'r-','LineWidth',1.5)
hold on
phi_t=sum(QQ)/sum(Q0);
VV=sum(QQ.*vm);
phi=q*VV/(sum(Q0))^2;
semilogx(eta(2:NS),phi(2:NS)/phi_t^2,'bs')

% axis([1e-4 1e2 0 1.05])


toc
