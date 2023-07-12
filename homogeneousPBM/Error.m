% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear all
clc

set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'DefaultAxesFontName','Times New Roman') 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

Sections=[
10
20
30
40
50
];

MC30PHI01=[
8.315
2.215
1.028
0.580
0.367    
];

MC32PHI01=[
3.494
0.913
0.416
0.235
0.150    
];

MC30PHI1=[
5.889
1.485
0.675
0.381
0.243    
];

MC32PHI1=[
0.116
0.012
0.012
0.007
0.004    
];

MC30PHI10=[
5.866
1.483
0.674
0.380
0.241
];

MC32PHI10=[
0.105
0.014
0.013
0.008
0.005    
];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

subplot(3,1,1),semilogy(Sections,MC30PHI01,'b-+','LineWidth',1.5,'MarkerSize',8)
hold on
semilogy(Sections,MC32PHI01,'k-s','LineWidth',1.5,'MarkerSize',8)

ytick=10.^(-1:1);
set(gca,'YTick',ytick)
set(gca,'FontSize',10.5,'LineWidth',1.5,'TickLength',[0.02 0.02],'YMinorTick','off')
axis([0 60 0.1 10])
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

subplot(3,1,2),semilogy(Sections,MC30PHI1,'b-+','LineWidth',1.5,'MarkerSize',8)
hold on
semilogy(Sections,MC32PHI1,'k-s','LineWidth',1.5,'MarkerSize',8)

ytick=10.^(-3:1);
set(gca,'YTick',ytick)
set(gca,'FontSize',10.5,'LineWidth',1.5,'TickLength',[0.02 0.02],'YMinorTick','off')
axis([0 60 0.001 10])

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

subplot(3,1,3),semilogy(Sections,MC30PHI10,'b-+','LineWidth',1.5,'MarkerSize',8)
hold on
semilogy(Sections,MC32PHI10,'k-s','LineWidth',1.5,'MarkerSize',8)

ytick=10.^(-3:1);
set(gca,'YTick',ytick)
set(gca,'FontSize',10.5,'LineWidth',1.5,'TickLength',[0.02 0.02],'YMinorTick','off')
axis([0 60 0.001 10])



set(gcf,'Units','centimeters','Position',[10 0.5 8.5 21])
% set(gca,'Position',[0.14 0.13])







