clear all;close all;clc;
MSE=[0.01 0.027 0.05  0.07 0.1 0.15 0.2];
MSE_dB=10*log10(MSE);
figure
ax = gca; % current axes
ax.FontSize = 14; grid on
Linewidth = 2;
marker_size = 5;

load('ADDop.mat')
plot(MSE_dB,ADD,'bo--','linewidth',Linewidth,'markersize',marker_size )
hold on 
load('ADDrd.mat')
plot(MSE_dB,ADD,'md-','linewidth',Linewidth,'markersize',marker_size )
xlabel(' Attack energy (dB)')
ylabel('Average detection delay')
%xlim([-46.4 -15.9])
legend('optimum attack','random attack'); 
grid on 
set(gca, 'fontsize', 14);