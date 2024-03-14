clear all;close all;clc;
%% IEEE57
figure
ax = gca; % current axes
ax.FontSize = 14; grid on
Linewidth = 2;
marker_size = 5;
clear ADD 
clc
load('ADD57.mat')
plot(PFA57(1,:),ADD57(1,:),'rd--','linewidth',Linewidth,'markersize',marker_size )
hold on 
plot(PFA57(2,:),ADD57(2,:),'b.--','linewidth',Linewidth,'markersize',marker_size )
plot(PFA57(3,:),ADD57(3,:),'cs-','linewidth',Linewidth,'markersize',marker_size )
xlabel('PFA')
ylabel('ADD')
xlim([0 0.15])
legend('sparsity = 2','sparsity = 6',  ' sparsity = 15'); 
title('IEEE bus 57')
grid on 
set(gca, 'fontsize', 14);