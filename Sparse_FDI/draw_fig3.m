clear all;close all;clc;
%% IEEE14
figure
ax = gca; % current axes
ax.FontSize = 14; grid on
Linewidth = 2;
marker_size = 5;
load('ADD_14_6.mat')
plot(PFA,ADD,'rd--','linewidth',Linewidth,'markersize',marker_size )
hold on
clear ADD;clear PFA;
load('ADD_14_5.mat')
plot(FA,0.5*ADD,'b.--','linewidth',Linewidth,'markersize',marker_size )
clear ADD;clear FA;
load('ADD_14_2.mat')
plot(PFA,ADD,'cs--','linewidth',Linewidth,'markersize',marker_size )
xlim([0,0.15])
xlabel('Probability of false alarm')
ylabel('Average detection delay')
%xlim([0 0.15])
legend('sparsity = 6','sparsity = 5',' sparsity = 2'); 
%title('IEEE bus 14')
grid on 
set(gca, 'fontsize', 10);