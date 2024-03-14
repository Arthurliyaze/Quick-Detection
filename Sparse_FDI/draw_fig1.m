% this code draw fig1 in the paper
clear all;close all;clc;
load('pos.mat');
sigma=[.01 .04 .06 .1:.1:1];
figure
Linewidth=1.5;
semilogx(sigma,pos(7,:),'m--','linewidth',Linewidth)
hold on
semilogx(sigma,pos(6,:),'b<-','linewidth',Linewidth)
semilogx(sigma,pos(5,:),'g-','linewidth',Linewidth)
semilogx(sigma,pos(4,:),'bd--','linewidth',Linewidth)
semilogx(sigma,pos(3,:),'k+-','linewidth',Linewidth)
semilogx(sigma,pos(2,:),'bo','linewidth',Linewidth)
semilogx(sigma,pos(1,:),'r','linewidth',Linewidth)
grid on
xlabel('Probability of false positive')
ylabel('Probability of stopping')
ylim([-0.1,1.1])
legend('7th iter','6th iter','5th iter','4th iter','3rd iter','2nd iter','1st iter')
grid on
set(gca, 'fontsize', 10);