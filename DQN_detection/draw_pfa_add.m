% this code plots the ADD and PFA figure of Samrat's and Yaze's detector
% Yaze Li University of Arkansas
clear all; close all; clc;
%% Load data
load('sam_dos.mat'); 
add1 = ADD;
pfa1 = FAR;
load('yaze_dos.mat');
add2 = ADD;
pfa2 = PFA;
clear ADD FAR PFA
load('dis_dos.mat');
add3 = ADD;
pfa3 = PFA;
clear ADD PFA
%% Draw
figure;
plot(pfa1,add1,'bo-','LineWidth',1);
grid on;
hold on;
plot(pfa2,add2,'r*-','LineWidth',1);
plot(pfa3,add3,'k^-','LineWidth',1);
xlim([0,0.2])
ylim([0,2.1])
xlabel('PFA')
ylabel('ADD(seconds)')
legend('Normalized Rao-CUSUM','DQN (continuous state)','DQN (discrete state)')
axes('Position',[.25 .4 .25 .25])
box on
plot(pfa2,add2,'r*-','LineWidth',1);
grid on
ylim([0.019,0.024])
xlim([0.048,0.067])
axes('Position',[.6 .4 .25 .25])
box on
plot(pfa3,add3,'k^-','LineWidth',1);
grid on
xlim([0.071,0.10])
ylim([0.029,0.041])