% this code plots learning curves for phi = 0.1 and different m
% Yaze Li, University of Arkansas
clear all; %close all; clc;
%% Load data
Data1 = importdata('run-DQN_9-tag-episode_reward.csv');     %m=4
Data2 = importdata('run-DQN_10-tag-episode_reward.csv');    %m=2
Data3 = importdata('run-DQN_11-tag-episode_reward.csv');    %m=1
%r1 = max(0,Data1.data);
%r2 = max(0,Data2.data);
%r3 = max(0,Data3.data);
r1 = Data1.data;
r2 = Data2.data;
r3 = Data3.data;
%% Moving aaverage as tensorboard
y1 = smoothdata(r1(:,3),'gaussian',20);
y2 = smoothdata(r2(:,3),'gaussian',20);
y3 = smoothdata(r3(:,3),'gaussian',20);
%% Plot
fig = figure;
subplot(3,1,1)
plot(r3(:,2),r3(:,3),'-','color',[0.9020    0.9020    0.9020],'LineWidth',1);
grid on;
hold on;
plot(r3(:,2),y3,'k-','LineWidth',1);
grid on;
ylim([0,1])
title('w = 1')

subplot(3,1,2)
plot(r2(:,2),r2(:,3),'-','color',[0.9804    0.7647    0.7647],'LineWidth',1);
grid on;
hold on;
plot(r2(:,2),y2,'r-','LineWidth',1);
grid on;
ylim([0,1])
title('w = 2')

subplot(3,1,3)
plot(r1(:,2),r1(:,3),'-','color',[0.8431    0.8706    0.9804],'LineWidth',1);
grid on;
hold on;
plot(r1(:,2),y1,'b-','LineWidth',1);
grid on;
ylim([0,1])
title('w = 4')

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'step');
ylabel(han,'episode reward');