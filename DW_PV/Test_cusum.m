% Test the normalized rao-cusum detector
clear; close all; clc;
load system_data.mat;
%rng(1,'philox')
replay_theoratical;
% seed 3 negative bias
% seed 2 small bias
% seed 1 positive bias
m = 4;
n = 5*2000;
%%
trials = 1000;
delay = zeros(trials,1);
attack_begin = 4.5;
beta = 0.067;
T3 = zeros(trials,n);
tic
for i = 1:trials
    i
Y1 = chi2rnd(m,1,n/2);
Y2 = chi2rnd(m+delSigma,1,n/2);
Y = [Y1,Y2];

% figure
% histogram(Y,100,'Normalization','pdf')
% hold on
% xl=0:0.1:20;
% pdf1 = chi2pdf(xl,m);
% plot(xl,pdf1,'LineWidth',1.5)
% legend('Sample','Chi-square pdf (degree m)')
% title('Histogram of M')

T1 = zeros(1,n);
T2 = zeros(1,n);
%T3 = zeros(trials,n);
for k = 1:n-1
    %M = mean(Y(1:k));
    T1(k+1) = max(0,T1(k)+(Y(k+1)-m)/sqrt(2*m));
    %T2(k+1) = T2(k)+(Y(k+1)-m)/sqrt(2*m);
    %T3(k+1) = max(0,T2(k)+(Y(k+1)-1.2*m)/sqrt(2*1.2*m));
end
    T2 = T1./(1:length(T1))+0.04;
    T3(i,:) = T2;
    detect = find(T2(3501:end)>=beta, 1 ) + 3500;
    attack  = (attack_begin-2)/pmuTime;
    delay(i) = detect-attack;
end
toc
%%
indices = find(delay>0);
ADD = mean(delay(indices))
PFA = 1-length(indices)/trials
%%
% figure
% plot(1:n,T1,'b','LineWidth',2)
% grid on
% hold on
% plot(1:n,T2,'r','LineWidth',2)
% %plot(1:n,T3,'k','LineWidth',2)
% %legend('Theoratical mean as bias','Sample mean as bias','1.15 Theoratical mean')
% legend('Theoratical mean as bias','No max')
%%
% close all
% figure
% plot(2+pmuTime:pmuTime:7,0.04+T1./(1:n),'b','LineWidth',2)
% xlim([3,7])
% grid on;
% xlabel('time (s)')
% ylabel('T')
% axes('Position',[.25 .5 .25 .25])
% box on
% plot(2+pmuTime:pmuTime:7,T1./(1:n),'b','LineWidth',2)
% grid on
% ylim([0.0,0.5])
% xlim([4.49,4.52])