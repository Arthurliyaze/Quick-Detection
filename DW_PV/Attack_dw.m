%Test the dw algorithm
%   Detailed explanation goes here
clc; clear;close all;
load system_data.mat;
% For replay
%load data1min_dw.mat;
%Attack_signal = [output.t(1:14001) output.vfre([1:9000,begin:begin+5000])...
    %output.vmag([1:9000,begin:begin+5000])];
%replay_theoratical;

Attack_signal = [output.t(1:14001) output.vfre([1:9000,1:1+5000])...
    output.vmag([1:9000,1:1+5000])];

Horizon = 5;
numSteps = Horizon/sampleTime+1;
time = sampleTime*(0:numSteps-1);
time = time';
attack_begin = 4.5;
attack_end = Horizon;
attack_trigger = zeros(numSteps,1);
attack_trigger(attack_begin/sampleTime+1:attack_end/sampleTime) = 1;
Attack_time = timeseries(attack_trigger,time);
%output = sim('pv_converter_dw_replay',Horizon);
load data5s_destable.mat;


%%
x_all = [output.id output.iq output.vd output.vq]';

% input u [vid, viq] 2x1
u_all = [output.vid output.viq]';

% output y [vfre, vmag] 2x1
y_all = [output.vfre output.vmag]';

% noise v [v1 v2] 2x1
v_all = [output.v1 output.v2]';

% dw e [e1 e2] 2x1
e_all = [output.e1 output.e2]';

% measurement z [z1 z2] 2x1
z_all = output.z';

x = x_all(:,begin:end);
u = u_all(:,begin:end);
y = y_all(:,begin:end);
v = v_all(:,begin:end);
e = e_all(:,begin:end);
z = z_all(:,begin:end);

t = output.t(begin:end);
clear x_all u_all y_all v_all e_all z_all

close all
figure
subplot(2,1,1)
plot(t(1:5000),y(1,1:5000))
hold on
grid on
ylim([58,62])
xlim([4.46,4.56])
plot(t(5001:end),z(1,5001:end),'color',[0.9290, 0.6940, 0.1250])
%plot(t(5001:end),y(1,5001:end))
xlabel('time (s)')
ylabel('\omega (Hz)')
legend('Actual signal before attack','Actual signal after attack')

subplot(2,1,2)
plot(t(1:5000),y(2,1:5000))
grid on
ylim([326.35,326.85])
xlim([4.46,4.56])
xlabel('time (s)')
ylabel('|V_g| (v)')
hold on
plot(t(5001:end),z(2,5001:end),'color',[0.9290, 0.6940, 0.1250])
%plot(t(5001:end),y(2,5001:end))
%legend('Actual signal before attack','False data','Actual signal after attack')
legend('Actual signal before attack','Actual signal after attack')


%% small deviation over the equilibrium state
delx = x-repmat(x_e',1,length(t));
delu = u-repmat(u_e',1,length(t));
delz = z-repmat(y_e',1,length(t));

out = lsim(kalmf,[delu;delz],0:Ts:(length(t)-1)*Ts,delx(:,1));
delxe = out(:,p+1:p+n)';


%% Calculate statistics
g_all = delxe(:,2:end)-Ad*delxe(:,1:end-1)-B*delu(:,1:end-1);
[U,V] = eig(g_all*g_all');
Whiten = U';
g = Whiten([3 4],:)*g_all;
l = 100; %slide window length
S = zeros(4,4,length(t)-l+1); %Sample covariance
% M = zeros(1,length(t));
% T1 = zeros(1,length(t)-1);
% T1(1) = 0;
% N = zeros(1,length(t));
% T2 = zeros(1,length(t)-1);
% T2(1) = 0;
Y = zeros(1,length(t));
T3 = zeros(1,length(t)-1);
T3(1) = 0;

load data1min_dw.mat Sigma;
% for i = l:length(t)-1
%     S(:,:,i) = cov(100*[g(:,i-l+1:i); e(:,i-l+1:i)]');
% end
% Sigma = mean(S,3); %Asymptotic sample variance 

for i = 5:length(t)-1
%     M(i) = g(:,i)'/Sigma1*g(:,i); %Chi-2
%     T1(i) = max(0,T1(i-1)+(M(i)-2)/sqrt(2*2));
%     r = [g(:,i-l+1:i); 1000*e(:,i-l+1:i)];
%     S = cov(r');
%     N(i) = l*(-log(det(S/Sigma2))+trace(S/Sigma2)-n);%Chi-n(n+1)/2
%     d = n*(n+1)/2;
%     %d = 12;
%     T2(i) = max(0,T2(i-1)+(N(i)-d)/sqrt(2*d));

    s = 100*[g(:,i);e(:,i-1)];
    Y(i) = s'/Sigma*s;
    % Use sample mean as bias
    M = 4;
    %M = mean(Y(1:i));
    %M = 4*(1-2/(9*4))^3;
    T3(i) = max(0,T3(i-1)+(Y(i)-M)/sqrt(2*M));

end


%% Check M,N value
% figure
% plot(t,M)
% title('Distance measure M');
% xlabel('time(s)')
% 
% figure
% plot(t,N)
% title('Distance measure N');
% xlabel('time(s)')

figure
plot(t,Y)
title('Distance measure Y');
xlabel('time(s)')

%% Check M,N distribution
% figure
% histogram(M,100,'Normalization','pdf')
% hold on
% xl=0:0.1:16;
% pdf1 = chi2pdf(xl,2);
% plot(xl,pdf1,'LineWidth',1.5)
% legend('Histogram of M','Chi-square pdf (degree 2)')
% title('Distribution of M')
% 
% figure
% histogram(N,50,'Normalization','pdf')
% hold on
% xl=0:0.1:60;
% pdf2 = chi2pdf(xl,d);
% plot(xl,pdf2,'LineWidth',1.5)
% legend('Histogram of N','Chi-square pdf (degree 10)')
% title('Distribution of N')

figure
histogram(Y,100,'Normalization','pdf')
hold on
xl=0:0.1:16;
pdf3 = chi2pdf(xl,4);
plot(xl,pdf3,'LineWidth',1.5)
legend('Histogram of M','Chi-square pdf (degree 4)')
title('Distribution of Y')

%% Hotelling t test for statistic 1
% l = 8; %slide window length
% T_2 = zeros(1,length(t)-l);
% for i = 1:length(t)-l
%     g_bar = mean(g_w(:,i:i+l-1),2);
%     S = 1/(l-1)*(g_w(:,i:i+l-1)-g_bar)*(g_w(:,i:i+l-1)-g_bar)';
%     T_2(i) = l*g_bar'/S*g_bar;
% end
% 
% F = (l-m)/m*T_2/(l-1);
% %close all
% % figure
% % plot(F,'o')
% figure
% histogram(F,200,'Normalization','pdf')
% hold on
% xl=0:0.1:10;
% pdf2 = fpdf(xl,m,l-m);
% plot(xl,pdf2)
% xlim([0,22])
% legend('Sample','F(2,6) pdf')
% title('Histogram of Hotelling t^2 statistic')

%% Check statistic
% figure
% plot(t(1:end-1),T1,'b','LineWidth',1.5)
% xlim([4,5])
% grid on;
% title('Statistic T_1');
% xlabel('time')
% ylabel('Detector 1 statistics')
% 
% figure
% plot(t(1:end-1),T2,'b','LineWidth',1.5)
% xlim([4,5])
% grid on;
% title('Statistic T_2');
% xlabel('time')
% ylabel('Detector 2 statistics')
close all
T = T3./(1:length(T3));
figure
plot(t(1:end-1),T,'b','LineWidth',2)
xlim([3,5])
grid on;
xlabel('time (s)')
ylabel('T')
axes('Position',[.25 .5 .25 .25])
box on
plot(t(1:end-1),T,'b-','LineWidth',2);
grid on
ylim([0.0,0.5])
xlim([4.45,4.55])