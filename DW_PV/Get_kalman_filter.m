% Get the Kalman filter
clc; clear; close all;
load data1min.mat
load est_matrix.mat

%% Design Kalman Filter for MIMO Plant
Ad = A_e;
B = B_e;
C = C_e;
BG = [B eye(n)];
DH = zeros(p,p+n);
Ts = pmuTime;
sys = ss(Ad,BG,C,DH,Ts); % Discrete
sys.InputName = {'u1','u2','w1','w2','w3','w4'};
sys.OutputName = {'y1','y2'};
known = [1 2];
sensors = [1 2];

W = 1e-6*eye(n); % Process noise covariance
V = 5e-7*eye(p); % Measurement noise covariance
E = 1e-6*eye(m); % Water mark covariance

% Discrete-Time Estimation
[kalmf,L,~,Mx,~] = kalman(sys,W,V,[],sensors,known);
out = lsim(kalmf,[delu;dely],0:Ts:(length(t)-1)*Ts,delx(:,1));
delxe = out(:,p+1:p+n)';

K = Mx;
[P,~,~,~] = idare(Ad',C',W,V,[],[]);
R = C*P*C'+V;

clear BG DH A_e B_e C_e

%% Check state estimation results
close all
i=3;
figure
plot(t,delx(i,:));
hold on;
grid on;
plot(t,delxe(i,:));
legend('True','Estimated')
xlim([2,5])
title('Small deviation of I_d over the equilibrium state');
xlabel('time')
ylabel('A')

%%
%save system_data.mat