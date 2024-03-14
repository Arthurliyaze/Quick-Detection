% Calculate the KL diverge
clc; clear; close all;
load system_data.mat;
load data1min_dw.mat;
q = 2;
U = Whiten([3 4],:)';
mu0 = zeros(q+m,1);
Sigma0 = Sigma/100;
Sigma_e = Sigma0([3 4], [3 4]);
% Define the KL diverge function
kl = @(mu,sigma) 0.5*(mu'/Sigma0*mu + trace(sigma/Sigma0)+...
    log(det(Sigma0)/det(sigma)) -m-q);

%% FDI_D
a = [0.05;-0.05];
mu1 = [U'*K*a;0;0];
D1 = kl(mu1,Sigma0);
%% FDI_N
Sigma_a = 100*diag(3e-6*[10,1]);
Sigma_2 = Sigma0+blkdiag(U'*K*Sigma_a*K'*U,zeros(m,m));
%D2 = 0.5*(trace(Sigma_2/D)+log(det(D)/det(Sigma_2))-m-q);
D2 = kl(mu0,Sigma_2);

%% Replay
Ae = (Ad+B*L')*(eye(4)-K*C);
a = Ae';
b = zeros(4,2);
Q = B*E*B';
[X,~,~] = idare(a,b,Q,[],[],[]);
Sigma_3 = Sigma0+[2*U'*K*C*X*C'*K'*U, -U'*K*C*B*Sigma_e;...
    -Sigma_e*B'*C'*K'*U, zeros(m,m)];
D3 = kl(mu0,Sigma_3);

%% Destabilization
load destable_xhat.mat;
Ap = [diag([1.5 1.5]) zeros(2,2)];
Pa = B*Ap*P*Ap'*B' + Ad*P*Ap'*B' + B*Ap*P*Ad';
mu4 = [U'*K*C*B*Ap*delxe(:,end);0;0];
Sigma_4 = Sigma0 + blkdiag(U'*K*C*Pa*C'*K'*U,zeros(m,m));
D4 = kl(mu4,Sigma_4);
