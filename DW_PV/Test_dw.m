%Test the dw algorithm
%   Detailed explanation goes here
clc; clear;close all;
output = sim('pv_converter_dw',0.5);

%%
% state x [id, iq, vdg, vqg] 4x1
x_all = [output.id output.iq output.vd output.vq]';

% input u [vid, viq] 2x1
u_all = [output.vid output.viq]';

% output y [vfre, vmag] 2x1
y_all = [output.vfre output.vmag]';

% noise v [v1 v2] 2x1
v_all = [output.v1 output.v2]';

% dw e [e1 e2] 2x1
e_all = [output.e1 output.e2]';

begin = 2001;
x = x_all(:,begin:end);
u = u_all(:,begin:end);
y = y_all(:,begin:end);
v = v_all(:,begin:end);
e = e_all(:,begin:end);

t = output.t(begin:end);
clear x_all u_all y_all v_all e_all

%% equilibrium state
id_e = -150;
iq_e = 0;
vd_e = -321;
vq_e = 0;
vid_e = -0.865;
viq_e = -0.103;
fre_e = 60;
vmag_e = 321;
x_e = [id_e iq_e vd_e vq_e];
u_e = [vid_e viq_e];
y_e = [fre_e vmag_e];

% small deviation over the equilibrium state
delx = x-repmat(x_e',1,length(t));
delu = u-repmat(u_e',1,length(t));
dely = y-repmat(y_e',1,length(t));

% Calculate system matrix
% Calculate C
C = (dely-v)/delx;

% Calculate Ad B
delx_p = delx(:,2:end);
delxu = [delx(:,1:end-1);delu(:,1:end-1)];
AB = delx_p/delxu;
Ad = AB(:,1:4);
B = AB(:,5:6);

% Design Kalman Filter for MIMO Plant
BG = [B eye(4)];
DH = zeros(2,6);
Ts = -1;
sys = ss(Ad,BG,C,DH,Ts); % Discrete
sys.InputName = {'u1','u2','w1','w2','w3','w4'};
sys.OutputName = {'y1','y2'};

%known = [1 2];
%sensors = [1 2 3 4];
W = 1e-6*eye(4);
V = 5e-7*eye(2);
%N = 0;

% Discrete-Time Estimation
[kalmf,L,~,Mx,~] = kalman(sys,W,V);%,N,sensors,known);
out = lsim(kalmf,[delu;dely-v],0:length(t)-1,delx(:,1));
%delye = out(:,1:2)';
delxe = out(:,3:6)';

clear AB BG DH

%% Calculate two statistics
% Calculate residual
r = delx(:,2:end)-Ad*delx(:,1:end-1)-B*delu(:,1:end-1);

G = Mx;
[P,~,~,~] = idare(Ad',C',W,V,[],[]);
M = zeros(4,4,length(t)-1);
N = zeros(2,4,length(t)-1);
m = zeros(1,length(t)-1);
n = zeros(1,length(t)-1);
for i = 1:length(t)-1
    M(:,:,i) = r(:,i)*r(:,i)'-G*(C*P*C'+V)*G';
    N(:,:,i) = e(:,i)*r(:,i)';
    m(i) = trace(sum(abs(M),3))/i;
    n(i) = sum(sum(abs(N),3),'all')/i;
end

close all
figure
plot(t(2:end),m);
grid on;
title('The absolute trace of statistic 1 when no attack');
xlabel('time(s)')

figure
plot(t(2:end),n);
grid on;
title('The absolute sum of statistic 2 when no attack');
xlabel('time(s)')
