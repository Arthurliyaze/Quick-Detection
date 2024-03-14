% QCD preformance analysis
clc; clear;close all;
%rng(1,'philox')
load system_data.mat;
load data1min_dw.mat Sigma;
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
t = output.t(begin:end);

clear x_all u_all y_all v_all e_all

%%
trials = 500;
delay = zeros(trials,1);
attack_begin = 4.5;
beta = 0.0595;
T3 = zeros(trials,length(t)-1);
tic
for k = 1:trials
    k
    z = z_all(:,begin:end);
    RandData_n = randn(2,length(z));
    noise =  1e3*V'*RandData_n;
    z = z + noise;

    delx = x-repmat(x_e',1,length(t));
    delu = u-repmat(u_e',1,length(t));
    delz = z-repmat(y_e',1,length(t));

    out = lsim(kalmf,[delu;delz],0:Ts:(length(t)-1)*Ts,delx(:,1));
    delxe = out(:,p+1:p+n)';

    g_all = delxe(:,2:end)-Ad*delxe(:,1:end-1)-B*delu(:,1:end-1);
    g = g_all([3 4],:);

    Y = zeros(1,length(t));
    T1 = zeros(1,length(t)-1);
    T1(1) = 0;

    for i = 5:length(t)-1
        s = 100*[g(:,i);e(:,i-1)];
        Y(i) = s'/Sigma*s;
        % Use sample mean as bias
        M = 4;
        T1(i) = max(0,T1(i-1)+(Y(i)-M)/sqrt(2*M));
    end

    T2 = T1./(1:length(T1));
    T3(k,:) = T2;
    detect = find(T2(3501:end)>=beta, 1 ) + 3500;
    attack  = (attack_begin-2)/pmuTime;
    delay(k) = detect-attack;
end
toc
%%
indices = find(delay>0);
ADD = mean(delay(indices))
PFA = 1-length(indices)/trials