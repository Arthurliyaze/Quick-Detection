%Simulate the pv system and get simulation data
%% Set up
clc; clear;close all;
sampleTime = 1e-6;
pmuTime = 0.5*1e-3; % 2kHz
Horizon = 60; % 1min
numSteps = Horizon/sampleTime+1;
time = sampleTime*(0:numSteps-1);
time = time';
%output = sim('pv_converter',Horizon);
output = sim('pv_converter_dw',Horizon);

%% Get Simulink Output
% System size
m=2;
n=4;
p=2;

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

%% Check stable time
close all
figure
plot(y(1,:))
%xlim([1,200])

%%
begin = 2/pmuTime+1;  %after the system is stable
x = x_all(:,begin:end);
u = u_all(:,begin:end);
y = y_all(:,begin:end);
v = v_all(:,begin:end);
e = e_all(:,begin:end);

t = output.t(begin:end); % time index
clear x_all u_all y_all v_all e_all

%% equilibrium state
x_e = mean(x,2)';
x_e(1) = -150;
x_e(2) = 0;
u_e = mean(u,2)';
y_e = mean(y,2)';

% small deviation over the equilibrium state
delx = x-repmat(x_e',1,length(t));
delu = u-repmat(u_e',1,length(t));
dely = y-repmat(y_e',1,length(t));

delx(3:4,:) = delx(3:4,:)*1000;
delu = delu*100;
dely = dely*1000;

%%
%save data1min.mat