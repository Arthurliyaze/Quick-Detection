function [T_whiten, x_k1_predict, Mk1, ak1, bk1] = step(current_step, attack_time, idx_inj, attack_mag, x_k_predict, Mk, ak, bk)
%STEP This function uses Israel and Samrat's code to calculate the statistic of Rao test
%and update the predicted state vector
% Yaze Li University of Arkansas

%   current_step:   The current time step during the training
%   attack_time:    The attack time [0,500]
%   idx_inj(1x2):        The attacked measurement indexes
%   attack_mag(1x2):     The attack magnitude
%   x_k_predict(1xN):    The predicted state vector from previous step
%   Mk:             The error covariance matrix from previous step
%   ak, bk(1xN):         The prediction parameter from previous step
%   T_whiten:            The statistic of Rao test
%   x_k1_predict(1xN):   The predicted state vector of current step
%   Mk1(1xN):            The error covariance matrix of current step
%   ak1, bk1(1xN):       The prediction parameter of current step

%% Load system data
load('Measurements_data.mat'); %#ok<*LOAD> % this matrix contains: z Chol_R Vm_true  del_true v n_meas n_samples nbus
load('Equations_inputs.mat');% this matrix contains: fbus_id tbus_id  nvi npi nqi npf nqf nbus G B qi bpq ppi

ref_idx = 1;
ang_ref = del_true(ref_idx,:); %#ok<*NODEF>
%% Change vector shape
attack_mag = attack_mag';
x_k_predict = x_k_predict';
ak = ak';
bk = bk';
%% initialization of the Dynamic model
N = 2*nbus-1;% # of state variables
Qk =1e-6*eye(N);% prediction noise covariance matrix
%Chol_Q = chol(Qk);

%% tolerance for Dynamic state estimation
tol_max =  1e-4;
max_iter = 2000;

%% Parameter Fk and Gk model identification
beta = 0.001 ;%0.001; %0.01 0.0009(Good predictions)
alpha = 0.95;%0.95; %0.99;  alpha = 1./(1+beta) %

%% noisy measurements at time k
k = current_step;
% Choose true state
v_meas = v(k);
ang_ref_bus = ang_ref(k);
% Calculate true measurement
z_true  = measurement_fcn_h(Vm_true(:,k),fbus_id,v_meas,npi,nqi,npf,nqf,nbus,G,del_true(:,k),B,qi,tbus_id,bpq,ppi );

RandData_n = randn(n_meas,1);
noise =  Chol_R'*RandData_n;
%noise = 0;
zk = z_true + noise; % noisy measurement

% Add attack vector
if (k >= attack_time)
    zk(idx_inj)= zk(idx_inj)+ attack_mag.*sqrt(diag(Rk(idx_inj, idx_inj)));
end

% Choose predicted state
ang_pred = [ang_ref_bus;x_k_predict(1:nbus-1)];
V_pred = x_k_predict(nbus:end);
% Calculate predicted measurement
z_pred  = measurement_fcn_h(V_pred,fbus_id,v_meas,npi,nqi,npf,nqf,nbus,G,ang_pred,B,qi,tbus_id,bpq,ppi );
% Calculate the Jacobian matrix
H_pred = Jacobian(V_pred,fbus_id,nvi,npi,nqi,npf,nqf,nbus,G,ang_pred,B,qi,tbus_id,bpq,ppi);
% Calculate the residual vector
v =  zk - z_pred;
% Calculate the covariance matrix of the residual vector
S = H_pred*Mk*H_pred' + Rk;
%% whitening
[U,D,~]= svd(S);
% W = inv(U*sqrt(D)); % Whitening matrix
% y = W*v;
y = (U*sqrt(D))\v; % whitening residual
T_whiten = y'*y;
%% Dynamic State estimation
[x_k_estimate,Sigma_k,~,~,~]= Estimation(zk,Mk,Rk,x_k_predict,v_meas,ref_idx,ang_ref_bus,tol_max,max_iter,N,k);

%% Parameter Fk and Gk model identification
F = alpha*(1+beta);
Fk = diag(F*ones(N,1));

if k==1
    x_k_predict = x_k_estimate ;
end

Gk = (1+beta)*(1-alpha)*x_k_predict-beta*ak + (1-beta)*bk ;
ak1 = alpha*x_k_estimate  +(1-alpha)*x_k_predict;
bk1 = beta*(ak1 -ak)+(1-beta)*bk;

%% Updata state prediction at k = k+1
x_k1_predict = Fk*x_k_estimate+ Gk;
Mk1 = Fk*Sigma_k*Fk.' + Qk; % update covariance matrix of xk_predict
%%  Change output to row vector
v = v';
x_k1_predict = x_k1_predict';
ak1 = ak1';
bk1 = bk1';
end