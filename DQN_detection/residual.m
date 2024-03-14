% this code tests the statistic of Rao test calculation function
% Yaze Li University of Arkansas
clear all; close all; clc
%% Load system data
load('Measurements_data.mat'); % this matrix contains: z Chol_R Vm_true  del_true v n_meas n_samples nbus ;
load('Equations_inputs.mat');% this matrix contains: fbus_id tbus_id  nvi npi nqi npf nqf nbus G B qi bpq ppi;

%%
N = 2*nbus-1;
n_samples = 100;
idx_inj = [4 5];% attacked measurement indexes
attack_time = 69; 
attack_mag =  [-1500 1000];
bus_idx = 2; %
%% Initial values at time = 1

V_in = ones(nbus,1);
ang_in = zeros(nbus-1,1);
x_k_estimate = [ang_in;V_in]';
Mk =  zeros(N,N);
x_k_predict = x_k_estimate;
ak = zeros(1,N);
bk = zeros(1,N);

%% noisy measurements at time k;
T = zeros(1,n_samples);
for k = 1:n_samples
    current_step = k;
    [T(k),x_k1_predict,Mk1,ak1,bk1] = step(current_step, attack_time, idx_inj, attack_mag, x_k_predict, Mk, ak, bk);
    x_k_predict = x_k1_predict;
    Mk = Mk1;
    ak = ak1;
    bk = bk1;
end