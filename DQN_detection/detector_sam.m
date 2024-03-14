% this code calculates the ADD and PFA of Samrat's detector
% Yaze Li University of Arkansas
clear all; close all; clc;
%% Load system data
load('Measurements_data.mat'); % this matrix contains: z Chol_R Vm_true  del_true v n_meas n_samples nbus ;
load('Equations_inputs.mat');% this matrix contains: fbus_id tbus_id  nvi npi nqi npf nqf nbus G B qi bpq ppi;

%% initialization of the Dynamic model
N = 2*nbus-1;% # of state variables
Qk =1e-6*eye(N);% prediction noise covariance matrix
n_samples = 100;

%% tolerance for estimation
tol_max =  1e-4;
max_iter = 2000;

%%
ref_idx = 1;
ang_ref = del_true(ref_idx,:);
trials = 1000;
attack_time = zeros(trials,1);
detect_time = zeros(trials,1);
delay = NaN(trials,1);
fa = 0;

thres_my_test = 5e2; % Israel
thres_sam = 1e8;

Tk = zeros(3,trials);

for j1 = 1:trials
    j1 = j1
    attack_num = randi([1,n_meas],1,1);
    idx_inj = randperm(n_meas,attack_num);
    attack_mag = 0;
    attack_time(j1) = randi([1,50],1,1);
    
    %%
    Sigma_k = zeros(N,N);
    V_in = ones(nbus,1);
    ang_in = zeros(nbus-1,1);
    x_k_estimate = [ang_in;V_in]; %[zeros(nbus-1,1);zeros(nbus,1)];% [0.99*del_true(2:end,1); 0.99*Vm_true(:,1)] ; %
    Mk =  zeros(N,N);
    %% Parameter Fk and Gk model identification
    beta = 0.001 ;%0.001; %0.01 0.0009(Good predictions)
    alpha = 0.95;%0.95; %0.99;  alpha = 1./(1+beta) %
    
    x_k_predict = x_k_estimate ; %zeros(N,1); % CHANGE
    a_k_min_1 = zeros(N,1);
    b_k_min_1 = zeros(N,1);
    
    %%
    % Sam --
    test_stat_sam = 0;
    detect_indicator = 0;
    
    for k = 1:n_samples
        Sample = k;
        
        %% noisy measurements at time k;
        v_meas = v(k);
        ang_ref_bus = ang_ref(k);
        z_true  = measurement_fcn_h(Vm_true(:,k),fbus_id,v_meas,npi,nqi,npf,nqf,nbus,G,del_true(:,k),B,qi,tbus_id,bpq,ppi );
        %     z_true = z(:,k);
        RandData_n = randn(n_meas,1);
        noise =  Chol_R'*RandData_n;
        
        zk = z_true + noise; % Change this later
        
        if (k >= attack_time(j1) && k <= 100)
            zk( idx_inj)= zk( idx_inj)+ attack_mag.*sqrt(diag(Rk(idx_inj, idx_inj)));
            %zk( idx_inj)= 0; % DOS attack
        end
        %%  Detection
        
        if (k~= 1)
            
            ang_pred = [ang_ref_bus;x_k_predict(1:nbus-1)];
            V_pred = x_k_predict(nbus:end);
            
            z_pred  = measurement_fcn_h(V_pred,fbus_id,v_meas,npi,nqi,npf,nqf,nbus,G,ang_pred,B,qi,tbus_id,bpq,ppi );
            H_pred = Jacobian(V_pred,fbus_id,nvi,npi,nqi,npf,nqf,nbus,G,ang_pred,B,qi,tbus_id,bpq,ppi);
            
            Nk = H_pred*Mk*H_pred';
            Res_var = Nk + Rk;
            res =  zk - z_pred;
            %% BD and sudden change detection
            
            sigma_N = sqrt(diag(Nk));
            sigma_R = sqrt(diag(Rk));
            std_dev_res =sigma_N +  sigma_R; % sqrt(diag(Nk)+diag(Rk)); % sigma_N +  sigma_R
            
            thres1 = 3.5*(std_dev_res); %gamma = 3.5;
            test1 = (abs(res)> thres1);
            clear idx;
            idx = find(test1 == 1);
            %% whitening
            [U1,Eeigens1,U2]= svd( Res_var);
            Whiten_mat = inv(U1*sqrt(Eeigens1));
            
            y = Whiten_mat*res;
            T_whiten = y'*y;
            
            % Sam --
            test_stat_sam = max(0, (test_stat_sam + T_whiten - n_meas)/sqrt(2*n_meas));
            detect_indicator = test_stat_sam > thres_sam ;
            % ------
            
            if detect_indicator == 1
                detect_time(j1) = k;
                if detect_time(j1) >= attack_time(j1)
                    delay(j1) = detect_time(j1) - attack_time(j1);
                else
                    fa = fa + 1;
                end
                if(~isempty(idx))
                    zk(idx) = z_pred(idx);
                end
                break
            end
            
        if k == attack_time(j1)-1
            Tk(1,j1) = T_whiten;
        end
        if k == attack_time(j1)
            Tk(2,j1) = T_whiten;
        end
        if k == attack_time(j1)+1
            Tk(3,j1) = T_whiten;
        end
        end
        
        %% ------------- Dynamic State estimation ----------------------
        %------------------------------------------------------------------------
        
        %% update h and H
        
        [x_k_estimate,Sigma_k,iter,tol,flag]= Estimation(zk,Mk,Rk,x_k_predict,v_meas,ref_idx,ang_ref_bus,tol_max,max_iter,N,k);
        %         flag = flag
        
        %% measurements after estimation
        ang_est = [ang_ref_bus;x_k_estimate(1:nbus-1)];
        V_est = x_k_estimate(nbus:end);
        
        zk_est= measurement_fcn_h(V_est,fbus_id,v_meas,npi,nqi,npf,nqf,nbus,G,ang_est,B,qi,tbus_id,bpq,ppi );
        
        %% Parameter Fk and Gk model identification
        F = alpha*(1+beta);
        
        if k==1
            x_k_predict = x_k_estimate ;
            x_k  = x_k_estimate ;
        end
        
        Gk = (1+beta)*(1-alpha)*x_k_predict-beta*a_k_min_1 + (1-beta)*b_k_min_1 ;
        a_k = alpha*x_k_estimate  +(1-alpha)*x_k_predict;
        b_k = beta*(a_k -a_k_min_1)+(1-beta)*b_k_min_1;
        
        a_k_min_1 = a_k;
        b_k_min_1 = b_k;
        
        Fk = diag(F*ones(N,1));
        
        %% state prediction at k = k+1;
        x_k_predict = Fk*x_k_estimate+ Gk;
        Mk = Fk*Sigma_k*Fk.' + Qk; % update covariance matrix of xk_predict
        
    end
    
end

%%
add = nanmean(delay);
pfa = fa/trials;
%%
close
figure
histogram(Tk(1,:));
hold on;
histogram(Tk(2,:));
histogram(Tk(3,:));
legend('pre','attack','post')