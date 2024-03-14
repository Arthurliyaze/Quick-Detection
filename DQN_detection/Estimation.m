
function [x_k_estimate,Sigma_k,iter,tol,flag]= Estimation(zk,Mk,Rk,x_k_predict,v,ref_idx, ang_ref_bus,tol_max,max_iter,N, k)
% Equations_inputs.mat includes: fbus_id, tbus_id,...
% nvi, npi, nqi, npf, nqf, nbus, G, B ,qi bpq ppi;
load('Equations_inputs.mat');% this matrix contains: fbus_id tbus_id  nvi npi nqi npf nqf nbus G B qi bpq ppi;

del= zeros(nbus,1); % initial bus angles at time 0;
del(2:end) = x_k_predict(1:nbus-1);
del(ref_idx) = ang_ref_bus;
V = x_k_predict(nbus:end);

x_k_estimate = x_k_predict;

flag = 1;
iter = 0;
tol = 5;
% Kk = zeros(N,length(zk));

while(tol > tol_max)%5e-4
    
    hk = measurement_fcn_h(V,fbus_id,v,npi,nqi,npf,nqf,nbus,G,del,B,qi,tbus_id,bpq,ppi);
    %Measurement Jacobian, H..--------------------------------------------
    H = Jacobian(V,fbus_id,nvi,npi,nqi,npf,nqf,nbus,G,del,B,qi,tbus_id,bpq,ppi);
    
    if (k == 1)
        %% State filtering (estimation) at k = k+1;
        Sigma_k = inv(H'*inv(Rk)*H);
        %Kk = Sigma_k*H'*inv(Rk)  ; % Mk*H.'*inv(H*Mk*H.'+Rk);        % update covariance matrix of xk_predict
        dE = Sigma_k*(H'*inv(Rk)*(zk-hk));
        
    else
        %% State filtering (estimation) at k = k+1;
        Sigma_k = inv(H'*inv(Rk)*H+inv(Mk));
        %Kk = Sigma_k*H'*inv(Rk)  ; % Mk*H.'*inv(H*Mk*H.'+Rk);        % update covariance matrix of xk_predict
        dE = Sigma_k*(H'*inv(Rk)*(zk-hk)- inv(Mk)*(x_k_estimate-x_k_predict) );
    end
    x_k_estimate = x_k_estimate + dE;
    
    iter = iter + 1;
    tol = max(abs(dE));
    
    del(2:end) = x_k_estimate(1:nbus-1);
    V = x_k_estimate(nbus:end);
    
    if iter> max_iter
        break;
        flag = 0; % did not converge
    end
    
end
%%
% Sigma_k =(eye(N)-Kk*H)*Mk;

