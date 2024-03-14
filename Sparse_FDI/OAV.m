function [a,index_s] = OAV( upsilion, var_x, snr, H, s )
%Optimum Attack Vector from adversary's perspective
%   a is the optimum attack vector
%   index_s is the index of measurement being attacked
%   upsilion is the enegry of the attack
%   var_x is the state variable variance
%   snr is signal-to-noise ratio
%   H is the measurement Jacobian matrix
%   s is the spasity
%rng(1);
[m,n]=size(H);
a=zeros(m,1);
%var_x=1;%state variable variance set as 1
var_e=var_x/(10^(snr/10));%measurement noise variance
Sigma_x=var_x*eye(n);
Sigma_e=var_e*eye(m);
Sigma_z=var_x*H*H'+Sigma_e;
Sigma_zr=inv(Sigma_z);

%index_s=1:s;%for test
index_s=sort(randsample(m,s))';
lambda_s=Sigma_zr(:,index_s);%columns of invz
Ks=Sigma_x*H'*lambda_s;
phi_s=Sigma_zr(index_s,index_s);
[V,D]=eig(phi_s,Ks'*Ks);
[~,index]=min(diag(D));
umin=V(:,index);
a_s=sqrt(upsilion)/norm(Ks*umin)*umin;
a(index_s)=a_s;
end