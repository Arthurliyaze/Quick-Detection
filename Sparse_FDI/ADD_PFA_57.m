clear all;close all;clc;
%this is code for ADD bus57
%% generate A
load('H_IEEE57.mat');%load H matrix
[m,n]=size(H);
snr=10;%signal-to-noise ratio set as 10
var_x=1;%state variable variance set as 1
var_e=var_x/(10^(snr/10));%measurement noise variance
Sigma_x=var_x*eye(n);
Sigma_e=var_e*eye(m);
Sigma_z=var_x*H*H'+Sigma_e;
Sigma_zr=inv(Sigma_z);
[~,D,V]=svd(Sigma_z);
U=V';
A=sqrt(inv(D))*U;
upsilion=0.0217;%the enegry of the attack vector sigma_a^2
%% compute ADD
sigma=0.01;
sparsity=[2,6,15];
beta=0.01:0.01:0.15;
p0=0.1;
%threshold=-log(beta*p0);

log_alpha = [30:4:50 60:2:70];
alpha_vec = exp(log_alpha);
%PFA = 1./alpha_vec;
threshold =log(alpha_vec/p0);

theta=geornd(p0)+1;
theta_hat=0;
for i=1:3
    loop=i
    s=sparsity(i);
for j=1:length(threshold)
    B=threshold(j);
for iter=1:1000
[a,index_s]= OAV(upsilion,var_x,10,H,s);%optimum attack vector
z=zeros(m,theta+100);
J=zeros(1,theta+100);
for l=1:theta+100
    %get measurement z
    x=randn(n,1)*sqrt(var_x);
    e=randn(m,1)*sqrt(var_e);
    if l<theta
        z(:,l)=H*x+e;
    else
        z(:,l)=H*x+a+e;
    end
    J=zeros(1,l);
    for k=1:l
        L=l-k+1;
        y=A*sum(z(:,k:l),2)/L;
        w=sum(z(:,k:l),2)/L;
        [a_est,t,It]=OMP(y,A,sigma,l,k);
        J(k)=L*(w'*Sigma_zr*a_est-0.5*a_est'*Sigma_zr*a_est);
    end
    GLR=max(J);
    if GLR>=B
        theta_hat=l;
        break;
    end
end
error(iter)=theta_hat-theta;
end
ADD(i,j)=mean(error(error>=0));
PFA(i,j)=length(error(error<0))/1000;
end
end