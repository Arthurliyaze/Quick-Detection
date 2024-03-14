clear all;close all;clc;
%check the GLR code of Isreal
load('H_IEEE14_mat.mat');%load H matrix
%% generate A
[m,n]=size(H);
snr=10;%signal-to-noise ratio set as 10
var_x=1;%state variable variance set as 1
var_e=var_x/(10^(snr/10));%measurement noise variance
Sigma_x=var_x*eye(n);
Sigma_e=var_e*eye(m);
Sigma_z=var_x*H*H'+Sigma_e;
Sigma_zr=inv(Sigma_z);
upsilion=0.0217;%the enegry of the attack vector sigma_a^2
k=1;%the starting instant for CUSUM set as 1
l=500;
%s=5;
theta=10;
%% attack vector and z
num_GLR=zeros(1,4);
for s=1:4
%s=4;
for iter=1:100
    [a,index_s]= OAV(upsilion,var_x,10,H,s);%optimum attack vector
    as=a(index_s);
    z=zeros(m,500);
    for ite=1:l%define observed measurement vector
        x=randn(n,1)*sqrt(var_x);
        e=randn(m,1)*sqrt(var_e);
        if ite<theta
            z(:,ite)=H*x+e;
        else
            z(:,ite)=H*x+a+e;
        end
    end
    %% GLR
    [Max_GLR,optimal_GLR_support ] = GLR_sparsity(z',Sigma_zr,s);%isreal check

    index_q=optimal_GLR_support;
    %error(iter)=sum((a_hat-a(index_s)).^2);
    %num=num+isequal(It,index_s);
    num_GLR(s)=num_GLR(s)+nnz(ismember(index_q,index_s));
end
num_GLR(s)=num_GLR(s)/(s*100);
end