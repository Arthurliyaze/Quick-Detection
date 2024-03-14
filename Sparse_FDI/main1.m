clear all;close all;clc;
%this code generates the matrix pos for fig1
%% generate A
load('H_IEEE14_mat.mat');%load H matrix
[m,n]=size(H);
snr=10;%signal-to-noise ratio set as 10
s=5;
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
k=1;%the starting instant for CUSUM set as 1
l=500;
%% compute the pos value
sigma=[.01 .04 .06 .1:.1:1];%the probability of false positive
pos=zeros(7,length(sigma));%pos at each iteration
for j=1:length(sigma)
    result=zeros(1,10000);
    for trial=1:10000
        %p=0.1;
        %theta=geornd(p);%attack time
        theta=10;
        [a,index_s]= OAV(upsilion,var_x,10,H,s);%optimum attack vector
        z=zeros(m,500);
        for i=1:l%define observed measurement vector
            x=randn(n,1)*sqrt(var_x);
            e=randn(m,1)*sqrt(var_e);
            if i<theta
                z(:,i)=H*x+e;
            else
                z(:,i)=H*x+a+e;
            end
        end
        y=A*sum(z(:,k:l),2)/(l-k+1);
        [~,t,~]=OMP(y,A,sigma(j),l,k);
        result(trial)=t;
    end
    [count,value]=hist(result,1:8);
    pos(:,j)=cumsum(count(1:7)')/10000;
end