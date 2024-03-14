clear all;close all;clc;
%this code generates the matrix for fig2
%% generate A
load('H_IEEE14_mat.mat');%load H matrix
[m,n]=size(H);
snr=10;%signal-to-noise ratio set as 10
s=1:5;%actual sparsity changes
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
%% compute the support recovery rate
sigma=0.05;%the probability of false positive
srr=zeros(3,length(s));%support recovery rate for 3 method
for j=2%1:length(s)
    number=zeros(3,1);%count number that recovered
    for trial=1:5000
        theta=10;
        [a,index_s]= OAV(upsilion,var_x,10,H,s(j));%optimum attack vector
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
        [~,index1]=OMPK(y,s(j),A);%support from OMPK
        [~,~,index2]=OMP(y,A,sigma,l,k);%from OMP
        if length(index2)>=s(j)
            index2=index2(1:s(j));
        end 
        %
        number(1)=number(1)+isequal(index1,index_s);%count OMPK
        number(2)=number(2)+isequal(index2,index_s);%count OMP
        %number(3)=number(3)+isequal(index3,index_s);%count GLR
    end
    srr(:,j)=number/5000;
end