clear all;close all;clc;
%check the OMP code with y and A provided in main1.m
load('H_IEEE14_mat.mat');%load H matrix
load('z.mat');
%% generate A
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
k=1;%the starting instant for CUSUM set as 1
l=500;
s=5;
theta=10;
sigma=0.15;
%% generate attack vector and y
num_OMP=zeros(1,4);
%for s=1:4
%for iter=1:10000
[a,index_s]= OAV(upsilion,var_x,10,H,s);%optimum attack vector
as=a(index_s);
%{
z=zeros(m,500);
for ite=1:l%define observed measurement vector
    r1=randn(n,1);
    r2=randn(m,1);
    %save random.mat r1 r2
    x=r1*sqrt(var_x);
    e=r2*sqrt(var_e);
    if ite<theta
        z(:,ite)=H*x+e;
    else
        z(:,ite)=H*x+a+e;
    end
end
%}
y=A*sum(z(:,k:l),2)/(l-k+1);
%% OMP with unknown s
r=zeros(m,m);%%%%%%%%%%%%%% r(:,t)=r_t-1 %%%%%%%%%%%%%%%
L=l-k+1;
It=[];%I0 empty
r(:,1)=y;%r0=y
t=0;%iteration
j=1:m;%index set
i=zeros(1,m);%column index
while 1
    t=t+1;
    [~,index_t]=max(abs(r(:,t)'*A(:,j)));%find i_t
    it=j(index_t);
    j(index_t)=[];%delete i_t from j
    It=sort([It,it]);%It={i_t}U{I_t-1}
    AIt=A(:,It);
    Pt=AIt*inv(AIt'*AIt)*AIt';
    Ptv=eye(m)-Pt;
    Ct=Ptv(1:(m-t),:);%first m-t row
    Sigma_t=Ct*Ct'/L;
    Tt=y'*Ct'*inv(Sigma_t)*(Ct*y);%test value
    lambda_t=2*gammaincinv(sigma,(m-t)/2,'upper');%threshold
    r(:,t+1)=Ptv*y;
    if Tt<lambda_t%t>5%
        break;
    end
end
a_hat=(AIt'*AIt)\AIt'*y;
a_est=zeros(54,1);
a_est(It)=a_hat;
%{
if length(It)>s
    It=It(1:s);
end
num=num+isequal(It,index_s);
num_OMP(s)=num_OMP(s)+nnz(ismember(It,index_s));
error(iter)=sum((a_hat-a(index_s)).^2);
%}
%end
%num_OMP(s)=num_OMP(s)/(s*10000);
%hist(error,20);
%end