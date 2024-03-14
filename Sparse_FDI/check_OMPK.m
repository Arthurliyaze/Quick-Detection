clear all;close all;clc;%to be done
%check the OMPK code with A and y provided
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
[~,D,V]=svd(Sigma_z);
U=V';
A=sqrt(inv(D))*U;
upsilion=0.0217;%the enegry of the attack vector sigma_a^2
k=1;%the starting instant for CUSUM set as 1
l=500;
theta=10;
%% attack vector and y
num_OMPK=zeros(1,5);
for s=1:5
for iter=1:1
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
    y=A*sum(z(:,k:l),2)/(l-k+1);
    %% OMPK s=5
    r=zeros(m,s+1);%%%%%%%%%%%%%% r(:,t)=r_t-1 %%%%%%%%%%%%%%%
    It=[];%I0 emptyC
    r(:,1)=y;%r0=y
    t=0;%iteration
    j=1:m;%index set
    %Pt=zeros(m,m);
    while t<s
        t=t+1;
        [~,index_t]=max(abs(r(:,t)'*A(:,j)));%find i_t
        it=j(index_t);
        j(index_t)=[];%delete i_t from j
        It=sort([It,it]);%It={i_t}U{I_t-1}
        AIt=A(:,It);
        Pt=AIt*inv(AIt'*AIt)*AIt';
        r(:,t+1)=(eye(m)-Pt)*y;
    end
    a_hat=inv(AIt'*AIt)*AIt'*y;
    %error(iter)=sum((a_hat-a(index_s)).^2);
    %num=num+isequal(It,index_s);
    num_OMPK(s)=num_OMPK(s)+nnz(ismember(It,index_s));
end
num_OMPK(s)=num_OMPK(s)/(s*1);
end
%% show error
%hist(error,20);