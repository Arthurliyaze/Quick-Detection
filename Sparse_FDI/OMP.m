function [a_est,t,It] = OMP( y,A,sigma,l,k )%unfinished
%OMP with unknown s
%   a_est is the estimator of a
%   t is the iteration that stops
%   It is the estimate support of attack vector 
%   y is m*1 vector used in CUSUM
%   A'A is (sigma_z)^(-1) m*m
%   sigma is the probability of false positive
%   l is the sampling instant
%   k is the starting instant for CUSUM

m=length(y);%get m
a_est=zeros(m,1);
r=zeros(m,m);%%%%%%%%%%%%%% r(:,t)=r_t-1 %%%%%%%%%%%%%%%
L=l-k+1;
It=[];%I0 empty
r(:,1)=y;%r0=y
t=0;%iteration
j=1:m;%index set
while t<8
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
a_est(It)=a_hat;