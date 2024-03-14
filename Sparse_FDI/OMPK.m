function [a_hat,It] = OMPK( y,s,A )
%OMP with Known s
%   a_hat is the estimator of a
%   It is the estimate support of attack vector 
%   y is m*1 vector used in CUSUM
%   s is the spasity
%   A'A is (sigma_z)^(-1) m*m

m=length(y);%get m
r=zeros(m,s+1);%%%%%%%%%%%%%% r(:,t)=r_t-1 %%%%%%%%%%%%%%%
It=[];%I0 emptyC
r(:,1)=y;%r0=y
t=0;%iteration
j=1:m;%index set
i=zeros(1,s);%column index
%Pt=zeros(m,m);
while t<s
    t=t+1;
    [~,i(1,t)]=max(abs(r(:,t)'*A(:,j)));%find i_t
    j(i(1,t))=[];%delete i_t from j
    It=sort([i(1,t),It]);%It={i_t}U{I_t-1}
    AIt=A(:,It);
    Pt=AIt*inv(AIt'*AIt)*AIt';
    r(:,t+1)=(eye(m)-Pt)*y;
end
a_hat=(AIt'*AIt)\AIt'*y;