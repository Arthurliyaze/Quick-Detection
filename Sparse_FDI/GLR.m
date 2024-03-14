function [index_q] = GLR(z,Sigma_zr,s,l)
%GLR
%   z is the observed measurement matrix m*500
%   Sigma_zr is inv of Sigma_z
%   s is the sparsity assume we now it
%   l is the time

m=length(z);
pattern=nchoosek(1:m,s);
[Qs,~]=size(pattern);
%for k=1:l
k=1;
    sumz=sum(z(:,k:l),2);%sum z_i from k->l
    %for s=1:m
    for q=1:Qs%search all pattern
        aq=pattern(q,:);
        Lambda_q=Sigma_zr(:,aq);
        Phi_q=Sigma_zr(aq,aq);
        eta(q)=0.5*sumz'*Lambda_q*inv(Phi_q)*Lambda_q'*sumz/(l-k+1);
    end
[~,Nq]=max(eta); 
index_q=pattern(Nq);
end

