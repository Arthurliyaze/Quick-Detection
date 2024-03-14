 function [h] = measurement_fcn_h(V,fbus,v,npi,nqi,npf,nqf,nbus,G,del,B,qi,tbus,bpq,ppi )

%% initialization of h(x)
h1 = v;
h2 = zeros(npi,1);
h3 = zeros(nqi,1);
h4 = zeros(npf,1);
h5 = zeros(nqf,1);
%%
% P = 0;
for i = 1:npi
    %     m = fbus(ppi(i));
    m = ppi(i);
    for k = 1:nbus
        Pik = V(m)*V(k)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
        h2(i) = h2(i) + Pik;
        
    end
end
%%
for i = 1:nqi
    m = qi(i);
    for k = 1:nbus
        h3(i) = h3(i) + V(m)*V(k)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
    end
end

for i = 1:npf

    m = fbus(i);
    n = tbus(i);
    
    h4(i) = -V(m)^2*G(m,n) - V(m)*V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
end

for i = 1:nqf

    m = fbus(i);
    n = tbus(i);
    h5(i) = -V(m)^2*(-B(m,n)+bpq(m,n)) - V(m)*V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
 
end

h = [h1; h2; h3; h4; h5];

