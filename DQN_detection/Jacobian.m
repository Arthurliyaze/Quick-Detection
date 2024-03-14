function [H] = Jacobian(V,fbus,nvi,npi,nqi,npf,nqf,nbus,G,del,B,qi,tbus,bpq,ppi)

%% initialization of Jacobian matrix H
H11 = zeros(nvi,nbus-1);
H22 = zeros(npi,nbus);
H31 = zeros(nqi,nbus-1);
H32 = zeros(nqi,nbus);
H41 = zeros(npf,nbus-1);
H42 = zeros(npf,nbus);
H51 = zeros(nqf,nbus-1);
H52 = zeros(nqf,nbus);
% H11 - Derivative of V with respect to angles.. All Zeros
% H11 = zeros(nvi,nbus-1);



% H12 - Derivative of V with respect to V..
H12 = zeros(nvi,nbus);
for k = 1:nvi
    for n = 1:nbus
        if n == k
            H12(k,n) = 1;
        end
    end
end

% H21 - Derivative of Real Power Injections with Angles..
H21 = zeros(npi,nbus-1);
for i = 1:npi
    %m = fbus(ppi(i));
    m = ppi(i);
    for k = 1:(nbus-1)
        if k+1 == m
            for n = 1:nbus
                H21(i,k) = H21(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
            end
            H21(i,k) = H21(i,k) - V(m)^2*B(m,m);
        else
            H21(i,k) = V(m)* V(k+1)*(G(m,k+1)*sin(del(m)-del(k+1)) - B(m,k+1)*cos(del(m)-del(k+1)));
        end
    end
end

% H22 - Derivative of Real Power Injections with V..
% H22 = zeros(npi,nbus);
for i = 1:npi
    %m = fbus(ppi(i));
    m = ppi(i);
    for k = 1:(nbus)
        if k == m
            for n = 1:nbus
                H22(i,k) = H22(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
            H22(i,k) = H22(i,k) + V(m)*G(m,m);
        else
            H22(i,k) = V(m)*(G(m,k)*cos(del(m)-del(k)) + B(m,k)*sin(del(m)-del(k)));
        end
    end
end

% H31 - Derivative of Reactive Power Injections with Angles..
% H31 = zeros(nqi,nbus-1);
for i = 1:nqi
    %m = fbus(qi(i));
    m = qi(i);
    for k = 1:(nbus-1)
        if k+1 == m
            for n = 1:nbus
                H31(i,k) = H31(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
            H31(i,k) = H31(i,k) - V(m)^2*G(m,m);
        else
            H31(i,k) = V(m)* V(k+1)*(-G(m,k+1)*cos(del(m)-del(k+1)) - B(m,k+1)*sin(del(m)-del(k+1)));
        end
    end
end

% H32 - Derivative of Reactive Power Injections with V..
% H32 = zeros(nqi,nbus);
for i = 1:nqi
    %m = fbus(qi(i));
    m = qi(i);
    for k = 1:(nbus)
        if k == m
            for n = 1:nbus
                H32(i,k) = H32(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
            H32(i,k) = H32(i,k) - V(m)*B(m,m);
        else
            H32(i,k) = V(m)*(G(m,k)*sin(del(m)-del(k)) - B(m,k)*cos(del(m)-del(k)));
        end
    end
end

% H41 - Derivative of Real Power Flows with Angles..
% H41 = zeros(npf,nbus-1);
for i = 1:npf
    %m = fbus(pf(i));
    %n = tbus(pf(i));
    m = fbus(i);
    n = tbus(i);
    for k = 1:(nbus-1)
        if k+1 == m
            H41(i,k) = V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
        else if k+1 == n
                H41(i,k) = -V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
            else
                H41(i,k) = 0;
            end
        end
  
    end
    
    
end

% H42 - Derivative of Real Power Flows with V..
% H42 = zeros(npf,nbus);
for i = 1:npf
    %m = fbus(pf(i));
    %n = tbus(pf(i));
    m = fbus(i);
    n = tbus(i);
    for k = 1:nbus
        if k == m
            H42(i,k) = -V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n))) - 2*G(m,n)*V(m);
        else if k == n
                H42(i,k) = -V(m)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            else
                H42(i,k) = 0;
            end
        end

    end
end

% H51 - Derivative of Reactive Power Flows with Angles..
% H51 = zeros(nqf,nbus-1);
for i = 1:nqf
    %m = fbus(qf(i));
    %n = tbus(qf(i));
    m = fbus(i);
    n = tbus(i);
    for k = 1:(nbus-1)
        if k+1 == m
            H51(i,k) = -V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
        else if k+1 == n
                H51(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            else
                H51(i,k) = 0;
            end
        end

    end
end

% H52 - Derivative of Reactive Power Flows with V..
% H52 = zeros(nqf,nbus);
for i = 1:nqf
    %m = fbus(qf(i));
    %n = tbus(qf(i));
    m = fbus(i);
    n = tbus(i);
    for k = 1:nbus
        if k == m
            H52(i,k) = -V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n))) - 2*V(m)*(-B(m,n)+ bpq(m,n));
        else if k == n
                H52(i,k) = -V(m)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
            else
                H52(i,k) = 0;
            end
        end

    end
end

%% Measurement Jacobian, H..--------------------------------------------
H = [H11 H12; H21 H22; H31 H32; H41 H42; H51 H52];

end