% Estimate the state-space model of the PV data
clc; clear;close all;
load data1min.mat

%% Estimate A,B
data1 = iddata(delx',delu',pmuTime);
A = zeros(n,n);
B = zeros(n,m);
C = eye(n);
D = zeros(n,m);
K = zeros(n,n);
x0 = delx(:,1);
Ts = pmuTime;
init_sys1 = idss(A,B,C,D,K,x0,Ts);
init_sys1.Structure.C.Free = false;
init_sys1.Structure.D.Free = false;
init_sys1.Structure.K.Free = false;
sys1 = ssest(data1,init_sys1);
%% Estimate C
data2 = iddata(dely',delu',pmuTime);
C = zeros(p,n);
D = zeros(p,m);
K = zeros(n,p);
init_sys2 = idss(sys1.A,sys1.B,C,D,K,x0,Ts);
init_sys2.Structure.A.Free = false;
init_sys2.Structure.B.Free = false;
init_sys2.Structure.K.Free = false;
sys2 = ssest(data2,init_sys2);
%%
A_e = sys2.A;
B_e = sys2.B;
C_e = sys2.C;
%save est_matrix.mat A_e B_e C_e