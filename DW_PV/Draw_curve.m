% Plot PFA-ADD curve
clc; clear; close all;
%%
% fdi_d
add1 = [8.1037 9.0176 10.0262 11.0849 12.2416 13.4444 15.6553 18.8999];
pfa1 = [0.190 0.092 0.047 0.034 0.015 0.01 0.002 0];

%fdi_n
add2 = [90.2316 91.4088 95.2234 96.644 98.9666 108.1858 118.6273 142.6240];
pfa2 = [0.171 0.134 0.087 0.059 0.041 0.015 0.002 0];

%replay
add3 = [69.2341 77.2751 83.5837 94.0831 99.9059 109.32 128.0304 146.2623];
pfa3 = 1e-3*[184 142 104 73 44 28 12 1];

%destable
add4 = [14.914 15.1585 15.5263 15.8086 16.2065 16.5996 16.8353 18.518];
pfa4 = 1e-3*[186 142 88 44 22 16 4 0];

%%
figure
plot(pfa1,add1,'bo-','LineWidth',2);
grid on
hold on
plot(pfa4,add4,'r*-','LineWidth',2);
xlabel('PFA')
ylabel('ADD')
legend('Deterministic FDI','Destabilize')

figure
plot(pfa2,add2,'bo-','LineWidth',2);
grid on
hold on
plot(pfa3,add3,'r*-','LineWidth',2);
xlabel('PFA')
ylabel('ADD')
legend('Noise FDI','Replay')