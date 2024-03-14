
% data2a.m file was missing from the data folder though it is the file used
% in the example for steady state simulations of the PSTmanual 
% so I, Israel, copied the data from thhe manual to recreate this
% file...not sure if it has the exact content of the original data2a.m file


% bus data format
% bus:
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
% bus_type - 1, swing bus
% - 2, generator bus (PV bus)
% - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max pu
% col15 v_min pu


bus = [...
1   1.03     18.5   7.00  1.61  0.00  0.00  0.00  0.00  1  99.0  -99.0  22.0  1.1  .9; 
2   1.01     8.80   7.00  1.76  0.00  0.00  0.00  0.00  2  5.0   -2.0   22.0  1.1  .9; 
3   0.9781  -6.1    0.00  0.00  0.00  0.00  0.00  0.00  3  0.0   0.0   230.0  1.5  .5; 
4   0.95    -10     0.00  0.00  9.76  1.00  0.00  0.00  3  0.0   0.0   115.0 1.05  .95; 
10  1.0103  12.1    0.00  0.00  0.00  0.00  0.00  0.00  3  0.0   0.0   230.0  1.5  .5; 
11  1.03    -6.8    7.16  1.49  0.00  0.00  0.00  0.00  2  5.0  -2.0    22.0  1.1  .9; 
12  1.01    -16.9   7.00  1.39  0.00  0.00  0.00  0.00  2  5.0  -2.0    22.0  1.1  .9; 
13  0.9899  -31.8   0.00  0.00  0.00  0.00  0.00  0.00  3  0.0   0.0   230.0  1.5  .5; 
14  0.95    -38     0.00  0.00  17.67 1.00  0.00  0.00  3  0.0   0.0   115.0 1.05  .95;  
20  0.9876   2.1    0.00  0.00  0.00  0.00  0.00  0.00  3  0.0   0.0   230.0  1.5  .5; 
101 1.05    -19.3   0.00  8.00  0.00  0.00  0.00  0.00  2  99.0  -99.0 230.0  1.5  .5; 
110 1.0125  -13.4   0.00  0.00  0.00  0.00  0.00  0.00  3  0.0   0.0   230.0  1.5  .5; 
120 0.9938  -23.6   0.00  0.00  0.00  0.00  0.00  0.00  3  0.0   0.0   230.0  1.5  .5 ]; 
 
% line data format
% line:
% col1 from bus
% col2 to bus
% col3 resistance(pu)
% col4 reactance(pu)
% col5 line charging(pu)
% col6 tap ratio
% col7 tap phase
% col8 tapmax
% col9 tapmin
% col10 tapsize

line = [...
1   10  0.0     0.0167   0.00    1.0  0. 0.  0.  0.; 
2   20  0.0     0.0167   0.00    1.0  0. 0.  0.  0.; 
3    4  0.0     0.005     0.00   1.0  0. 1.2 0.8 0.05; 
3   20  0.001   0.0100   0.0175  1.0  0. 0.  0.  0.; 
3   101 0.011   0.110    0.1925  1.0  0. 0.  0.  0.; 
3   101 0.011   0.110    0.1925  1.0  0. 0.  0.  0.; 
10  20  0.0025  0.025    0.0437  1.0  0. 0.  0.  0.; 
11  110 0.0     0.0167   0.0     1.0  0. 0.  0.  0.; 
12  120 0.0     0.0167   0.0     1.0  0. 0.  0.  0.; 
13   14 0.0     0.005    0.00    1.0  0. 1.2 0.8 0.05; 
13  101 0.011   0.11     0.1925  1.0  0. 0.  0.  0.; 
13  101 0.011   0.11     0.1925  1.0  0. 0.  0.  0.; 
13  120 0.001   0.01     0.0175  1.0  0. 0.  0.  0.; 
110 120 0.0025  0.025    0.0437  1.0  0. 0.  0.  0.]; 