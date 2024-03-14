# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 01:56:28 2020
This code tests the step function in matlab
@author: yazeli
"""

import numpy as np
import matlab.engine
#%% Parameters
N = 13    # number of buses
n = 2*N-1 # length of states
m = 55    # length pf measurement
M = 4     # slide window length
T = 500   # sample of measurement
idx_inj = np.array([4,5])   # attacked bus index
attack_mag = np.array([-1500,1000])
x_predict = np.concatenate((np.zeros((N-1,1)),np.ones((N,1))),axis=None).reshape((1,n))
Mk = np.zeros((n,n))
a = np.zeros((1,n))
b = np.zeros((1,n))
current_step = 0
tau = 250
n_sample = 10
res = np.zeros((m,n_sample))
#%% Convert np datatype to matlab
m_idx_inj = matlab.double(idx_inj.tolist())
m_attack_mag = matlab.double(attack_mag.tolist())

#%%
eng = matlab.engine.start_matlab()
for i in range(n_sample):
    m_x_predict = matlab.double(x_predict.tolist())
    m_M = matlab.double(Mk.tolist())
    m_a = matlab.double(a.tolist())
    m_b = matlab.double(b.tolist())
       
    m_res,m_x_predict,m_M,m_a,m_b = eng.step(i+1, tau, m_idx_inj, m_attack_mag, m_x_predict, m_M, m_a,m_b,nargout = 5)
    res[:,i] = np.float32(m_res)
    x_predict = np.float32(m_x_predict)
    Mk = np.float32(m_M)
    a = np.float32(m_a)
    b = np.float32(m_b)