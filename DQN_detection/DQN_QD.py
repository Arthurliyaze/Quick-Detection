# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 11:17:49 2020
Project: Quick detection of FDI by DQN, discrete state
@author: yazeli
Version: 3 
"""
#%% Importing Libraries and Dependencies
from silence_tensorflow import silence_tensorflow
silence_tensorflow()

import numpy as np
import os
clear = lambda: os.system('cls')
import matlab.engine
import time

from gym import Env
from gym.spaces import Discrete,MultiDiscrete
#%%
from stable_baselines.deepq.policies import FeedForwardPolicy
from stable_baselines.common.vec_env import DummyVecEnv
from stable_baselines import DQN
#%%
eng = matlab.engine.start_matlab()

#%% Parameters
N = 13    # number of buses
n = 2*N-1 # length of states
m = 55    # length pf measurement
M = 4     # slide window length
T = 100   # sample of measurement

phi = 1 #  0.01 0.02 0.05 0.1 0.2 0.5 1


#%% Define environment
class DetectEnv(Env):
  """
  A customized environment for training and testing
  """
  metadata = {'render.modes': ['human']}
  
  def __init__(self, train):
      """
      Environment Parameters
      dimensions of action space (1): declare an attack or not
      dimensions of state space (4): slide window of Rao statistic
      """      
      self.action_space = Discrete(2)
      self.observation_space =  MultiDiscrete([4,4,4,4])
      self.ep_length = T
      self.tau_detect = 0
      self.train = train
      self.reset()
      
      
      
  def reset(self):
      self.current_step = 1
      self.x_predict = np.concatenate((np.zeros((N-1,1)),np.ones((N,1))),axis=None).reshape((1,n))
      self.Mk = np.zeros((n,n))
      self.a = np.zeros((1,n))
      self.b = np.zeros((1,n))
      self.Tw = np.zeros(T)
      self.dTw = np.zeros(T)
      self.attack_num = np.random.randint(low = 1, high = n, size = 1)
      self.idx_inj = np.random.choice(n, size = self.attack_num, replace = False)+1
      self.attack_mag = np.random.randint(low = -2000, high = 2000, size = self.attack_num)
      if not self.train: # testing
          self.tau = np.random.randint(low=T/2,high=T)
      else: # training          
          self.tau = np.random.randint(low=4,high=T/2)
      self.state = np.zeros(M)
      return self.state
  
  def step(self, action):
      
          
      # Convert np datatype to matlab
      m_idx_inj = matlab.int32(self.idx_inj.tolist())
      m_attack_mag = matlab.double(self.attack_mag.tolist())
      m_x_predict = matlab.double(self.x_predict.tolist())
      m_M = matlab.double(self.Mk.tolist())
      m_a = matlab.double(self.a.tolist())
      m_b = matlab.double(self.b.tolist())
      
      
      m_t,m_x_predict,m_M,m_a,m_b = eng.step_dos(self.current_step, self.tau, m_idx_inj, m_attack_mag, m_x_predict, m_M, m_a, m_b, nargout = 5)
      self.Tw[self.current_step-1] = np.float32(m_t)
      self.x_predict = np.float32(m_x_predict)
      self.Mk = np.float32(m_M)
      self.a = np.float32(m_a)
      self.b = np.float32(m_b)
      
      # Discretize observations
      if self.Tw[self.current_step-1] > 32:
          self.dTw[self.current_step-1] = self.dTw[self.current_step-1]+1
      if self.Tw[self.current_step-1] > 64:
          self.dTw[self.current_step-1] = self.dTw[self.current_step-1]+1
      if self.Tw[self.current_step-1] > 128:
          self.dTw[self.current_step-1] = self.dTw[self.current_step-1]+1
      
      ob = np.zeros(M)
      for i in range(M):
           if self.current_step-M+i >= 0:
               ob[i] = self.dTw[self.current_step-M+i]
      
      if self.current_step <= self.tau:
          reward = 1/self.tau
      else:
          reward = -phi
              
      if action == 1:
          reward = 0
          self.tau_detect = self.current_step
          self.current_step = self.ep_length
      
      # Transition
      self.current_step += 1
      done = self.current_step > self.ep_length
      if done:
          self.reset
      else:
          self.state = ob
      return self.state, reward, done, {}

  def render(self, mode='human', close=False):
      pass

#%% Define policy
# Custom MLP policy of two layers of size [32,32]
class CustomDQNPolicy(FeedForwardPolicy):
    def __init__(self, *args, **kwargs):
        super(CustomDQNPolicy, self).__init__(*args, **kwargs,
                                           layers=[32, 32],
                                           layer_norm=True,
                                           feature_extraction="mlp")
#%% Define the model
log_dir = "./dqn_qd/"
# command for tensorboard:
# tensorboard --logdir="./dqn_qd"
#Define the training environment and model
# train = DummyVecEnv([lambda: DetectEnv(train=True)])
# model = DQN(CustomDQNPolicy, train, gamma=1, verbose=0, tensorboard_log=log_dir)
# #%% Train the model
# model.learn(total_timesteps= 100000, tb_log_name='dDQN')
# model.save("ddqn_m2phi0.1")
#%% Testing
point = 1
trial = 100
add = np.zeros(point)
pfa = np.zeros(point)
tic = time.clock()
for i in range(point):
    dd = np.ones(trial)*np.nan
    fa = 0
    actions = np.zeros(100)
    obsers = np.zeros((M,100))
    ep_reward = 0
    tau = np.zeros(trial)
    tau_detect = np.zeros(trial)
    for j in range(trial):
        clear()
        if j%10 == 0:
            print('trial = ',j)
        test = DummyVecEnv([lambda: DetectEnv(train=False)])
        model = DQN.load("ddqn_m2phi1", env=test, policy=CustomDQNPolicy)
        obs = test.reset()
        dones = False
        k = 0
        tau[j] = test.envs[0].tau
        while not dones:
            action, _states = model.predict(obs)
            obs, rewards, dones, info = test.step(action)
            actions[k] = action
            obsers[:,k] = obs
            ep_reward = ep_reward+rewards
            k = k+1
            test.render()
        tau_detect[j] = test.envs[0].tau_detect    
        if tau_detect[j] >= tau[j]:
              dd[j] = tau_detect[j]- tau[j]
        else:
              fa = fa+1
    add[i] = np.nanmean(dd)
    pfa[i] = fa/trial
toc = time.clock()
minutes = (toc-tic)/60