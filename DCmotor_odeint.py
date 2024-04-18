#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 11:02:56 2023

@author: luca
"""

import sys
sys.modules[__name__].__dict__.clear()

import time
startTime0 = time.time()

import pandas as pd
import numpy as np
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from scipy import signal
import scipy
from cmcrameri import cm
import sympy as sym
sym.init_printing(use_latex=True)

Folder_dir = os.getcwd()
Folder_data = Folder_dir
Folder_dir_save = Folder_dir + '\\Images\\' 


# mpl.rcParamsDefault     # this command allows to get the default settings of rc Params
# mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.size"] = 12 # default = 'medium'

mpl.rcParams["axes.labelsize"] = 12 # default = 'medium'
mpl.rcParams["axes.labelpad"] = 1 # default = 4

mpl.rcParams["xtick.labelsize"] = 12 # default = 'medium'
mpl.rcParams["ytick.labelsize"] = 12 # default = 'medium'
mpl.rcParams["xtick.major.pad"] = 1 # default = 3.4
mpl.rcParams["ytick.major.pad"] = 1 # default = 3.4

mpl.rcParams["axes.linewidth"] = 0.5 # default = 0.8

#from matplotlib import rc

mpl.rcParams["text.usetex"] = True

plt.close("all")
FigSize = (7, 9)

R = 0.1
C = 0.00001
L = 0.002
K = 0.01
G = 0.3
M = 0.01
J = 0.0001

KP = 1
KI = 100

f = [1, 1] 

def InputDef(time):
    # u = [v, Wm]
    i_ref0 = 1
    WM_0 = 1
    
    if time > 0.5 and time < 1.5:
            u = np.array([i_ref0, 0 * WM_0]) 
    elif time >= 1.5 and time < 20:
            u = np.array([i_ref0, WM_0]) 
    elif time >= 25 and time < 60:
            u = np.array([2 * i_ref0, 0 * WM_0]) 
    elif time >= 60:
            u = np.array([i_ref0, WM_0]) 
    else:
            u = np.array([0 * i_ref0, 0 * WM_0]) 

    # ck = R_tot

    return u

    
def PMDC_mot(f, time, R, M, L, J, K, G):
    # y = [i, wm]

    # A * dy + B * y = C * u

    # u = [i_ref, W_M]
    u = InputDef(time)

    A = np.array([[L, 0],
                [0, J]])
    
    B = np.array([[R, M],
                [-K, G]])
    
    D = np.array([[1, 0],
                [0, -1]])
    
    invF = np.linalg.inv(A)
    # dy + inv(A) * B * y = inv(A) * C * u

    Anew = invF@B
    Bnew = invF@D
    # dy + Anew * y = Bnew * u

    # dy = - Anew * y + Bnew * u
    dfdt = - Anew @ f + Bnew @ u

    return dfdt

# Initial conditions
# y0 = [i0, wm0]
w0 = 0
i0 = 0

f0 = [i0, w0] 

# Simulation time
time_final = 10
dt = 0.001
time = np.arange(0.0, time_final, dt)

from scipy.integrate import odeint
#  Compute solution
sol = odeint(PMDC_mot, f0, time, args = (R, L, K, G, M, J))

inp = np.zeros([time.shape[0], 2])
for idx_t, val_t in enumerate(time):
    inp[idx_t, :] = InputDef(val_t)
    
plt.figure()
plt.plot(time, sol[:,0])
plt.figure()
plt.plot(time, sol[:,1])
# plt.figure()
# plt.plot(time, sol[:,2])

plt.figure()
plt.plot(time, inp[:,0])
plt.figure()
plt.plot(time, inp[:,1])

