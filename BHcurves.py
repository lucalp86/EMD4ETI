# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 00:48:05 2022

@author: lucap
"""

import numpy as np
import matplotlib.pyplot as plt

mu0 = 4 * np.pi * 10**-7
H = np.linspace(0, 10**5, 1000)
Bsat = 1.5
muR = 1000

B = mu0 * H + 2 * Bsat / np.pi * np.arctan(mu0 * np.pi / (2 * Bsat) * (muR - 1) * H)

plt.figure()
plt.plot(H, B)
plt.xlabel('H [A/m]')
plt.ylabel('B [T]')

import scipy.io
mat = scipy.io.loadmat('BHcurves_students.mat')

H_Hoganas = mat['H_Hoganas']
B_Hoganas = mat['B_Hoganas']

plt.figure()
plt.plot(H_Hoganas, B_Hoganas)
plt.xlabel('H [A/m]')
plt.ylabel('B [T]')

H_M19 = mat['H_M19']
B_M19 = mat['B_M19']

plt.figure()
plt.plot(H_M19, B_M19)
plt.xlabel('H [A/m]')
plt.ylabel('B [T]')

plt.figure()
plt.plot(H_M19, B_M19, color = 'red', label = 'M19')
plt.plot(H_Hoganas, B_Hoganas, color = 'blue', label = 'Hoganas')
plt.plot(H, B, color = 'green', label = 'atan')
plt.legend(loc = 'best')
plt.xlim([0, 80*10**3])
plt.xlabel('H [A/m]')
plt.ylabel('B [T]')


plt.figure()
plt.plot(H_M19, np.gradient(B_M19[:,0], H_M19[:,0], axis = 0)/mu0, color = 'red', label = 'M19')
plt.plot(H_Hoganas, np.gradient(B_Hoganas[:,0], H_Hoganas[:,0], axis = 0)/mu0, color = 'blue', label = 'Hoganas')
plt.plot(H, np.gradient(B, H, axis = 0)/mu0, color = 'green', label = 'atan')
plt.xlim([0, 20*10**3])
plt.ylim([0, 12.5**3])

# plt.plot(H_Hoganas, B_Hoganas, color = 'blue', label = 'Hoganas')
# plt.plot(H, B, color = 'green', label = 'atan')
# plt.legend(loc = 'best')
# plt.xlabel('H [A/m]')
# plt.xlabel('B [T]')

