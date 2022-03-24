# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 21:49:23 2022

@author: lucap
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.close('all')

def pendulum(alfa, time, G, D, L, M):

    alfa1 = alfa[0]
    alfa2 = alfa[1]

    dalfa1_dt = alfa2
    dalfa2_dt = - (D / M) * alfa2 - (G / L) * np.sin(alfa1)
    
    dalfa_dt = [dalfa1_dt, dalfa2_dt]

    return dalfa_dt

# Parameters
G = 9.81     # [m / s**2]
M = 1        # [kg]
L = 1        # [m]
D = 0.3        # [N / m]

# Initial conditions
alfa1_0 = 90 * np.pi / 180
alfa2_0 = 0 * np.pi / 180

alfa0 = [alfa1_0, alfa2_0] 

# Simulation time
time_final = 20
dt = 0.05
time = np.arange(0.0, time_final, dt)

#  Compute solution
alfa = odeint(pendulum, alfa0, time, args = (G, D, L, M))

x = L * np.sin(alfa[:, 0])
y = -L * np.cos(alfa[:, 0])

fig = plt.figure()
plt.subplot(2, 1, 1)
plt.plot(time, alfa[:, 0] * 180 / np.pi, 'r')
plt.xlabel('Time [s]')
plt.ylabel('Angle [deg]')

plt.subplot(2, 1, 2)
plt.plot(time, alfa[:, 1] * 180 / np.pi, 'b')
plt.xlabel('Time [s]')
plt.ylabel('Angular speed [deg/s]')

fig = plt.figure()
plt.subplot(1, 1, 1)
plt.plot(alfa[:, 0] * 180 / np.pi, alfa[:, 1] * 180 / np.pi, 'b')
plt.ylabel('Angular speed [deg/s]')
plt.xlabel('Angle [deg]')

kplot = 1.2

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, autoscale_on = False, xlim = (-L * kplot, L * kplot), ylim = (-L * kplot, L * kplot))
ax.set_aspect('equal')


line, = ax.plot([], [], 'o-', lw = 2, color = 'red')
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def initialisation():
    
    line.set_data([], [])
    
    time_text.set_text('')
    return ax, time_text

def animate(i):
    xt = [0, x[i]]
    yt = [0, y[i]]
 
    line.set_data(xt, yt)
    time_text.set_text(time_template % (i*dt))
    return ax, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(time)),
                              interval = 25, blit = True, init_func = initialisation)

# ani.save('double_pendulum.mp4', fps=15)
