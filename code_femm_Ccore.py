# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 10:32:01 2019

@author: Luca Papini
"""

# ==========================================================================
# Copyright (C) 2018 Dr. Luca Papini
#
# This file is part of Pyfemm, Python code for design and analysis of 
# electrical machines
#
# Pyfemm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pyfemm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see 
#
#             http://www.gnu.org/licenses/
#
# =============================================================================
# ================================ PYFEMM =====================================
# ================================ DESIGN =====================================
# ********************** Active Magnetic Bearing (AMB) ************************
# ************************* Inner Rotor, w/o Sleeve  **************************
# ********************* Single Layer Distributed winding **********************
# =============================================================================
#
# =============================================================================
# Program:      Py_MDfemm.m
# Sub-program:  Py_MDfemm_AMB.m
# Author:       Luca Papini (lpapini)
# Date:         04/11/2019
# Version:      0.1.1
# Requirement:  Pyton (femm,numpy,matplotlib libraries), FEMM, 
#               Material Library, Winding Modeller
# References:   - Py-FEMM Manual
#
# Revision History:
#      Date       Version    Author      Description
#  - 04/11/2019:  0.1.1      lpapini     Main code for femm in Python
#                                        -Added circular shape edge slot
#                                        -Added timing
#                                        -Skewing modelling

#
# Missing parts:
#    Status     Description                               Add in Version
#               Magnetization Generalization 
#               (generic input span/magnetization/sub-element)                            
#      v        -  Radial                                             0.1.5
#      v        -  Parallel                                           0.1.5
#      x        - Halbach
#      x        - Halbach generalized
#               Different magnetization 
#      v        - Radial                                             0.1.2
#      v        - Parallel                                           0.1.2
#      x        - Halbach
#      x        - Halbach generalized
#               Different Magnet shapes
#      x        - bread loaft
#      x        - fully shaped
#               Rotor structure (RotorType variable)
#      v        - Inset PM                                           0.1.2
#      v        - CP                                                 0.1.5
#      v        - CP/InsetPM                                         0.1.5
#      x        - CP/InsetPM with daftail rotor 
#      x        - Syhncronous Reluctance
#      x        Generalize winding structure
#      x        Winding Connection (weastone bridge, Y/D)                             
#      x        Add non-linear materials
#      v        Skewing                                              0.1.4
#
# Bugs:
# Status        Description                               Fixed in Version
#
# To Test:
# Status        Description                             Tested with Version
#      x        PM width=1                                    
#      x        Slot/pole combination                         
#        - v     - Integer ~=1                                       0.1.1
#        - x     - Integer =1                                 
#        - v     - Fractional ~=1                                    0.1.4
#        - x     - Fractional =1                              
#      x        Rotor topologies
#        - v     - Surface Mounted (RotorType=0)                     0.1.1
#        - v     - Inset (RotorType=1)                               0.1.4
#        - v     - Consequent Pole (RotorType=2/3)                   0.1.5
#
# Legend:       x-> open (not done, to do)     v-> closed (done, implemented)
# =============================================================================
import time
startTime0 = time.time()

import femm 
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from matplotlib import rc
import pandas as pd

mpl.rcParams["text.usetex"] = True

def draw_part(RZ, Pos):
    coord = []
    for key in RZ.keys():    
        coord.append(RZ[key])
        for cd in np.arange(2):
            coord[-1][cd] += Pos[cd]
        
    
    for pt in np.arange(len(RZ)):  
        if pt == np.arange(len(RZ))[-1]:
            pt_nxt = np.arange(len(RZ))[0]
        else:
            pt_nxt = pt + 1
            
        femm.mi_drawline(coord[pt][0], coord[pt][1], coord[pt_nxt][0], coord[pt_nxt][1])
    
def create_component(const, name, mesh_size = 1, group = 0, dir_mag = 0):
     
    femm.mi_addblocklabel(const[0], const[1])
    femm.mi_selectlabel(const[0], const[1])
    femm.mi_setblockprop(name, 0, mesh_size, 0, dir_mag, 0, 0)
    #femm.mi_setblockprop(name, 0, mesh_size, 0, 'theta', 0, 0)
    femm.mi_setgroup(group)
    femm.mi_clearselected()

def build_part(RZ, Pos, const, name, group = 0, mesh_size = 1):

    # ----- Draw Component 
    draw_part(RZ, Pos)
    
    # ----- Create Component 
    create_component(const, name, mesh_size, group)
#    create_component(const, name = name, mesh_size = mesh_size, group = group):


# def build_part(RZ, Pos, const, name, group = 0, mesh_size = 1):

#     # ----- Draw Component 
#     coord = []
#     for key in RZ.keys():    
#         coord.append(RZ[key])
#         for cd in np.arange(2):
#             coord[-1][cd] += Pos[cd]
        
    
#     for pt in np.arange(len(RZ)):  
#         if pt == np.arange(len(RZ))[-1]:
#             pt_nxt = np.arange(len(RZ))[0]
#         else:
#             pt_nxt = pt + 1
            
#         femm.mi_drawline(coord[pt][0], coord[pt][1], coord[pt_nxt][0], coord[pt_nxt][1])
    
#     # ----- Create Component 
#     femm.mi_addblocklabel(const[0], const[1])
#     femm.mi_selectlabel(const[0], const[1])
#     femm.mi_setblockprop(name, 0, mesh_size, 0, 0, 0, 0)
#     femm.mi_setgroup(group)
#     femm.mi_clearselected()
    
# mpl.rcParamsDefault     # this command allows to get the default settings of rc Params
mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.size"] = 12 # default = 'medium'

mpl.rcParams["axes.labelsize"] = 12 # default = 'medium'
mpl.rcParams["axes.labelpad"] = 1 # default = 4

mpl.rcParams["xtick.labelsize"] = 12 # default = 'medium'
mpl.rcParams["ytick.labelsize"] = 12 # default = 'medium'
mpl.rcParams["xtick.major.pad"] = 1 # default = 3.4
mpl.rcParams["ytick.major.pad"] = 1 # default = 3.4

mpl.rcParams["axes.linewidth"] = 0.5 # default = 0.8

plt.close("all")


Folder_dir = os.getcwd()

Colors = np.array([[1,0,0],[0,0,1],[0,0.6,0],[1,0.5,0],[0.05,0.52,0.78],[0.17,0.5,0.34],[0.85,0.16,0],[0.08,0.17,0.55],[0.1,0.31,0.21],[0.75,0,0.75],[0,0.8,0.8],[0.2,0.8,0]])
Colors = np.tile(Colors,[3,1])

pi = np.pi

# ----- Open femm application
femm.openfemm() 

# ----- Create new document
#  0 for a magnetics problem, 
#  1 for an electrostatics problem, 
#  2 for a heat flow problem, 
#  3 for a current flow problem

femm.newdocument(0);

# femm.hideconsole() 
femm.main_resize(600,600) 

mu0 = 4*pi*10**-7          # Vacuum magnetic permeability [H/m]

Res20_Cu = 1.68*10**-8        # Resistivity copper [ohm*m] @ 20 [C]
alphaT_Cu = 0.004041          # Temperature-resistivity coefficient copper [1/C]

Res20_GF = 1.2*10**-8         # Resistivity glass fiber [ohm*m] @ 20 [C]
alphaT_GF = 0.0004            # Temperature-resistivity coefficient glass fiber [1/C]

active_length = 20

# ----- Problem definition
#model_type = 0     # 0 -> planar      1 -> axialsymmetric
#femm.mi_probdef(0 ,'millimeters','axi' ,10**-9, active_length ,30)
femm.mi_probdef(0,'millimeters','planar' ,10**-9, active_length ,30)

femm.showpointprops
femm.main_restore

# =========================================================================
# ================ Create Model Material (linear)
# =========================================================================

# ------ Air gap
Mu_X_AG = 1              # Air gap fluid X-Relative magnetic permeability
Mu_Y_AG = Mu_X_AG        # Air gap fluid Y-Relative magnetic permeability
Sigma_AG = 0             # Conductivity air gap fluid in [MS/m]

# ------ Stator Conductor 
Temp_ConductStat = 180                                           # Operative temperature Stator Conductor  
Rho_ConductStat = Res20_Cu * (1 + alphaT_Cu * (Temp_ConductStat - 20))   # Resistivity Conductor stator [ohm*m]

Mu_X_CuStat = 1                                    # Stator Copper X-Relative magnetic permeability
Mu_Y_CuStat = Mu_X_CuStat                          # Stator Copper Y-Relative magnetic permeability
Sigma_CuStat = 1/Rho_ConductStat * 10**-6             # Conductivity Stator Copper in [MS/m]

# ------ Stator Iron
NL_Stat = 0                 # 0 -> Linear     1 -> Non Linear
Mu_X_FeStat = 5000         # Stator iron X-Relative magnetic permeability
Mu_Y_FeStat = Mu_X_FeStat    # Stator iron Y-Relative magnetic permeability
#Mu_X_FeStat=float('Inf')
Sigma_FeStat = 0             # Conductivity Stator iron in [MS/m]

# ------ Rotor Iron
NL_Rot = 0                  # 0 -> Linear     1 -> Non Linear
# Mu_X_FeRot=Inf           # Rotor iron X-Relative magnetic permeability
Mu_X_FeRot = Mu_X_FeStat           # Rotor iron X-Relative magnetic permeability
# Mu_FeRot=Inf does not work (yet) with PMs
Mu_Y_FeRot = Mu_X_FeRot      # Rotor iron Y-Relative magnetic permeability
Sigma_FeRot = 0              # Conductivity Rotor iron in [MS/m]

# ------ Shaft
Mu_X_Shaft = 1     #10**10    # Shaft X-Relative magnetic permeability
Mu_Y_Shaft = Mu_X_Shaft       # Shaft Y-Relative magnetic permeability
Sigma_Shaft = 0               # Conductivity Shaft in [MS/m]

# ------ Permanent Magnet
Mu_X_PM = 1.05         # Stator iron X-Relative magnetic permeability
Mu_Y_PM = Mu_X_PM    # Stator iron Y-Relative magnetic permeability
Sigma_PM = 0             # Conductivity Stator iron in [MS/m]
BR_PM = 1.4
Hc_PM = - BR_PM / Mu_X_PM / mu0

## ------ Create FEMM material
femm.mi_addmaterial('Air', Mu_X_AG, Mu_Y_AG, 0, 0, Sigma_AG, 0, 0, 1, 0, 0, 0)

femm.mi_addmaterial('Stator_Coil', Mu_X_CuStat, Mu_Y_CuStat, 0, 0, Sigma_CuStat, 0, 0, 1, 0, 0, 0)

femm.mi_addmaterial('Stator_Fe', Mu_X_FeStat, Mu_Y_FeStat, 0, 0, Sigma_FeStat, 0, 0, 1, 0, 0, 0)

femm.mi_addmaterial('Mover_Fe', Mu_X_FeRot, Mu_Y_FeRot, 0, 0, Sigma_FeRot, 0, 0, 1, 0, 0, 0)

#femm.mi_addmaterial('Mover_Fe', Mu_X_Shaft, Mu_Y_Shaft, 0, 0, Sigma_Shaft, 0, 0, 1, 0, 0, 0)

femm.mi_getmaterial('M-15 Steel')
femm.mi_getmaterial('416 Stainless Steel')

#femm.mi_addmaterial('PM', Mu_X_PM, Mu_Y_PM, Hc_PM, 0, Sigma_PM, 0, 0, 1, 0, 0, 0)

# femm.mi_addmaterial('JmagSteel', Mu_X_FeRot, Mu_Y_FeRot, 0, 0, Sigma_FeRot, 0, 0, 1, 0, 0, 0)
# data_filename = 'B_H_JmagDataLast2.csv'
# Jmag_steel = pd.read_csv(Folder_dir + '\\Materiali\\' + data_filename).to_numpy()
# for pt in range(len(Jmag_steel)):
#     femm.mi_addbhpoint('JmagSteel', Jmag_steel[pt][1] ,Jmag_steel[pt][0])

# femm.mi_addbhpoint('JmagSteel', Jmag_steel[:, 1] ,Jmag_steel[:, 0])

# ====== 
# plt.figure()
# plt.plot(Jmag_steel[:, 0] ,Jmag_steel[:, 1])
# ====== 


# plt.figure()
# plt.plot(np.gradient(Jmag_steel[:, 0]))

# plt.figure()
# plt.plot(np.gradient(Jmag_steel[:, 1]))

# femm.mi_addmaterial('AISI_C12_L14', Mu_X_FeRot, Mu_Y_FeRot, 0, 0, Sigma_FeRot, 0, 0, 1, 0, 0, 0)
# data_filename = 'B_H_AISIC12L14_1.csv'
# AISI = pd.read_csv(Folder_dir + '\\Materiali\\' + data_filename).to_numpy()

# for pt in range(len(Jmag_steel)):
#     femm.mi_addbhpoint('AISI_C12_L14', AISI[pt][1] ,AISI[pt][0])


# =========================================================================
# ================ Geometry Input Data 
# =========================================================================
NonLinear = 1
active_length = 100

# ---- Mover dimensions
eps_mover = 5
width_mover = 20

# ---- Stator & Coil dimensions
air_gap = 0.5
tooth_width = width_mover / 4
opening_width = width_mover - 2 * tooth_width
tooth_heigth = tooth_width

coil_heigth = tooth_heigth/3*2

# -------- Coil and Current settings
turns = 110
current = np.linspace(0, 20, 11)
# current = [10]

# -------- Lagrangian variable settings
chi = np.linspace(0.4, -2, 11)
# chi = [0]

# -------- Mesh settings
SizeMesh_Air = 0.5
SizeMesh_Mover = 1
SizeMesh_Stator = 1
SizeMesh_Coil = 2
SizeMesh_Air_out = 8

# ------- Post Processing
Npoint_Bsampling = 100
PositionY = 0

# MaterialStator = 'Stator_Fe'
# MaterialMover = 'Mover_Fe'

MaterialStator = 'M-15 Steel'
MaterialMover = 'M-15 Steel'

# =========================================================================
# ================ Design Geometry 
# =========================================================================

# ===== Group definition
# -----   0: stator
# -----   1: motion elements
# -----   2: winding

# ********************* Mover
femm.mi_drawline(-width_mover/2, 0, +width_mover/2, 0)
femm.mi_drawline(-width_mover/2, -eps_mover, +width_mover/2, -eps_mover)
femm.mi_drawline(+width_mover/2, 0, +width_mover/2, -eps_mover)
femm.mi_drawline(-width_mover/2, 0, -width_mover/2, -eps_mover)

# ----- Create mover region
femm.mi_addblocklabel(0,-eps_mover/2)
femm.mi_selectlabel(0,-eps_mover/2)
femm.mi_setblockprop(MaterialMover, 0, SizeMesh_Mover, 0, 0, 0, 0)
femm.mi_setgroup(1)
femm.mi_clearselected
    
femm.mi_selectsegment(0, 0)
femm.mi_selectsegment(0, -eps_mover)
femm.mi_selectsegment(eps_mover/2, -eps_mover/2)
femm.mi_selectsegment(-eps_mover/2, -eps_mover/2)
femm.mi_setgroup(1)
femm.mi_clearselected

# ********************* Stator Core

femm.mi_drawline(-width_mover/2, +air_gap, (-width_mover/2+tooth_width), +air_gap)
femm.mi_drawline(+width_mover/2, +air_gap, (+width_mover/2-tooth_width), +air_gap)

femm.mi_drawline(-width_mover/2, +air_gap, -width_mover/2, (+air_gap+tooth_heigth+tooth_width))
femm.mi_drawline(+width_mover/2, +air_gap, +width_mover/2, (+air_gap+tooth_heigth+tooth_width))

femm.mi_drawline((-width_mover/2+tooth_width), +air_gap, (-width_mover/2+tooth_width), (+air_gap+tooth_heigth))
femm.mi_drawline((+width_mover/2-tooth_width), +air_gap, (+width_mover/2-tooth_width), (+air_gap+tooth_heigth))


femm.mi_drawline((-width_mover/2+tooth_width), (+air_gap+tooth_heigth), (+width_mover/2-tooth_width), (+air_gap+tooth_heigth))
femm.mi_drawline(-width_mover/2, (+air_gap+tooth_heigth+tooth_width), +width_mover/2, (+air_gap+tooth_heigth+tooth_width))

# ----- Create mover region
femm.mi_addblocklabel(0, +air_gap+tooth_heigth+tooth_width/2)
femm.mi_selectlabel(0, +air_gap+tooth_heigth+tooth_width/2)
femm.mi_setblockprop(MaterialStator, 0, SizeMesh_Stator, 0, 0, 0, 0)
femm.mi_setgroup(0)
femm.mi_clearselected


# ********************* Stator Coil
# ------ Plus side
femm.mi_drawline((-width_mover/2+tooth_width), +air_gap+tooth_heigth-coil_heigth, (+width_mover/2-tooth_width), +air_gap+tooth_heigth-coil_heigth)

femm.mi_addcircprop('Coil', current[0], 1)

# ----- Create coil plus region
femm.mi_addblocklabel(0,+air_gap+tooth_heigth-coil_heigth/2)
femm.mi_selectlabel(0,+air_gap+tooth_heigth-coil_heigth/2)
femm.mi_setblockprop('Stator_Coil', 0, SizeMesh_Coil, 'Coil', 0, 3, turns)
femm.mi_setgroup(2)
femm.mi_clearselected

# ------ Minus side
femm.mi_drawline((-width_mover/2+tooth_width), +air_gap+tooth_heigth+tooth_width+coil_heigth, (+width_mover/2-tooth_width), +air_gap+tooth_heigth+tooth_width+coil_heigth)
femm.mi_drawline((-width_mover/2+tooth_width), +air_gap+tooth_heigth+tooth_width, (-width_mover/2+tooth_width), +air_gap+tooth_heigth+tooth_width+coil_heigth)
femm.mi_drawline((+width_mover/2-tooth_width), +air_gap+tooth_heigth+tooth_width, (+width_mover/2-tooth_width), +air_gap+tooth_heigth+tooth_width+coil_heigth)

# ----- Create coil minus region
femm.mi_addblocklabel(0,+air_gap+tooth_heigth+tooth_width+coil_heigth/2)
femm.mi_selectlabel(0,+air_gap+tooth_heigth+tooth_width+coil_heigth/2)
femm.mi_setblockprop('Stator_Coil', 0, SizeMesh_Coil, 'Coil', 0, 3, -turns)
femm.mi_setgroup(2)
femm.mi_clearselected

# ********************* Air-Gap
femm.mi_drawline(-width_mover/2, 0, -width_mover/2, +air_gap)
femm.mi_drawline(+width_mover/2, 0, +width_mover/2, +air_gap)

# ----- Create mover region
femm.mi_addblocklabel(0,+air_gap)
femm.mi_selectlabel(0,+air_gap)
femm.mi_setblockprop('Air', 0, SizeMesh_Air, 0, 0, 0, 0)
femm.mi_clearselected

# ********************* Outer Air
femm.mi_drawline(-2*width_mover, 2*width_mover, -2*width_mover, -2*width_mover)
femm.mi_drawline(2*width_mover, 2*width_mover, 2*width_mover, -2*width_mover)
femm.mi_drawline(2*width_mover, 2*width_mover, -2*width_mover, 2*width_mover)
femm.mi_drawline(2*width_mover, -2*width_mover, -2*width_mover, -2*width_mover)

# ----- Create mover region
femm.mi_addblocklabel(-2*width_mover+1,2*width_mover-1)
femm.mi_selectlabel(-2*width_mover+1,2*width_mover-1)
femm.mi_setblockprop('Air', 0, SizeMesh_Air_out, 0, 0, 0, 0)
femm.mi_clearselected

femm.mi_saveas('Electromagnet.fem')
femm.mi_zoomnatural()

# ================= SOLUTION OF THE PROBLEM
# Field = zeros(14, Npoint_Bsampling,size(current,2));
FluxLinkage = np.zeros([len(current), len(chi)])   
Energy = np.zeros([len(current), len(chi)])   
Coenergy = np.zeros([len(current), len(chi)])   
Inductance = np.zeros([len(current), len(chi)])   

for idx_c, val_c in enumerate(current):

    femm.mi_modifycircprop('Coil', 1, current[idx_c])
#    femm.mi_modifycircprop('Coil', val_c, 1)

    for idx_x, val_x in enumerate(chi):

        femm.mi_selectgroup(1)
        femm.mi_movetranslate(0, val_x)
        femm.mi_clearselected()


        femm.mi_analyze()
    
        # ================= POST-PROCESSING
        femm.mi_loadsolution()

#     x_Bsampling = linspace(-width_mover/2, width_mover/2, Npoint_Bsampling);

#     % ----- Post processing / Field sampling 
#     y_pos = air_gap/2;

#     for f = 1 : size(x_Bsampling,2)
#         x_pos = x_Bsampling(1,f);

#         Field(1:14, f, c)=mo_getpointvalues(x_pos,y_pos);
#         %	1: A_z           2: B_x      3: B_y
#         %   4: Sigma (Conductivity)      5: En (Energy)
#         %   5: H_x           6: H_y
#         %   7: Jeddy (Eddy current losses)      
#         %   9: Jsource (Source current density)
#         %   10: Mu_x         11: Mu_y
#         %   12: P_ohm (Ohmic losses)     13: P_hyst (Hysteresis losses)

#     end
    
        # ----- Post processing / Circuit solution 
        CircuitSolution = femm.mo_getcircuitproperties('Coil') 
        
        # Current[idx_c, idx_x] = CircuitSolution[0]
        # Voltage[idx_c, idx_x] = CircuitSolution[1]
 
        FluxLinkage[idx_c, idx_x] = CircuitSolution[2]
        Inductance[idx_c, idx_x] =  FluxLinkage[idx_c, idx_x] / val_c
        
        # ----- Post processing / Energy    
        femm.mo_selectblock(0, air_gap/2)
        Coenergy[idx_c, idx_x] = femm.mo_blockintegral(17)   # o 2??   # Pag 93 of the Manuale for the other values
        Energy[idx_c, idx_x] = FluxLinkage[idx_c, idx_x] * val_c - Coenergy[idx_c, idx_x]
        
        femm.mi_selectgroup(1)
        femm.mi_movetranslate(0, -val_x)
        femm.mi_clearselected()

 
if len(chi) != 1:
    plt.figure()
    for idx_c, val_c in enumerate(current):
        plt.plot(chi, FluxLinkage[idx_c, :]) 
    plt.xlabel(r'$\chi\,[mm]$')
    plt.ylabel(r'$\varphi(i,\chi)\,[Wb]$')
            
if len(current) != 1:
    plt.figure()
    for idx_x, val_x in enumerate(chi):
        plt.plot(current, FluxLinkage[:, idx_x]) 
    plt.xlabel(r'$i\,[A]$')
    plt.ylabel(r'$\varphi(i,\chi)\,[Wb]$')           

# plt.figure()
# for idx_c, val_c in enumerate(current):
#     plt.plot(chi, Coenergy[idx_c, :]) 
#     plt.plot(chi, Energy[idx_c, :], linestyle = '--') 
# plt.xlabel(r'$\chi\,[mm]$')
# plt.ylabel(r'$\varphi(i,\chi)\,[Wb]$')
            

# plt.figure()
# for idx_x, val_x in enumerate(chi):
#     plt.plot(current, Coenergy[:, idx_x]) 
# plt.xlabel(r'$i\,[A]$')
# plt.ylabel(r'$\varphi(i,\chi)\,[Wb]$')           



