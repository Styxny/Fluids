# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 08:55:12 2018

@author: jstickney7

written in Python 3

Adv. Aerodynamics HW3 AE - 6015
Vortex lattice method fro swept wing
"""
import numpy as np
import matplotlib.pyplot as plt

"""
Define Variables
"""
AR = 5                          # Aspect Ratio
chord = 1                       # Sweep Angle in degrees
n = 3                           # Number of rows of panels used
m = 8                           # Number of columns of panels used
panel = n*m                     # number of panels
span = AR                       # AR is to span for chord =1
hspan = span/2
sweep = 45                      # Sweep angle in degrees
sw_r = np.deg2rad(sweep)        # Sweep angle in radians
UINF = 1                        # x axis flow velocity
p_wid = span/m                  # width of each panel
point_start = p_wid/2           # ycor start point for the calculation
p_len = chord/n                 # length of each panel
R_c = 0.001                     # Core radius
rho = 0.0023769                 # Density of air in slug/ft^3
Mach = 0.5                      # Mach Value
beta = np.sqrt(1-Mach**2)       # Compressibility Correction Factor

yc = np.arange(0,hspan+p_wid,p_wid)
yc = np.append(yc,[-x for x in yc if x!= 0])
yc = np.sort(yc)
ycc = yc

    
    
point = np.arange(point_start,hspan,p_wid)
point = np.append(point,-point)
point = np.sort(point)
pt = point
#print(point)
        

def Vortex_lat(a):
    AoA = np.deg2rad(a)             # Angle of Attack in degrees
    # Create Panel Class-----------------------------------------------------------
    class short_V:
        xcor = [0,0]
        ycor = [0,0]
    
    class long_V:
        xcor1 = [0,1000]
        xcor2 = [0,1000]
        ycor = [0,0]
    
    class Panel:
        cor = [0,0]
        
    
    # Generate Vortex objects------------------------------------------------------
    sV = []
    for obj in range(panel):
        obj = short_V()
        sV.append(obj)
    
    lV = []
    for obj in range(panel):
        obj = long_V()
        lV.append(obj)
    
    panels = []
    for obj in range(panel):
        obj = Panel()
        panels.append(obj)
    
#    print(len(lV))
    # Assign Ycor to each panel----------------------------------------------------
    #------------------------------------------------------------------------------
    # Vorticies--------------------------------------------------------------------
    for j in range(n):
        for i in range(m):
            sV[i+j*m].ycor = [yc[i],yc[i+1]] 
            lV[i+j*m].ycor = [yc[i],yc[i+1]]
            
            
    # Assign Xcor to vorticies-----------------------------------------------------
    for j in range(n):
        for i in range(m):
            xa = (1/beta)*(abs(sV[i+j*m].ycor[0]) * np.tan(sw_r) + p_len * (j + 0.25))
            xb = (1/beta)*(abs(sV[i+j*m].ycor[1]) * np.tan(sw_r) + p_len * (j + 0.25))
            sV[i+j*m].xcor = [xa,xb]
            lV[i+j*m].xcor1[0] = xa
            lV[i+j*m].xcor2[0] = xb
            
    # Generate Calculation points for each panel-----------------------------------
        
    for j in range(n):
        for i in range(m):
            panX = (1/beta)*    (abs(point[i]) * np.tan(sw_r) + p_len * (j + 0.75)) 
            panels[i+j*m].cor = [panX,point[i]]

#    plt.figure()
#    for i in range(panel):
#        plt.plot(panels[i].cor[0],panels[i].cor[1],'ro')
#    plt.show()
    
    # Create the A matrix----------------------------------------------------------
    A = np.zeros((panel,panel))
    
    # Create the values for the a matrix-------------------------------------------
    for j in range(panel):
        for i in range(panel):
           Px = panels[j].cor[0]
           Py = panels[j].cor[1]
           vxa = sV[i].xcor[0]                  # Short Vortex start xcor
           vya = sV[i].ycor[0]                  # Short Vortex start ycor
           vxb = sV[i].xcor[1]                  # Short Vortex end xcor
           vyb = sV[i].ycor[1]                  # Short Vortex end ycor
           lxb = lV[i].xcor1[1]                 # 1st Long Vortex end xcor
           lyb = lV[i].ycor[0]                  # 1st Long Vortex ycor
           lxb_2 = lV[i].xcor2[1]               # 2nd Long Vortex end xcor
           lyb_2 = lV[i].ycor[1]                # 2nd Long Vortex ycor
           r_1 = np.array([Px-vxa,Py-vya])
           r_2 = np.array([Px-vxb,Py-vyb])   
           r_3 = np.array([Px-lxb,Py-lyb])
           r_4 = np.array([Px-lxb_2,Py-lyb_2])
           mr_1 = np.linalg.norm(r_1)               # Magnitude of r_1
           mr_2 = np.linalg.norm(r_2)               # Magnitude of r_2
           mr_3 = np.linalg.norm(r_3)               # Magnitude of r_3
           mr_4 = np.linalg.norm(r_4)               # Magnitude of r_4
           dot12 = np.dot(r_1,r_2)
           dot13 = np.dot(r_1,r_3)
           dot24 = np.dot(r_2,r_4)
           # Short Finite Vortex Contribution--------------------------------------
           pong1 = ((mr_1+mr_2)*(1-(dot12/(mr_1*mr_2))))/ \
               ((mr_1*mr_2)**2 - dot12**2 + (R_c**2)*(mr_1**2 + mr_2**2 - 2*dot12))
           A_1 = (1/(4*np.pi))*np.cross(r_1,r_2)* pong1
           
           # Left Side semi-infinite Vortex Contribution---------------------------
           pong2 = ((mr_1+mr_3)*(1-(dot13/(mr_1*mr_3))))/ \
               ((mr_1*mr_3)**2 - dot13**2 + (R_c**2)*(mr_1**2 + mr_3**2 - 2*dot13))
           A_2 = (1/(4*np.pi))*np.cross(r_3,r_1)* pong2
           
           # Right side semi-infintie Vortex Contribution--------------------------
           
           pong3 = ((mr_4+mr_2)*(1-(dot24/(mr_4*mr_2))))/ \
               ((mr_4*mr_2)**2 - dot24**2 + (R_c**2)*(mr_4**2 + mr_2**2 - 2*dot24))
           A_3 = (1/(4*np.pi))*np.cross(r_2,r_4)* pong3
           
           
           # Giving index of A a value---------------------------------------------
           A[j,i] = A_1 + A_2 + A_3
#    print(np.shape(A))      
    b = [-UINF*np.sin(AoA) for x in range(panel)]
    Gamma = np.dot(np.linalg.inv(A),b)
#    drag = [rho*abs(b[x])*Gamma[x]*p_wid for x in range(panel)]
#    Drag = sum(drag)
    L_us = rho * UINF *Gamma
    S = chord*hspan*2
    L_F = L_us * p_wid
    L = sum(L_F)
#    CD = Drag/(0.5*rho*UINF**2*S)
#    print(A)
    C_L = L/(0.5 * rho * UINF**2 *S)
    return C_L, A


a = [0,2,4,6,8]
CL_cal, A =Vortex_lat(a)
print(CL_cal)
plt.figure()
plt.plot(a,CL_cal,'r', label='Vortex Lattice CL')
plt.xlim(0,7)
plt.xlabel('Angle of Attack in Degrees')
plt.ylabel('Lift coefficient')
plt.legend(loc='upper left')
plt.ylim(0,0.5)
plt.show()
#    