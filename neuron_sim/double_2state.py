#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Attempted Implementation of Double Two-State Opsin Model. 

Abandoned in favor of Four-State Model.

Unsuccessful & incomplete, but left here if needed for future reference

Created on Sun Nov 13 19:46:16 2022

@author: majawitter
"""

import math
import matplotlib.pyplot as plt
from neuron import Neuron
from izhikevich import Izhikevich
from vpython import *

# =============================================================================
# Start using RSRS final params (double reciprocal addition – Equation 16), 
# note other option of PP params (double product – Equation 15)
# both RSRS and PP also have intermediate params listed
# =============================================================================

# [p1, p2, p3]
TO_I_params = [1.81, 1.17, 0.021] #PP: [1.93, 0.88, 0.030]

# [p1, p2, p3]
TO_V_params = [23.14, -0.39, 13.19] #PP: [0.63, -88.67, 8.37]

# [p1, p2, p3, p4, p5, p6]
TR_I_params = [10, 0.56, -1.58, 0.87, 1.96, 0.11] #PP: [6.73, 0.5, 1.98, 0.11, -1.28, 0.88]

# [p1, p2, p3]
TR_V_params = [99.74, -38.69, 12.02] #PP: [1.66, -64.54, 28.55]

# [p1, p2]
O_inf_params = [3.38, 0.62] #PP: [3.44, 0.68]

# [p1, p2, p3]
R_inf_params = [1.96, 0.12, 0.77] #PP: [2.25, 0.065, 0.75]

# [p1, p2, p3]
G_params = [10.77, 1.25, 44.52] #PP: [9.10, 1.27, 41.47]
# p1 = 1 is another param solution

g_ChR2 = 1 # 10.77 (RSRS) or 9.10 (PP) is another param solution

E_ChR2 = 0

def TO_I(I):
    p1, p2, p3 = TO_I_params[0], TO_I_params[1], TO_I_params[2]
    return (p3/(1 + (math.exp(p1/p2)*math.pow(I, (1/p2*math.log(10))))))

def TR_I(I):
    p1, p2, p3, p4, p5, p6 = TR_I_params[0], TR_I_params[1], TR_I_params[2], TR_I_params[3], TR_I_params[4], TR_I_params[5]
    #y = math.pow(I, (-1/p4*math.log(10)))
    c1 = (p2/(1 + (math.exp(p3/p4)*math.pow(I, (-1/p4*math.log(10))))))
    c2 = ((1-p2)/(1 + (math.exp(p5/p6)*math.pow(I, (-1/p6*math.log(10))))))
    return (p1*(1.0 - c1 - c2))

def TO_V(V):
    p1, p2, p3 = TO_V_params[0], TO_V_params[1], TO_V_params[2]
    return (p1/(1 + math.exp(-(V - p2)/p3)))

def TR_V(V):
    p1, p2, p3 = TR_V_params[0], TR_V_params[1], TR_V_params[2]
    return (p1/(1 + math.exp(-(V - p2)/p3)))

def TX_IV(TX_I, TX_V):
    #RSRS
    return math.pow((math.pow(TX_I, -1) + math.pow(TX_V, -1)), -1)
    #PP
    #return (TX_I * TX_V)

def O_inf(I):
    p1, p2 = O_inf_params[0], O_inf_params[1]
    return (1/(1 + (math.exp(p1/p2)*math.pow(I, (-1/p2*math.log(10))))))

def R_inf(I):
    p1, p2, p3 = R_inf_params[0], R_inf_params[1], R_inf_params[2]
    return (1.0 - (p3/(1 + (math.exp(p1/p2)*math.pow(I, (-1/p2*math.log(10)))))))

def G(V):
    p1, p2, p3 = G_params[0], G_params[1], G_params[2]
    Vd = V - E_ChR2
    n = p1*(1-(p2*(math.exp(-(Vd)/p3))))
    return (n/Vd)

def heaviside(x):
    if x > 0:
        return 1
    else:
        return 0

def O_on(t, t_on, t_off, Oi, O0, I, V):
    if (heaviside(t-t_on) == 0) or (heaviside(t_off-t) == 0):
        return 0
    else:
        TOI = TO_I(I)
        TOV = TO_V(V)
        TOIV = TX_IV(TOI, TOV)
        return (heaviside(t-t_on)*heaviside(t_off-t)*(Oi-((Oi-O0)*math.exp(-(t-t_on)/(TOIV)))))

def O_off(t, t_off, O_on, V):
    if (heaviside(t-t_off) == 0):
        return 0
    else:
        TO0 = TO_I(0)
        TOV = TO_V(V)
        TO0V = TX_IV(TO0, TOV)
        return (O_on*math.exp(-(t-t_off)/t_off)*heaviside(t-t_off))

def R_on(t, t_on, t_off, Ri, R0, I, V):
    if (heaviside(t-t_on) == 0) or (heaviside(t_off-t) == 0):
        return 0
    else:
        TRI = TR_I(I)
        TRV = TR_V(V)
        TRIV = TX_IV(TRI, TRV)
        x = (Ri-((Ri-R0)*math.exp(-(t-t_on)/(TRIV))))
        return (x*heaviside(t-t_on)*heaviside(t_off-t))

def R_off(t, t_off, R_on, V):
    if (heaviside(t-t_off) == 0):
        return 0
    else:
        TR0 = TR_I(0)
        TRV = TR_V(V)
        TR0V = TX_IV(TR0, TRV)
        x = (1-((1-R_on)*math.exp(-(t-t_off)/TR0V)))
        return x*heaviside(t-t_off)

def i_ChR2(V, O_on, O_off, R_on, R_off):
    return (g_ChR2*G(V)*(O_on+O_off)*(R_on+R_off)*(V-E_ChR2))

def irradiance(I, t, t_on, t_off):
    if (t < t_on) or (t > t_off):
        return 0
    else:
        return I

def findT_recov():
    Trecov = 0.001 #1ms
    t = 0.0
    dt = 0.00015
    V = -70
    I = 5500
    t_on1 = 0.25
    t_off1 = (t_on1 + 0.5)
    t_on2 = (t_off1 + Trecov)
    t_off2 = (t_on2 + 0.5)
    tF = 2.5
    O0 = 0
    R0 = 1
    
    i_s = []
    I_p1 = 0
    I_p2 = 0
    v_s = []

    izhikevich = Izhikevich(dt, 0.02, 0.2, -65, 8, 30.0)
    
    ratio = 0
    
    #while not math.isclose(ratio, (1-exp(-1)), rel_tol=(10**(-6)), abs_tol=(10**(-6))):
    for x in range (0,1):
        #print(t_on1, t_off1, t_on2, t_off2)
        # 0.25 0.75 0.751 1.251
        while t < t_off2:
            if t < t_off1:
                t_on = t_on1
                t_off = t_off1
            else:
                t_on = t_on2
                t_off = t_off2
            
            Oi = O_inf(I)
            Ri = R_inf(I)
            
            Oon = O_on(t, t_on, t_off, Oi, O0, I, V)
            Oon_t_off = O_on(t_off, t_on, t_off, Oi, O0, 0, V)
            Ooff = O_off(t, t_off, Oon_t_off, V)
            
            Ron = R_on(t, t_on, t_off, Ri, R0, I, V)
            Ron_t_off = R_on(t_off, t_on, t_off, Ri, R0, 0, V)
            Roff = R_off(t, t_off, Ron_t_off, V)
            
            i = i_ChR2(V, Oon, Ooff, Ron, Roff)
            izhikevich.nextTimeStep(abs(i))
            V = izhikevich.getMembranePotential()
            
            i_s.append(i)
            if (t < t_off1) and (abs(i)>abs(I_p1)):
                I_p1 = i
            elif (t > t_on2) and (abs(i)>abs(I_p2)):
                I_p2 = i
            #print(V, i)
            
            t += dt
            
        print('I_p1: ' + str(I_p1))
        print('I_p2: ' + str(I_p2))
        print(I_p2/I_p1)
        print((1-exp(-1)))
        print(math.isclose((I_p2/I_p1), (1-exp(-1)), rel_tol=(10**(-6)), abs_tol=(10**(-6))))
        
        ratio = I_p2/I_p1
        
        Trecov += 0.001
        t_on2 = t_off1 + Trecov
        t_off2 = t_on2 + 0.5
        izhikevich = Izhikevich(dt, 0.02, 0.2, -65, 8, 30.0)
        V = -70
        t = 0.0
        
        fig = plt.figure(figsize=(12,9), dpi=180)
        plt.plot(i_s)
        plt.show()
        

findT_recov()    

# =============================================================================
# Params to set in simulation:
#     I: light intensity (time dependant?)
#     V: membrane potential (time dependant)
#     t_on: onset of optical pulse
#     t_off: offset of optical pulse
#     O0, R0: initial values of O and R at t = t_on (0, 1 when fully dark adapted)
# =============================================================================

t = 0.0
dt = 0.00015
V = -60
I = 5500
t_on = 0.25
t_off = 0.75
tF = 1.5
O0 = 0
R0 = 1

i_s = []
v_s = []

izhikevich = Izhikevich(dt, 0.02, 0.2, -65, 8, 30.0)

# =============================================================================
# while t < tF:
#     print(t)
#     
#     Oi = O_inf(I)
#     Ri = R_inf(I)
#     
#     Oon = O_on(t, t_on, t_off, Oi, O0, I, V)
#     Oon_t_off = O_on(t_off, t_on, t_off, Oi, O0, 0, V)
#     Ooff = O_off(t, t_off, Oon_t_off, V)
#     
#     Ron = R_on(t, t_on, t_off, Ri, R0, I, V)
#     Ron_t_off = R_on(t_off, t_on, t_off, Ri, R0, 0, V)
#     Roff = R_off(t, t_off, Ron_t_off, V)
#     
#     i = i_ChR2(V, Oon, Ooff, Ron, Roff)
#     izhikevich.nextTimeStep(abs(i))
#     V = izhikevich.getMembranePotential()
#     
#     i_s.append(i)
#     print(V, i)
#     
#     t += dt
#     
# fig = plt.figure(figsize=(12,9), dpi=180)
# plt.plot(i_s)
# plt.show()
# =============================================================================



