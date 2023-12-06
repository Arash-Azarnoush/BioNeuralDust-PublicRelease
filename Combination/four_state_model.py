#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementation of the Williams et Al Four-State Model

Created on Sun Nov  6 17:33:32 2022

@author: majawitter
"""

import math
import matplotlib.pyplot as plt
#from izhikevich import Izhikevich

class FourStateModel:
    def __init__(self, gamma, Gd2, ep1, ep2, sigma_ret, w_loss, T_ChR2):
        self.gamma = gamma
        self.Gd2 = Gd2
        self.ep1 = ep1
        self. ep2 = ep2
        self.sigma_ret = sigma_ret
        self.w_loss = w_loss
        self.hc = 1.986446*(10**(-25))
        self.T_ChR2 = T_ChR2
        self.temp = 22
        self.temp_scale = True
        self.c1 = 1.0
        self.c2 = 0.0
        self.o1 = 0.0
        self.o2 = 0.0
        self.p = 0
        self.Is = []

    def getC1(self):
        return self.c1
    
    def getC2(self):
        return self.c2
    
    def getO1(self):
        return self.o1
    
    def getO2(self):
        return self.o2
    
    # rate constant for O1->C1 transition: 0.084, 0.108, 0.11
    def Gd1(self, v): #ms^(-1)
        gd1 = 0.075 + 0.043 * math.tanh((v+20)/(-20))
        if self.temp_scale:
            return gd1*(1.97**((self.temp-22)/10))
        else:
            return gd1

    # rate constant for C2->C1 transition: 0.004, 0.0004
    def Gr(self, v): #ms^(-1)
        gr = (4.34587 * (10**-5) * math.exp(-0.0211539274*v))
        if self.temp_scale:
            return gr*(2.56**((self.temp-22)/10))
        else:
            return gr
    
    def logphi(self, irradiance):
        if (irradiance>0):
            return math.log(1+(irradiance/0.024))
        else:
            return 0
    
    # rate constant for O1->O2 transition: e12d = 0.011, 0.0297, 0.03, 0.0
    def e12(self, logphi0): #ms^(-1)
        e12d = 0.011
        if self.temp_scale:
            return e12d*(1.1**((self.temp-22)/10)) + 0.005*logphi0
        else:
            return e12d + 0.005*logphi0

    # rate constant for O2->O1 transition: e21d = 0.008, 0.018, 0, 0.015
    def e21(self, logphi0): #ms^(-1)
        e21d = 0.008
        if self.temp_scale:
            return e21d*(1.95**((self.temp-22)/10)) + 0.004*logphi0
        else:
            return e21d + 0.004*logphi0
    
    # photon flux, number of photons per molecule per second
    def F(self, irr, wavelength): #ms^(-1)
        Ephoton = (1*(10**9)) * self.hc / wavelength #J
        flux = 1000*irr/Ephoton #W/m^2
        return (flux*self.sigma_ret/(self.w_loss*1000))

    # state-variable, time- and irradiance-dependent activation function for ChR2
    def So(self, irr):
        theta = 100*irr
        return 0.5 * (1 + math.tanh(120 * (theta - 0.1)))
    
    def reset(self):
        self.c1 = 1.0
        self.o1 = 0.0
        self.o2 = 0.0
        self.c2 = 0.0
        self.p = 0
        self.Is = []
    
    def getPhotocurrent(self, irr, V, dt, wavelength, g_ChR2, E_ChR2):   
        Gd1 = self.Gd1(V)
        logphi0 = self.logphi(irr)
        e12 = self.e12(logphi0)
        e21 = self.e21(logphi0)
        F = self.F(irr, wavelength)
        Gr = self.Gr(V)
        Gd2 = self.Gd2
        ep1 = self.ep1
        ep2 = self.ep2
        
        if self.temp_scale:
            Gd2 = Gd2*(1.77**((self.temp-22)/10))
            ep1 = ep1*(1.46**((self.temp-22)/10))
            ep2 = ep2*(2.77**((self.temp-22)/10))
        
        dP = (self.So(irr) - self.p)/self.T_ChR2
        dC1O1 = ep1*F*self.p*self.c1
        dO1C1 = Gd1*self.o1
        dO1O2 = e12*self.o1
        dO2O1 = e21*self.o2
        dO2C2 = Gd2*self.o2
        dC2O2 = ep2*F*self.p*self.c2
        dC2C1 = Gr*self.c2
        dC1  =  dC2C1 + dO1C1 - dC1O1
        dO1  =  dC1O1 + dO2O1 - dO1C1 - dO1O2
        dO2  =  dC2O2 + dO1O2 - dO2C2 - dO2O1
        dC2  =  dO2C2 - dC2O2 - dC2C1
        
        # apply changes multiplied by the time step
        self.c1 += dC1*dt
        self.o1 += dO1*dt
        self.o2 += dO2*dt
        self.c2 += dC2*dt
        self.p += dP*dt
                
        A = +10.6408
        B = -14.6408
        C = +42.7671
        I_ChR2 = g_ChR2 * ((A+B*math.exp(-V/C))) * (self.o1 + self.gamma * self.o2)
        
        self.Is.append(I_ChR2)
        
        return I_ChR2
    
    def getPlot(self):
        fig = plt.figure(figsize=(12,9), dpi=180)
        plt.plot(self.Is)
        plt.show()
