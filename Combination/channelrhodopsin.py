#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains constant/know ChR2 parameters (used across a variety of models)

Relies on ChR2 model class for current logic

Created on Thu Nov  3 20:23:44 2022

@author: majawitter
"""

import math
#from four_state_model import FourStateModel

class Channelrhodopsin:
    def __init__(self, model, dt, E_ChR2, g_ChR2, wavelength, holding_potential, Vclamp):
        self.model = model
        self.dt = dt
        self.I_ChR2 = 0.0 #mA/cm^2
        self.E_ChR2 = E_ChR2
        self.g_ChR2 = g_ChR2
        self.wavelength = wavelength
        self.holding_potential = holding_potential
        self.Vclamp = Vclamp
        
    def reset(self, model, wavelength, holding_potential):
        self.model = model
        self.t = 0.0
        self.I_ChR2 = 0.0 #mA/cm^2
        self.wavelength = wavelength
        self.holding_potential = holding_potential
        
    def getModel(self):
        return self.model
        
    def getHoldingPotential(self):
        return self.holding_potential
    
    def setHoldingPotential(self, holding_potential):
        self.holding_potential = holding_potential
    
    # =============================================================================
    # Calls on model for updated photocurrent, passes either the holding potential or
    # the neuron's membrane potential (depending on if Vclamp is set)
    # Note: purpose of voltage clamp & which voltage is intended to be passed to the 
    # ChR2 model is still unclear to me. Most papers have been using a voltage clamp though
    # and my results with it show better correspondance so that is what I have been using.
    # =============================================================================
    def nextTimeStep(self, v, irr):
        if self.Vclamp:
            V = self.holding_potential
        else:
            V = v
        self.I_ChR2 = self.model.getPhotocurrent(irr, V, self.dt, self.wavelength, self.g_ChR2, self.E_ChR2)
                    
    def getI(self):
        return self.I_ChR2