#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementation of the Izhikevich Neuron Model

Created on Thu Nov  3 14:10:44 2022

@author: majawitter
"""
import math
import matplotlib.pyplot as plt

class Izhikevich:
    
    def __init__(self, time_step, time_scale, sensitivity, potential_reset, recovery_reset, ap_threshold):
        self.dt = time_step
        
        # time scale of recovery variable u
        self.a = time_scale
        # sensitivity of u to the fluctuations in v
        self.b = sensitivity
        # after-spike reset value of v
        self.c = potential_reset
        # after-spike reset value of u
        self.d = recovery_reset
        
        # stable resting potential, -70mV
        self.v_rest_minus = (12.5 * self.b) - 62.5 - (12.5 * math.sqrt((self.b ** 2) - (10 * self.b) + 2.6))

        # action potential firing threshold
        self.v_thresh = ap_threshold
        
        # membrane potential, initially equal to the stable resting potential
        self.v = self.v_rest_minus

        # membrane recovery variable
        self.u = self.b*self.v
        
        self.I = 0.0

        self.vs = []
        self.spike_times = []
        
        self.spiked = False
        self.converged = False
        
    def reset(self):
        # membrane potential, initially equal to the stable resting potential
        self.v = self.v_rest_minus

        # membrane recovery variable
        self.u = self.b*self.v
        
        self.I = 0.0

        self.vs = []
        self.spike_times = []
        
        self.spiked = False
        self.converged = False
        
    def getRestingPotential(self):
        return self.v_rest_minus

    def updateMembranePotential(self, t, I):
        self.vs.append(self.v)
        
        dv = (0.04 * self.v ** 2) + (5.0 * self.v) + 140.0 - self.u + I
        du = self.a * ((self.b * self.v) - self.u)
            
        # apply changes multiplied by the time step
        self.v += dv * self.dt
        self.u += du * self.dt
        
        # detect spikes, and since we don't have other neurons to consider here we only mark down the time when a spike is detected
        if self.v >= self.v_thresh:
            #spike
            self.spike_times.append(t)
            self.spiked = True
            
            #reset
            self.v = self.c
            self.u = self.u + self.d
        
        # check if membrane potential has converged to resting potential
        percent_difference = (abs(abs(self.v) - abs(self.v_rest_minus)))/((abs(self.v) + abs(self.v_rest_minus))/2) * 100
        if percent_difference <= 0.5:
            self.converged = True
        else:
            self.converged = False
        
    def hasSpiked(self):
        return self.spiked
    
    def isSpiking(self):
        if self.spiked and (self.v <= self.v_rest_minus):
            self.spiked = False
        if (self.v >= -55): 
            return True
        else:
            return False
    
    def isConverged(self):
        return self.converged
    
    def setSpiked(self, val):
        self.spiked = val
        
    def getMembranePotential(self):
        return self.v
    
    def getRestingPotential(self):
        return self.v_rest_minus
    
    def getPlot(self):
        fig = plt.figure(figsize=(12,9), dpi=180)
        plt.plot(self.vs)
        plt.show()
    