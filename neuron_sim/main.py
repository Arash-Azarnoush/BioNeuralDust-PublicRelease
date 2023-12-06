#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main program to run Neuron Simulation

Created on Thu Nov  3 20:26:07 2022

@author: majawitter
"""

from neuron import Neuron
from izhikevich import Izhikevich
from channelrhodopsin import Channelrhodopsin
from four_state_model import FourStateModel
from vpython import *
import copy

# =============================================================================
# Simulation Variables
# =============================================================================
t = 0.0 # start time (ms)
dt = 0.02 #time step (ms) [0.15, 0.015, 0.001]
sim_time = 600.0

# simulate consecutive light pulses by setting single_pulse to False, not incorportated into control panel
init_pulse_times = [0, 20, 40]
init_pulse_widths = [6, 6, 6]
pulse_times = copy.copy(init_pulse_times)
pulse_widths = copy.copy(init_pulse_widths)
single_pulse = True

light_delay = 100 # ms
pulse_width = 400 # ms
light_intensity = 5.5 # mW/mm^2
wavelength = 470 # nm
if not single_pulse:
    light_delay = pulse_times.pop(0)
    pulse_width = pulse_widths.pop(0)

# =============================================================================
# Initialize Graphs
# =============================================================================
ChR2_graph = graph(xtitle="Time (ms)", ytitle="ChR2 Current (pA/pF)", width = 550, height = 250, xmin = 0, xmax = sim_time, align='left')
ChR2_plot = gcurve(color=color.red, graph = ChR2_graph)

neuron_graph = graph(xtitle="Time (ms)", ytitle="Membrane Potential (mV)", width = 550, height = 250, xmin = 0, xmax = sim_time, align='left')
neuron_plot = gcurve(color=color.blue, graph = neuron_graph)

light_graph = graph(xtitle="Time (ms)", ytitle="Light Intensity (mW/mm^2)", width = 550, height = 250, xmin = 0, xmax = sim_time, align='left')
light_plot = gcurve(color=color.yellow, graph = light_graph)

ChR2_plot.plot([0,0])
neuron_plot.plot([0,0])
light_plot.plot([0,0])

# =============================================================================
# Set Scene
# =============================================================================
scene.title = "Neuron Simulation \n\n"
scene.background=color.black
scene.width=640
scene.height=440
running = False

# =============================================================================
# Initializing Izhikevich Model
# =============================================================================
iz_a = 0.02
iz_b = 0.2
iz_c = -65
iz_d = 8
ap_threshold = 30.0
izhikevich = Izhikevich(dt, iz_a, iz_b, iz_c, iz_d, ap_threshold)

# =============================================================================
# Initializing Four-State ChR2 Model
# =============================================================================
gamma = 0.1
Gd2 = 0.05
ep1 = 0.8535
ep2 = 0.14
sigma_ret = 12*(10**-20)
w_loss = 1.3
T_ChR2 = 1.3
four_state_model = FourStateModel(gamma, Gd2, ep1, ep2, sigma_ret, w_loss, T_ChR2)

# =============================================================================
# Create Channelrhodopsin Object
# =============================================================================
E_ChR2 = 0.0
g_ChR2 = 2.0/5
holding_potential = -70.0 # mV
ChR2 = Channelrhodopsin(four_state_model, dt, E_ChR2, g_ChR2, wavelength, holding_potential, True)

# =============================================================================
# Create Neuron Object
# =============================================================================
soma_radius = 18
neuron = Neuron(soma_radius, izhikevich, ChR2)

# =============================================================================
# Display for current time & membrane potential
# =============================================================================
V_display = label( pos=vec(10,40,0), 
      pixel_pos=True, 
      align='left',
      text='Membrane Potential: ' + str(neuron.getMembraneModel().getMembranePotential()) + ' mV \nTime: ' + str(t) + ' ms' )

# =============================================================================
# Return the current Light Intensity
# =============================================================================
def getIrradiance():
    global pulse_times, pulse_widths, light_delay, pulse_width
    if (not single_pulse) and (len(pulse_times) > 0):
        if (t > (light_delay + pulse_width)):
            light_delay = pulse_times.pop(0)
            pulse_width = pulse_widths.pop(0)
    if (t < light_delay) or (t > (light_delay + pulse_width)):
        return 0
    else:
        return light_intensity

# =============================================================================
# Light Source - for visual representation of input
# =============================================================================
light_source = local_light(pos=vector(0,0,(neuron.soma.radius*-1)+100),
                           color=color.blue,
                           visible=False)

# =============================================================================
# Buttons
# =============================================================================
# start or stop the simulation
def toggleSim(b):
    global running
    running = not running
    if running:
        b.text = 'Stop'
    else:
        b.text = 'Start'

button_start = button(text='Start', pos=scene.title_anchor, bind=toggleSim)

# reset the simulation
def resetSim(b):
    global t, izhikevich, ChR2, neuron, four_state_model, running, ChR2_graph, ChR2_plot, light_graph, light_plot, neuron_graph, neuron_plot, pulse_times, pulse_widths, light_delay, pulse_width
    if running:
        running = not running
        button_start.text = 'Start'
    t = 0.0
    
    if not single_pulse:
        pulse_times = copy.copy(init_pulse_times)
        pulse_widths = copy.copy(init_pulse_widths)
        light_delay = pulse_times.pop(0)
        pulse_width = pulse_widths.pop(0)
    
    izhikevich.reset()
    four_state_model.reset()
    ChR2.reset(four_state_model, wavelength, holding_potential)
    neuron.reset(izhikevich, ChR2)
    V_display.text = 'Membrane Potential: ' + str(neuron.getMembraneModel().getMembranePotential()) + ' mV \nTime: ' + str(t) + ' ms'
    light_source.visible = False

    ChR2_plot.delete()
    ChR2_plot = gcurve(color=color.red, graph = ChR2_graph)
    
    neuron_plot.delete()
    neuron_plot = gcurve(color=color.blue, graph = neuron_graph)
    
    light_plot.delete()
    light_plot = gcurve(color=color.yellow, graph = light_graph)
    return

button_reset = button(text='Reset', pos=scene.title_anchor, bind=resetSim)

# =============================================================================
# Control Panel
# =============================================================================
def setSimTime(x):
    global sim_time, ChR2_graph, neuron_graph
    sim_time = x.number
    ChR2_graph.xmax = x.number
    neuron_graph.xmax = x.number
    light_graph.xmax = x.number
    resetSim(x)
    
def setLightDelay(x):
    global light_delay
    light_delay = x.number
    resetSim(x)
    
def setPulseWidth(x):
    global pulse_width
    pulse_width = x.number
    resetSim(x)
    
def setLightIntensity(x):
    global light_intensity
    light_intensity = x.number
    resetSim(x)
    
def setHoldingPotential(x):
    global holding_potential
    if x.number > 0:
        holding_potential = x.number*(-1)
    else:
        holding_potential = x.number
    resetSim(x)

def setWavelength(x):
    global wavelength
    wavelength = x.number
    resetSim(x)

scene.caption = "Simulation Time: \n"
sim_time_input = winput(bind=setSimTime,
                        prompt="Enter here",
                        type="numeric",
                        text=sim_time,
                        width=100)
scene.append_to_caption(' ms')

scene.append_to_caption('\nLight Delay: \n') 
light_delay_input = winput(bind=setLightDelay,
                        prompt="Enter here",
                        type="numeric",
                        text=light_delay,
                        width=100)
scene.append_to_caption(' ms')

scene.append_to_caption('\nPulse Width: \n') 
pulse_width_input = winput(bind=setPulseWidth,
                        prompt="Enter here",
                        type="numeric",
                        text=pulse_width,
                        width=100)
scene.append_to_caption(' ms')

scene.append_to_caption('\nLight Intensity: \n') 
light_intensity_input = winput(bind=setLightIntensity,
                        prompt="Enter here",
                        type="numeric",
                        text=light_intensity,
                        width=100)
scene.append_to_caption(' mW/mm^2')

scene.append_to_caption('\nHolding Potential: \n') 
holding_potential_input = winput(bind=setHoldingPotential,
                        prompt="Enter here",
                        type="numeric",
                        text=holding_potential,
                        width=100)
scene.append_to_caption(' mV')

scene.append_to_caption('\nWavelength: \n') 
wavelength_input = winput(bind=setWavelength,
                        prompt="Enter here",
                        type="numeric",
                        text=wavelength,
                        width=100)
scene.append_to_caption(' nm\n')

# =============================================================================
# Loop running the simulation
# =============================================================================
while True:
    while t < sim_time:
        if running:
            neuron.nextTimeStep(t, getIrradiance())
                
            if getIrradiance() > 0:
                light_source.visible = True
            else:
                light_source.visible = False
                
            V_display.text = 'Membrane Potential: ' + str(neuron.getMembraneModel().getMembranePotential()) + ' mV \nTime: ' + str(t) + ' ms'
            
            # update graphs in real time
            ChR2_plot.plot([t, neuron.getOpsinModel().getI()])
            neuron_plot.plot([t, neuron.getMembraneModel().getMembranePotential()])
            light_plot.plot([t,getIrradiance()])
            
            t += dt