from numpy import arange, zeros, array
import matplotlib.pyplot as plt
from photon_emission import PhotonEmission
from ser import SmoothEndoRet
from calcium import CalciumCluster
from M13 import M13Phage

from neuron import Neuron
from izhikevich import Izhikevich
from channelrhodopsin import Channelrhodopsin
from four_state_model import FourStateModel
from vpython import *
import copy

# =============================================================================
# Nano Bubble Component
# =============================================================================

# =============================================================================
# Initial simulation parameters for the nano bubble (change as needed)
# =============================================================================
Pout = 0.2E-3 #moles
Vrest = -0.050 #50mV
Vmax = -0.010

S = 1E-5 #10 microseconds   #NOTE: 1E-2 would be 10,000 micorsecond which is 10 miliseconds, running this takes about 2 to 3 minutes on my laptop so it is acceptable in my opinion (we still need to up it to atleast 100 miliseconds)
dt = 1E-9 #1 ns

# =============================================================================
#  Set variables and objects needed for simulation
# =============================================================================

#instantiate calcium model class
calcium_model = CalciumCluster(Vrest,Vmax,Pout)

#set initial ratio and extra/intracellular concentrations
(init_ratio, init_Pin, init_Pout) = calcium_model.updateConcentrations(Vrest)

#get the max extracellular concentration of ion (i.e, when Vt = Vmax)
Pout_max = calcium_model.getPoutMax()

#setup init photon emission variables
diffusion_coefficient = 1E-8
aequorin_radius = 1E-8
photon_emiter = PhotonEmission(init_Pout, Pout_max, diffusion_coefficient, aequorin_radius, dt)

#file path of the excel provided by Arash
filePath = 'M13_Voltage_V3.xlsx'

#instantiate bacteriophage model object
bacteriophage = M13Phage(filePath)

# =============================================================================
#  zap() method will run a single simulation for a given soundwave intensity
#  and a given number of ser channels.
# =============================================================================

#time array with S/dt elements (i.e. if S = 1μs and dt = 1ns, len(S_array) = 1000)
S_array = arange(0, S, dt)

#this function performs the entire simulation from a given input ultrasound intensity returns
#value arrays for the membrane voltage, extracellular concentration and light intensity over time
def zap(input_soundwave, channels):
    
    #set the input intensity on the M13 (determines the output source voltage)
    bacteriophage.setUltrasound(input_soundwave, channels)

    #get the new source voltage (i.e., what is being applied to the ser)
    Vb = bacteriophage.getVoltage()

    #instantiate SER class
    ser = SmoothEndoRet(channels, Vrest, Vmax, dt)
    
    #start with a voltage Vb
    ser.updateVoltage(Vb)
    
    #initial light intensity
    init_output_intensity = photon_emiter.getLightIntensity(init_Pout,dt)

    #initialize arrays to store values for each time increments 
    Vt_array = zeros(len(S_array))
    Pt_out_array = zeros(len(S_array))
    light_array = zeros(len(S_array))  
    
    #add initial values to first index of arrays
    Vt_array[0] = Vrest * 1000
    Pt_out_array[0] = init_Pout * 1000
    light_array[0] =  init_output_intensity/1000
    
    for t in range(1,len(S_array)):
        
        #update ser voltage w/ Leaky Integrate and Fire Model variation
        Vt = ser.nextTimeStep()

        #add new membrane voltage to array at corresponding time step (in mV)
        Vt_array[t] = Vt * 1000
        
        #update/get real concentrations according to new voltage
        ratio, Pin_t, Pout_t = calcium_model.updateConcentrations(Vt)

        #update light intensity according to new extracellular concentration
        output_intensity = photon_emiter.getLightIntensity(Pout_t,dt)

        #add new concentration to array (in mM)
        Pt_out_array[t] = Pout_t  * 1000

        #add new light intensity to array (in mW/mm^2)
        light_array[t] = output_intensity/1000
    
    return (Vt_array, 
            Pt_out_array, 
            light_array)#<----- THIS WOULD BE THE INPUT FOR THE OPSIN NEURON

#this function gets the highest blue light intensity values in the light_array
#it then takes the average of these values (which will be used for the neuron component)
def avg_light_intensity(blue_light_array):
    #sort the array, in ascending order
    temp = array(blue_light_array)
    temp[::-1].sort()

    #gets the number of indexes to use in the average calculation(top 5% for for this case)
    cut_off = 5 #(percentage) - change as needed
    treshold = (cut_off*len(temp))//100
    print(len(blue_light_array))
    print(treshold)

    #iterates through the top 10% values and calculates the average
    avg = 0
    for x in range(treshold):
        avg += temp[x]
    
    return avg/treshold

# =============================================================================
# Neuron Component
# =============================================================================

# =============================================================================
# Simulation Variables
# =============================================================================
t = 0.0 # start time (ms)
dt_neuron = 0.02 #time step (ms) [0.15, 0.015, 0.001]
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
membrane_graph = graph(xtitle="Time (μs)", ytitle="Membrane potential (mV)", align='left')
membrane_plot = gcurve(color=color.green)

extracellular_graph = graph(xtitle="Time (μs)", ytitle="Extracellular Concentration (mM)", align='left')
extracellular_plot = gcurve(color=color.green)

blue_graph = graph(xtitle="Time (μs)", ytitle="Blue light intensity (mW/mm^2)", align='left')
blue_plot = gcurve(color=color.green)

ChR2_graph = graph(xtitle="Time (ms)", ytitle="ChR2 Current (pA/pF)", width = 550, height = 250, xmin = 0, xmax = sim_time, align='left')
ChR2_plot = gcurve(color=color.red, graph = ChR2_graph)

neuron_graph = graph(xtitle="Time (ms)", ytitle="Membrane Potential (mV)", width = 550, height = 250, xmin = 0, xmax = sim_time, align='left')
neuron_plot = gcurve(color=color.blue, graph = neuron_graph)

light_graph = graph(xtitle="Time (ms)", ytitle="Light Intensity (mW/mm^2)", width = 550, height = 250, xmin = 0, xmax = sim_time, align='left')
light_plot = gcurve(color=color.yellow, graph = light_graph)

# =============================================================================
# Set Scene
# =============================================================================
scene.title = "Combined Framework\n\n"
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
izhikevich = Izhikevich(dt_neuron, iz_a, iz_b, iz_c, iz_d, ap_threshold)

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
ChR2 = Channelrhodopsin(four_state_model, dt_neuron, E_ChR2, g_ChR2, wavelength, holding_potential, True)

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
    global pulse_times, pulse_widths, light_delay, pulse_width, average_Intensity
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


###############################################
####           DRAW NANO-BUBBLE            ####
###############################################
bubble_radius=100
#radius = 1000 (or slider value)
#opacity = 65% transparent
nano_bubble=sphere(color=vector(0.5,0.7,0.9), #blue
                   pos=vector(-200,0,0),      
                   radius=bubble_radius,
                   opacity=0.35, 
                   shininess=5)

#blue_light = attach_light(nano_bubble, local_pos = vector(0,0,0), color=color.red, intensity = light_intensity)
# =============================================================================
# RUN SIMULATION TESTS - this section runs the simulations
# =============================================================================

#varaibles for the number of channels and the input soundwave
number_channels = 1000
soundwave_intensity = 2.0

#run zap() -> get output 
(Vt_arr, extra_conc_arr, light_arr) = zap(soundwave_intensity, number_channels)

#get M13's volatge (for graph labels)
v = bacteriophage.getVoltage()*1000

light_intensity = avg_light_intensity(light_arr)
print(light_intensity)

#displays the caculated average light intensity generated by the bubble
scene.append_to_caption('\nAverage Maximum Blue Light intensity: '+str(round(light_intensity,3)) + ' mW/mm^2\n') 

# Add points to the curve of the membrane potential, extracellular concentration, 
# and blue light intensity graphs using the caculates arrays by the zap function
for i in range(len(S_array)):
    membrane_plot.plot(pos=(S_array[i], Vt_arr[i]))
    extracellular_plot.plot(pos=(S_array[i], extra_conc_arr[i]))
    blue_plot.plot(pos=(S_array[i], light_arr[i]))

#initializes the last 3 graphs for the neuron
ChR2_plot.plot([0,0])
neuron_plot.plot([0,0])
light_plot.plot([0,0])

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
            
            t += dt_neuron