from numpy import arange, zeros
import matplotlib.pyplot as plt
from photon_emission import PhotonEmission
from ser import SmoothEndoRet
from calcium import CalciumCluster
from M13 import M13Phage


# =============================================================================
# Initial simulation parameters (change as needed)
# =============================================================================

Pout = 0.2E-3 #moles
Vrest = -0.050 #50mV
Vmax = -0.010

S = 1E-6 #2 microseconds
dt = 1E-9 #ns



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





# =============================================================================
# RUN SIMULATION TESTS - this section runs the simulations for every possible
#                        number of channels, and tests each with 3 different 
#                        ultrasound intensity inputs
# =============================================================================

#initialize graph figures
fig, axs = plt.subplots(3,2,sharex=True,sharey='row',figsize=(20, 12))
plt.rc('figure', titlesize=25)
fig.suptitle(f'Applying different voltages for {S*1E6}μs with dt = {dt*1E9}ns and \n initial extracellular concentration of {Pout*1E3} millimoles')

#axs[0] = voltages graphs
#axs[1] = extracellular concentrations graphs
#axs[2] = light intensities graphs


#dictionary with the channels to test and their different input soundwaves (w/ different colours for plotting)
channels_to_test = {1000: [(0.92,'b'),(0.06,'r'),(19.80,'g')],   
                     500: [(0.48,'b'),(0.06,'r'),(19.80,'g')]#,
                     #100: [(0.45,'b'),(0.06,'r'),(19.80,'g')]#,
                     #50: [(0.42,'b'),(0.06,'r'),(19.80,'g')],  #not including simulations for 50 or 1 channels
                     #1: [(0.45,'b'),(0.06,'r'),(19.80,'g')]
                     }


i=0 #this helps with parsing through the plot columns

#run simulations and graph everything
for channel, input_intensities in channels_to_test.items():

    for intensity, colour in input_intensities:

        #run zap() -> get output 
        (Vt_arr, extra_conc_arr, light_arr) = zap(intensity, channel)

        #get M13's volatge (for graph labels)
        v = bacteriophage.getVoltage()*1000

        #format label
        lbl = f'{intensity}mW/cm\u00b2 ({v:1.2f}mV)'

        #plot graphs
        axs[0][i].plot(S_array, Vt_arr, color=colour, label=lbl)
        axs[1][i].plot(S_array, extra_conc_arr, color=colour, label=lbl)
        axs[2][i].plot(S_array, light_arr, color=colour, label=lbl)

    #add legend
    axs[0][i].legend(loc="upper right")
    axs[1][i].legend(loc="upper right")
    axs[2][i].legend(loc="upper right")

    #format axis labels
    axs[0][i].set(ylabel='mV')
    axs[1][i].set(ylabel='mM')
    axs[2][i].set(ylabel='mW/mm\u00b2')
    axs[2][i].set(xlabel='time (μs)')

    #move to next plot column
    i+=1



#annotate columns by the number of channels
col_chan = list(channels_to_test.keys())
cols = ['{} channels'.format(col) for col in col_chan]

for ax, col in zip(axs[0], cols[:3]):
    ax.annotate(col, 
                xy=(0.5, 1), 
                xytext=(0, 11),
                xycoords='axes fraction', 
                textcoords='offset points',
                size=20, 
                ha='center', 
                va='baseline')
    
    
#annotate rows by their definitions (see "rows" array)
rows = ['Membrane potential (mV)\n', 'Extracellular concentration (mM)\n', 'Blue light intensity\n']
for ax, row in zip(axs[:,0], rows):
    ax.annotate(row, 
                xy=(0, 0.5), 
                xytext=(-ax.yaxis.labelpad - 7, 0),
                xycoords=ax.yaxis.label, 
                textcoords='offset points',
                size=13, 
                ha='right', 
                va='center', 
                rotation=90)


plt.show()

