from ast import literal_eval
from photon_emission import PhotonEmission
from ser_3D import SmoothEndoRet3D
from calcium import CalciumCluster
from M13_3D import M13Phage3D
from vpython import scene, attach_light, vector, arrow, color, box, curve, text, sphere, button, slider, wtext, sqrt, graph, rate, gcurve
from random import uniform


###############################################
####        SET PARAMETERS/VARIABLES       ####
###############################################
# get number of channels
num_channels = int(input("Enter amount of channels: "))

# get initial concentration of ca
input_Pout_str = input("Enter initial amount of calcium ions (in picomoles): ")
input_Pout_float = literal_eval(input_Pout_str)*1E-12

#sphere radius and scales
bubble_radius=1000
calcium_radius=0.23
calcium_scale=20
aequorin_radius=10
aequorin_scale=1.5

#time increments dt
dt = 1E-9 #ns

#init voltage to the ser
Vb = 0.04 #mV


###############################################
####              SET SCENE                ####
###############################################

scene.title = "Bio-Neural Dust System Simulation \n\n"
scene.background=color.black
scene.width=1380
scene.height=860
scene.align="left"
running = False
scene.select()



###############################################
####   FUNCTIONS CALLED BY CONTROL PANEL   ####
###############################################

def togglePlay(b):
    global running
    running = not running
    if running: 
        b.text = "Stop"
    else: 
        b.text = "Start"
 
def toggleVoltage(b):
    if b.text[:4]=="Stop":
        ser.updateVoltage(0)
        b.text = "Apply voltage"
    else: 
        ser.updateVoltage(0.04)
        b.text = "Stop voltage"
        
def updateSoundWaveSlider(s):
    
    #set new input soundwave
    bacteriophage.setUltrasound(s.value, num_channels)
    #get new output voltage
    new_V = bacteriophage.getVoltage()

    #set ser's new input voltage
    ser.updateVoltage(new_V)

    #update wtext
    wt_intensity.text = '{:1.2f}'.format(s.value)
    wt_m13_voltage.text = '{:1.2f}'.format(new_V*1000)
    
#these functions are no longer used but might be useful in the future
"""        
def setBubbleDiameter(s):
    wt.text = '{}'.format(s.value)
    nano_bubble.radius=(s.value/2)*1000
    updateBackground(nano_bubble.radius)
    global bubble_radius 
    bubble_radius=nano_bubble.radius
"""

###############################################
####       SET CONTROL PANEL OBJECTS       ####
###############################################

#for pausing and resuming animation in real-time
button(bind=togglePlay, text="Start", pos=scene.title_anchor)

#for pausing and resuming the input voltage
button(bind=toggleVoltage, text="Stop voltage", pos=scene.title_anchor)


scene.caption = ("<b>    -> Ultrasonic wave input:<b> \n") 
scene.append_to_caption('        Frequency: ') 
frequency_text=1.2965
wt_frequency = wtext(text='{}'.format(frequency_text))
scene.append_to_caption('MHz\n')

scene.append_to_caption('        Intensity: ') 
intensity_slider= slider(step=0.01, 
                       left=20,
                       min=0.05, 
                       max=20, 
                       value=0.9245, 
                       length=250, 
                       bind=updateSoundWaveSlider)
wt_intensity = wtext(text='{}'.format(intensity_slider.value))
scene.append_to_caption('mW/cm\u00b2\n')

scene.append_to_caption('<b>    -> M13 Bacteriophage:<b> \n') 
scene.append_to_caption('        Input voltage: ') 
m13_voltage=0.04
wt_m13_voltage = wtext(text='{:1.2f}'.format(m13_voltage*1000))
scene.append_to_caption('mV\n')

#scene.append_to_caption('=> Diffusion process: \n\n') 
scene.append_to_caption(' <b>   -> Smooth Endoplasmic Reticulum: <b>\n') 
scene.append_to_caption('        Membrane potential: ') 
test_length2=-70
wt_wavelength6 = wtext(text='{}'.format(test_length2))
scene.append_to_caption('mV\n')

scene.append_to_caption('        Channels: ') 
test_intensity2=num_channels
wt_wavelength7 = wtext(text='{}'.format(test_intensity2))
scene.append_to_caption('\n')
  
scene.append_to_caption('<b>    -> Calcium ions<b> \n') 
scene.append_to_caption('        Total concentration: ') 
calc=0
wt_wavelength8 = wtext(text='{:1.3f}'.format(calc*1E12))
scene.append_to_caption('pM\n')

scene.append_to_caption('        Inside SER: ')   
wt_cs4 = wtext(text='{:1.3f}'.format(0*1E12))
scene.append_to_caption('pM\n')
scene.append_to_caption('        Outside SER: ') 
wt_cs5 = wtext(text='{}'.format(0*1E12))
scene.append_to_caption('pM\n')


    
###############################################
####  FUNCTIONS CALLED DURING SIMULATION   ####
###############################################    

#function to invert the velocity of a molecule (used to make molecules bounce
#backwards when hitting a wall)
def invertVelocity(vel, weight):
   
    #change vector to opposite direction
    x_vel=vel.x*(-1)*weight
    y_vel=vel.y*(-1)*weight
    z_vel=vel.z*(-1)*weight
    
    #return new velocity
    return vector(x_vel,y_vel,z_vel)

#function to give a random velocity to a molecule (used to emulate a pseudo-
#brownian motion) (WIP)
def randomVelocity(vel_range):
   
    #change vector to opposite direction
    x_vel=weight(uniform(-vel_range,vel_range))
    y_vel=weight(uniform(-vel_range,vel_range))
    z_vel=weight(uniform(-vel_range,vel_range))
    
    #return new velocity
    return vector(x_vel,y_vel,z_vel)
        
#function to add a weight of 50 to the positive or negative number
#used when determining the random velocity each frame
def weight(x):
    if x < 0:
        return x - 50
    else:
        return x + 50

#function to determine whether a molecule is inside the nano-bubble or not
def isInBubble(calcium_coord,r):
    
    x=calcium_coord.x
    y=calcium_coord.y
    z=calcium_coord.z
    
    #if point is at edge of nano bubble
    try: 
        
        w=(pow(x+calcium_radius,2))+(pow(y+calcium_radius,2))+(pow(z+calcium_radius,2))
        
        if sqrt(w) <= bubble_radius:
            return True
        else: 
            return False
        
    except OverflowError:
        print("numbers too big plz stop")
        print(f'x: {x}\ny: {y}\nz: {z}\n')


#function to get a random position inside SER
def getRandInsideSER():
    while True:    
        x=uniform(3,298)
        y=uniform(-40,40)
        z=uniform(-40,40)
        xyz = vector(x,y,z)
        
        if (ser.isInside(xyz)):
            return xyz


#function to get a random position outside SER
def getRandOutsideSER():
    while True:    
        x=uniform(-700,700)
        y=uniform(-700,700)
        z=uniform(-700,700)
        xyz = vector(x,y,z)
    
        if (not ser.isInside(xyz)) and (isInBubble(xyz, calcium_radius)):
            return xyz


###############################################
####    DRAW AXIS AND BACKGROUND WALLS     ####
###############################################
 
orgn=(-1)*bubble_radius
orgn_vec=vector(orgn,orgn,orgn)

arw_x=arrow(pos=orgn_vec, axis=vector((bubble_radius*2)+300,0,0), color=color.red, shaftwidth=5, headwidth=20, round=True)
wall_x=box(pos=vector(0, 0,orgn), height=bubble_radius*2, width=6, axis=vector(bubble_radius*2,0,0), color=vector(1,0,0))
     
arw_y=arrow(pos=orgn_vec, axis=vector(0,(bubble_radius*2)+300,0), color=color.green, shaftwidth=5, headwidth=20, round=True)
wall_y=box(pos=vector(orgn, 0,0), width=bubble_radius*2, height=6, axis=vector(0,bubble_radius*2,0), color=vector(0,1,0))
arw_z=arrow(pos=orgn_vec, axis=vector(0,0,(bubble_radius*2)+300), color=color.blue, shaftwidth=5, headwidth=20, round=True)
wall_z=box(pos=vector(0,orgn,0), width=bubble_radius*2, length=bubble_radius*2, height=6, color=vector(0,0,1))

#function called by slider to scale background walls according to bubble radius
def updateBackground(rad):
    
    orgn=(-1)*rad
    orgn_vec=vector(orgn,orgn,orgn)
    
    arw_x.pos=orgn_vec
    arw_x.axis=vector((rad*2)+300,0,0)
    wall_x.pos=vector(0, 0,orgn)
    wall_x.axis=vector(rad*2,0,0)
    wall_x.height=rad*2
    
    arw_z.pos=orgn_vec
    arw_z.axis=vector(0,0,(rad*2)+300)
    wall_z.pos=vector(0, orgn,0)
    wall_z.axis=vector(rad*2,0,0)
    wall_z.width=rad*2
    wall_z.length=rad*2
    
    arw_y.pos=orgn_vec
    arw_y.axis=vector(0,(rad*2)+300,0)
    wall_y.pos=vector(orgn,0,0)
    wall_y.axis=vector(0,rad*2,0)
    wall_y.width=rad*2


###############################################
####     DRAW OBJECTS FOR SCALE DEMO       ####
###############################################

#300nm scale objects
arw_ser=curve(pos=[vector(0,-75,0), vector(300,-75,0)],
              radius=2, 
              color=color.red, 
              round=True)

arw_ser1=curve(pos=[vector(0,-85,0), vector(0,-65,0)],
               radius=2, 
               color=color.red, 
               round=True)

arw_ser1=curve(pos=[vector(300,-85,0), vector(300,-65,0)],
               radius=2, 
               color=color.red, 
               round=True)

ser_text=text(text='300nm',
              align='center',
              pos=vector(150,-100,0), 
              height=25)

#1mm scale objects
#arw_mm=curve(pos=[vector(-38000,-500000,0),vector(-38000,500000,0)],radius=10000, color=color.red, round=True)
#arw_mm=curve(pos=[vector(-80000,-500000,0),vector(10000,-500000,0)],radius=10000, color=color.red, round=True)
#arw_mm=curve(pos=[vector(-80000,500000,0),vector(10000,500000,0)],radius=10000, color=color.red, round=True)
#mm_text=text(text='1mm',pos=vector(-30000,-50000,0), height=40200)

#credit card object
#credit_card=box(pos=vector(-42880000, 0,0), height=760000, width=53980000, length=85600000, axis=vector(-1,0,0), color=color.gray(0.8))
#cc_tex2=text(text='Credit card',pos=vector(-83000000,550000,-15000000), up=vector(0,0,-1),height=6000000, color=color.gray(0.5))
#cc_tex3=text(text='1234 5678 9012 3456',pos=vector(-83000000,550000,8000000), up=vector(0,0,-1),height=6000000, color=color.gray(0.5))
#cc_tex3=text(text='Optimus Prime',pos=vector(-83000000,550000,17500000), up=vector(0,0,-1),height=6000000, color=color.gray(0.5))

#2μm scale object
arw_bubble=curve(pos=[vector(bubble_radius+30,-1000,0),
                      vector(bubble_radius+30,1000,0)],
                 radius=10,
                 color=color.red, 
                 round=True)

arw_bubble1=curve(pos=[vector(bubble_radius-20,-1000,0),
                       vector(bubble_radius+80,-1000,0)],
                  radius=10, 
                  color=color.red, 
                  round=True)

arw_bubble2=curve(pos=[vector(bubble_radius-20,1000,0),
                       vector(bubble_radius+80,1000,0)],
                  radius=10, 
                  color=color.red, 
                  round=True)

bubble_text=text(text='2μm \n(0.002mm)',
                 pos=vector(bubble_radius+40,150,0), 
                 height=300)


###############################################
####           DRAW NANO-BUBBLE            ####
###############################################

#radius = 1000 (or slider value)
#opacity = 65% transparent
nano_bubble=sphere(color=vector(0.5,0.7,0.9), #blue
                   pos=vector(0,0,0),      
                   radius=bubble_radius,
                   opacity=0.35, 
                   shininess=5)


###############################################
####  DRAW SMOOTH ENDOPLASMIC RETICULUM    ####
###############################################

# membrane resting potential and max potential
Vrest = -0.050 #-50mV
Vmax = -0.010  #-10mV

#instantiate class
ser = SmoothEndoRet3D(50,5, 300, num_channels, Vrest, Vmax, dt)

#start with a voltage Vb
ser.updateVoltage(Vb)



        
###############################################
####          DRAW CALCIUM IONS            ####
###############################################

#total calcium spheres
num_calcium=10000
calcium_cluster={}  #save all calcium objects to this array

#instantiate model class
calcium_model = CalciumCluster(Vrest,Vmax,input_Pout_float)

#set initial ratio and extra/intracellular concentrations
(init_ratio, init_Pin, init_Pout) = calcium_model.updateConcentrations(Vrest)

#update texts in control panel
wt_wavelength8.text = '{:1.3f}'.format((init_Pin+init_Pout)*1E12)
wt_cs4.text = '{:1.3f}'.format(init_Pin*1E12)
wt_cs5.text = '{:1.3f}'.format(init_Pout*1E12)

#get the max extracellular concentration of ion (i.e, when Vt = Vmax)
Pout_max = calcium_model.getPoutMax()


#equivalent concentration for 3d calcium spheres (1000 total)
Pin3 = round(num_calcium/(init_ratio + 1))
Pout3 = num_calcium - Pin3

#variables to keep track of the numbers during loop
num_in_ser=Pin3
num_out_ser=Pout3
external_indexes = []
internal_indexes = []

for i in range(Pin3):
    
    xyz_ve = getRandInsideSER()

    #append calcium ion to calcium cluster array 
    calcium_cluster[i]=(sphere(color=vector(0,1,0), #green
                               pos=xyz_ve,      
                               radius=calcium_radius*5, #30nm
                               opacity=1, 
                               shininess=1,
                               velocity=vector(0,0,0),
                               in_ser=True, #created this class variable to track if a molecule is inside SER
                               )
                        )  
            
    internal_indexes.append(i)
    

for j in range(Pout3):

    xyz = getRandOutsideSER()
    
    #append calcium ion to calcium cluster array 
    calcium_cluster[Pin3+j]=(sphere(color=vector(0,1,0), #green
                                  pos=xyz,      
                                  radius=calcium_radius*calcium_scale, #30nm
                                  opacity=1, 
                                  shininess=1,
                                  velocity=vector(uniform(-1,1), 
                                                  uniform(-1,1), 
                                                  uniform(-1,1)),
                                  in_ser=False, #created this class variable to track if a molecule is inside SER
                                  )
                           )  
    
    external_indexes.append(Pin3+j)
    
#these stay in memory so might as well rm them
del xyz_ve, xyz



###############################################
####       DRAW AEQUORIN MOLECULES         ####
###############################################

num_aequorin=30
aequorin_cluster = [] #save all aequorin objects in this array

#initial brightness of each aequorin molecule
init_brightness = round(init_Pout/Pout_max,2)

for _ in range(num_aequorin):
    
    xyz_a = getRandOutsideSER()

    aequ = sphere(color=vector(0.5,0.7,0.9),
                    pos=xyz_a,      
                    radius=aequorin_radius*aequorin_scale,
                    opacity=0.8, 
                    velocity=vector(uniform(-1,1), 
                                    uniform(-1,1),  #start with a random velocity between -1 and 1
                                    uniform(-1,1)),
                    in_ser=False, #created this class variable to track if a molecule is inside SER
                    shininess=5,
                    make_trail=True,
                    trail_type="curve",
                    retain=4,
                    trail_radius=0.1)
    
    blue_light=attach_light(aequ, offset=vector(0,0,-15), color=vector(0,0,init_brightness))
    aequorin_cluster.append((aequ,blue_light))

#these stay in memory so might as well rm them
del xyz_a   
    
#setup init photon emission variables
diffusion_coefficient = 1E-8
aequorin_radius = 1E-8
photon_emiter = PhotonEmission(init_Pout, Pout_max, diffusion_coefficient, aequorin_radius, dt)
init_output_intensity = photon_emiter.getLightIntensity(init_Pout,dt)



###############################################
####        DRAW M13 BACTERIOPHAGE         ####
###############################################

#file path of the excel provided by Arash
filePath = 'M13_Voltage_V3.xlsx'

#(M13 virus - will respond to ultrasonic wave and emit current through nanowire)
bacteriophage = M13Phage3D(filePath, 
                           radius=6, 
                           length=880, 
                           position=vector(-600,0,-300))



###############################################
####             SIMULATION LOOP           ####
###############################################   

#this will be called by rate() to cap animation to 60 loops/sec
loop_rate=150

#variable to track the number of iterations and time
loop_num = 0
time = 0
time_step = dt*1E9

#graph the number of calciums spheres inside the SER
g = graph(align = 'right', title='<b>Extracellular concentration</b>', width = 500, height = 200,
      xtitle='<i>nanoseconds</i>', ytitle='<i>picomoles</i>', ymin=0, ymax=Pout_max*1E12+1)
f1 = gcurve(color=color.cyan)

#graph the SER membrane potential over-time
g2 = graph(align = 'right', title='<b>Membrane potential</b>', width = 500, height = 200,
      xtitle='<i>nanoseconds</i>', ytitle='<i>millivolts<i>', ymin=-70, ymax= 0)
f12 = gcurve(color=color.cyan)

#graph the light intensity over-time
g3 = graph(align = 'right', title='<b>Blue light intensity</b>', width = 500, height = 200,
      xtitle='<i>nanoseconds</i>', ytitle='<i>mW/mm\u00b2<i>')
f13 = gcurve(color=color.cyan)
     
f12.plot([time, Vrest*1000])
f1.plot([time, init_Pout*1E12])
f13.plot([time, init_output_intensity/1000])

time = time + time_step

#this loop will check every molecule inside the nano-bubble and update
#their velocities according to certain conditions.
while True:
  
  if running:
      
      rate(loop_rate)
      
      #get the SER membrane's voltage for this time-step
      Vt = ser.nextTimeStep()

      #update/get real concentrations according to new voltage
      ratio, Pin_t, Pout_t = calcium_model.updateConcentrations(Vt)
      
      #equivalent concentration for 3d calcium spheres 
      Pin_3d = round((num_calcium)/(ratio + 1))
      Pout_3d = (num_calcium) - Pin_3d
      
      #if Pout_3d is larger than num_out_ser -> need to diffuse some
      diffusion_mode = (Pout_3d > num_out_ser) 
      
      #if Pout_3d is lower than num_out_ser -> need to absorb some
      absorbtion_mode = (Pout_3d < num_out_ser)
      
      wt_cs4.text = '{:1.3f}'.format(Pin_t*1E12)
      wt_cs5.text = '{:1.3f}'.format(Pout_t*1E12)
      wt_wavelength6.text = '{:1.3f}'.format(Vt*1000)

      intens = round(Pout_t/Pout_max,2)      
      
      #get output blue light intensity
      output_intensity = photon_emiter.getLightIntensity(Pout_t, dt)

      #iterate through every extracellular caclium ion
      if absorbtion_mode:
          num_to_absorb = num_out_ser - Pout_3d
          
          for i in range(num_to_absorb):
              temp = external_indexes.pop()
              
              calcium_ion = calcium_cluster[temp]
              
              ### set pos inside SER
              xyz = getRandInsideSER()
              
              calcium_ion.pos = xyz
              calcium_ion.in_ser=True
              
              calcium_ion.radius=calcium_radius
              calcium_ion.color=vector(1,0,0)
              
              calcium_ion.velocity=vector(0,0,0)


              num_out_ser-=1
              num_in_ser+=1
              
              internal_indexes.append(temp)
              
              
      elif diffusion_mode:
          
          num_to_diffuse = Pout_3d - num_out_ser
          cutter = -num_to_diffuse
          
          for i in range(num_to_diffuse):
              temp = internal_indexes.pop()
              
              calcium_ion = calcium_cluster[temp]
              
              ### set pos outside SER
              xyz = getRandOutsideSER()
              
              calcium_ion.pos = xyz
              calcium_ion.in_ser=False
             
              calcium_ion.radius=calcium_radius*calcium_scale
              calcium_ion.color = vector(1,1,0) #yellow
 
              num_out_ser+=1
              num_in_ser-=1
              
              external_indexes.append(temp)
          
                
      for i in external_indexes:
          
          calcium_ion = calcium_cluster[i]
          
          #if Ca is in bubble 
          if isInBubble(calcium_ion.pos, calcium_radius):
              
              #if Ca is in SER
              if ser.isInside(calcium_ion.pos):
                  #set new velocity  to opposite direction (aka bounce back outside)
                  calcium_ion.velocity = invertVelocity(calcium_ion.velocity, 1.2)
              else:
                  #set new random velocity
                  calcium_ion.velocity = randomVelocity(100)
                      
          else:
             #set new velocity  to opposite direction (aka bounce back inside)
             calcium_ion.velocity = invertVelocity(calcium_ion.velocity, 1.2)

          #move according to new velocity 
          calcium_ion.pos= calcium_ion.pos + calcium_ion.velocity*0.01
 
        
      #iterate through every aequorin molecule
      for aequorin in aequorin_cluster:
          
          if isInBubble(aequorin[0].pos, 15):
              
              #if is in SER
              if ser.isInside(aequorin[0].pos):
                  aequorin[0].velocity = invertVelocity(aequorin[0].velocity, 1.2)
              else:
                  aequorin[0].velocity = randomVelocity(100)
          else:
              aequorin[0].velocity = invertVelocity(aequorin[0].velocity, 1.2)
              
          aequorin[0].pos = aequorin[0].pos + aequorin[0].velocity*0.07
          aequorin[1].color = vector(0,0,intens)


      f12.plot([time, Vt*1000])
      f1.plot([time, Pout_t*1E12])
      f13.plot([time, output_intensity/1000])
      time = time + time_step
      loop_num=loop_num+1