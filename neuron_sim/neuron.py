#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Neuron Class contains all 3D graphics & visual simulation logic

Simulation of parameters is handled by the model classes

Created on Thu Nov  3 17:22:18 2022

@author: majawitter
"""

from vpython import cylinder, vector, sphere, compound, ellipsoid, cone, color, curve, pi, cos, sin, dot, cross, shapes, paths, extrusion, local_light
import random
import math
import numpy as np

class Neuron:
    
    def __init__(self, soma_radius, membrane_model, opsin_model):        
        self.prob_open = 0.0
        
        self.light_intensity = 0.0
        
        self.membrane_model = membrane_model
                
        self.opsin_model = opsin_model
        
        components = []
        
        myelin_list = []
        myelin_length = 20
        
        # =============================================================================
        # Creating Soma
        # =============================================================================
        self.soma = sphere(color=vector(0.95,  0.67,  0.73), 
                           pos=vector(0,0,0),
                           radius=soma_radius)
        components.append(self.soma)
        
        # =============================================================================
        # Creating Axon + Opsins on Axon
        # =============================================================================
        self.axon_connection = cone(color=vector(0.95,  0.67,  0.73), 
                                    pos=vector(soma_radius-1,0,0),
                                    axis=vector(1,0,0),
                                    radius=myelin_length/4,
                                    length=myelin_length)
        
        mi = ellipsoid(color=vector(0.99, 0.81, 0.33),
                              axis=vector(1,0,0),
                              pos=vector((self.axon_connection.pos.x+self.axon_connection.length/2+myelin_length/2),0,0),
                              length=myelin_length,
                              height=myelin_length/3,
                              width=myelin_length/3)
        
        myelin_list.append(mi)
        
        opsin_list = []
        xi = mi.pos.x - (mi.length/2)
        x_inc = mi.length/5
        
        a = mi.length/2
        b = mi.height/2
        h = mi.pos.x
        k = mi.pos.y
        
        for i in range (6):

            y = math.sqrt((b**2)*(1-(((xi-h)**2)/(a**2)))) + k
            theta = random.randint(0, 90)
            
            opsin_list.append(sphere(color=color.red,
                                     pos=vector(xi, y, 0),
                                     radius = 1))
            opsin_list.append(sphere(color=color.red,
                                     pos=vector(xi, -y, 0),
                                     radius = 1))
            xi += x_inc
        
        components.append(myelin_list[0])
        
        # next midpoint with current axis
        x1 = myelin_list[0].pos.x+myelin_length
        y1 = myelin_list[0].pos.y
        m = -0.2
        b = -1*m*x1
        x2 = x1-(myelin_length/2)
        y2 = m*x2+b
        v1 = x1-x2
        v2 = y1-y2
        u1 = v1/(math.sqrt(v1**2+v2**2))
        u2 = v2/(math.sqrt(v1**2+v2**2))
        xn = x1 - (myelin_length/2)*u1
        yn = y1 - (myelin_length/2)*u2
        # current endpoint, next start point
        xo = x1 - (myelin_length/2)
        yo = 0
        
        for i in range(1,8):
            myelin = ellipsoid(color=vector(0.99, 0.81, 0.33),
                                  axis=vector(1,m,0),
                                  pos=vector(x1-(xn-xo),y1-(yn-yo),0),
                                  length=myelin_length,
                                  height=myelin_length/3,
                                  width=myelin_length/3)
            myelin_list.append(myelin)
            
            # opsins
            x_norm = myelin.axis.norm().x
            y_norm = myelin.axis.norm().y
            num_ops = 5
            inc = myelin.length/num_ops
            
            a = myelin.length/2
            b = myelin.height/2
            h = myelin.pos.x
            k = myelin.pos.y
            
            theta = -1*(myelin.axis.norm().diff_angle(vector(1,0,0)))
            
            y_diff = (myelin.pos.y + (a*y_norm))-(myelin.pos.y - (a*y_norm))
            y_inc = y_diff/num_ops
            mult_i = -2

            for j in range (mult_i,mult_i+num_ops):
                x = myelin.pos.x + (j*inc*x_norm)
                y = myelin.pos.y +(j*y_inc)
                
                if i < 4:
                    aq = (b**2)*sin(theta)**2 + (a**2)*cos(theta)**2
                    bq = 2*(b**2)*(x-h)*cos(theta)*sin(theta) - 2*(a**2)*(x-h)*sin(theta)*cos(theta)
                    cq = (b**2)*((x-h)*cos(theta))**2 + (a**2)*((x-h)*sin(theta))**2 - ((a**2)*(b**2))
                    
                    ym1 = ((-bq + math.sqrt((bq**2)-4*aq*cq))/(2*aq)) + k
                    ym2 = ((-bq - math.sqrt((bq**2)-4*aq*cq))/(2*aq)) + k
                    
                    opsin_list.append(sphere(color=color.red,
                                             pos=vector(x,ym1,0),
                                             radius = 1))
                    opsin_list.append(sphere(color=color.red,
                                             pos=vector(x,ym2,0),
                                             radius = 1))
                else:
                    aq = (b**2)*cos(theta)**2 + (a**2)*sin(theta)**2
                    bq = 2*(b**2)*(y-k)*cos(theta)*sin(theta) - 2*(a**2)*(y-k)*sin(theta)*cos(theta)
                    cq = (b**2)*((y-k)*sin(theta))**2 + (a**2)*((y-k)*cos(theta))**2 - ((a**2)*(b**2))
                    
                    xm1 = ((-bq + math.sqrt((bq**2)-4*aq*cq))/(2*aq)) + h
                    xm2 = ((-bq - math.sqrt((bq**2)-4*aq*cq))/(2*aq)) + h
  
                    opsin_list.append(sphere(color=color.red,
                                             pos=vector(xm1,y,0),
                                             radius = 1))
                    opsin_list.append(sphere(color=color.red,
                                             pos=vector(xm2,y,0),
                                             radius = 1))
            
            components.append(myelin_list[i])
            
            # current endpoint, next start point
            xo = (x1-(xn-xo)) + (myelin_length/2)*u1
            yo = y1-(yn-yo) + (myelin_length/2)*u2
            # next midpoint with current axis
            x1 = (x1-(xn-xo)) + myelin_length*u1
            y1 = -(yn-yo) + myelin_length*u2
            m += -0.4
            b = y1 - m*x1
            x2 = x1-(myelin_length/2)
            y2 = m*x2+b
            # vector between the two points
            v1 = x1-x2
            v2 = y1-y2
            # normalized vector to scale other points/distances
            u1 = v1/(math.sqrt(v1**2+v2**2))
            u2 = v2/(math.sqrt(v1**2+v2**2))
            xn = x1 - (myelin_length/2)*u1
            yn = y1 - (myelin_length/2)*u2
            
        
        self.axon = myelin_list
                
        # =============================================================================
        # Creating Synapes with associated light visuals
        # =============================================================================
        vo = vector(xo,yo,0)
        synapse_list = []
        synapse_connection_list = []
        
        for dz in [0, 0.5, -0.8]:
            v1 = vo
            cv_points = [v1]
            y = -1
            z = 0
            if dz == 0:
                mx = -1
            else:
                mx = 1
            x = mx
            for i in range(10):
                v1 = vector(v1.x+x, v1.y+y, v1.z+z)
                y += -0.8
                x -= ((i-5)/8)*(mx)
                z += dz
                cv_points.append(v1)
            cv = curve(pos=cv_points, color=vector(0.95,  0.67,  0.73), radius=myelin_length/12)
            
            synapse = cone(color=vector(0.98, 0.11, 0.33),
                           pos=vector(v1.x+x/2, v1.y+y/2, v1.z+z/2),
                           axis=vector(-x,-y,-z),
                           radius=myelin_length/4,
                           length=myelin_length/2,
                           emissive=False)
            
            ap_light = local_light(pos=vector(v1.x+x, v1.y+y, v1.z+z),
                                   color=color.white,
                                   visible=False)
            
            synapse_connection_list.append(cv)
            synapse_list.append((synapse, ap_light))
            
        self.synapses = synapse_list
        
        self.synapse_connections = synapse_connection_list
                
        self.neuron_cell = compound(components)
        
        # =============================================================================
        # Creating Dendrites
        # =============================================================================
        z_list = [-1,-0.7,0,0.7,1]
        
        dendrite_list = []
        dend_points = set()
        branches = []
        
        for i in range(6):
            theta = i*2/6*pi
            for z in z_list:
                if (z == 0) and (i == 0):
                    continue
                x = (math.sqrt(1-(z**2))*cos(theta))
                y = (math.sqrt(1-(z**2))*sin(theta))
                if x == 0:
                    x = abs(x)
                if y == 0:
                    y = abs(y)

                if dend_points.issuperset({str((x,y,z))}):
                    continue
                dend_points.add(str((x,y,z)))
                
                r_dend = random.uniform(myelin_length/6, myelin_length/4)
                        
                dendrite_list.append(cone(color=vector(0.95,  0.67,  0.73), 
                                          pos=vector(x*(soma_radius-0.5),y*(soma_radius-0.5),z*(soma_radius-0.5)),
                                          axis=vector(x,y,z),
                                          length=r_dend*4,
                                          radius=r_dend))
                branch_list = []
                branch = curve(color=vector(0.95,  0.67,  0.73), 
                               pos=[vector(x*(soma_radius),y*(soma_radius),z*(soma_radius)),
                                    vector(x*(soma_radius+r_dend*5),y*(soma_radius+r_dend*5),z*(soma_radius+r_dend*5))],
                               radius=r_dend/2)
                branch_list.append(branch)
                
                px = x*(soma_radius+r_dend*5)
                py = y*(soma_radius+r_dend*5)
                pz = z*(soma_radius+r_dend*5)
                dist = (soma_radius+r_dend*5)
                r_dend = r_dend*0.75
                
                self.dendriteBranch(branch_list, random.randint(3, 4), vector(x,y,z), vector(x*(soma_radius),y*(soma_radius),z*(soma_radius)), vector(px, py, pz), dist+r_dend*5, r_dend/2, theta)
                    
        self.dendrites = dendrite_list
        
        # =============================================================================
        # Adding opsins on the Soma
        # =============================================================================
        z_list = [-1,-0.9,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,0.9,1]
        ops_points = set()
        
        for i in range(14):
            theta = i*2/14*pi
            for z in z_list:
                if (z == 0) and (i == 0):
                    continue
                x = (math.sqrt(1-(z**2))*cos(theta))
                y = (math.sqrt(1-(z**2))*sin(theta))
                if x == 0:
                    x = abs(x)
                if y == 0:
                    y = abs(y)

                if ops_points.issuperset({str((x,y,z))}):
                    continue
                ops_points.add(str((x,y,z)))
                
                opsin_list.append(sphere(color=color.red,
                                         pos=vector(x*(soma_radius),y*(soma_radius),z*(soma_radius)),
                                         radius = 1))
        
        self.opsins = opsin_list
        
    # =============================================================================
    # Generate a randomly branching dendrite
    # =============================================================================
    def dendriteBranch(self, branch_list, joints, parent_axis, parent_start, start_point, dist, rad, theta):
        branch = random.uniform(0, 1)
        if joints > 0:
            # main branch
            x_main = parent_axis.x + random.uniform(-0.05,0.05)
            y_main = parent_axis.y + random.uniform(-0.05,0.05)
            main_branch = curve(color=vector(0.95,  0.67,  0.73), 
                                pos=[start_point,
                                     vector(x_main*dist,y_main*dist,parent_axis.z*dist)],
                                radius=rad)
            
            branch_list.append(main_branch)
            branch_list = self.dendriteBranch(branch_list, joints-1, vector(x_main, y_main, parent_axis.z), start_point,
                            vector(x_main*dist,y_main*dist,parent_axis.z*dist),
                            dist+rad*5,rad*0.75, theta)
            
            if branch >= 0.5:
                # side branch
                x_side = random.uniform(max(-1,parent_axis.x-0.1),min(1,parent_axis.x+0.1))
                y_side = random.uniform(max(-1,parent_axis.y-0.1),min(1,parent_axis.y+0.1))
                z_side = random.uniform(max(-1,parent_axis.z-0.1),min(1,parent_axis.z+0.1))
                
                theta += random.uniform(-pi/5,pi/5)
                
                
                side_branch = curve(color=vector(0.95,  0.67,  0.73), 
                                    pos=[start_point,
                                         vector(x_side*dist,y_side*dist,z_side*dist)],
                                    radius=rad*0.8)
                
                branch_list.append(side_branch)
                return self.dendriteBranch(branch_list, joints-3, vector(x_side, y_side, z_side), start_point,
                                vector(x_side*dist,y_side*dist,z_side*dist),
                                dist+rad*5,rad*0.75,theta)
            else:
                return branch_list
                
        else:
            branch_list.append(curve(color=vector(0.95,  0.67,  0.73), 
                                     pos=[start_point,
                                          vector(parent_axis.x*dist,parent_axis.y*dist,parent_axis.z*dist)],
                                     radius=rad))
            return branch_list
    
    # =============================================================================
    # Shows a flash when an Action Potential is fired
    # =============================================================================
    def apVisual(self, on):
        for synapse in self.synapses:
            syn = synapse[0]
            light = synapse[1]
            if on:
                syn.color = color.white
                light.visible = True
            else:
                syn.color = vector(0.98, 0.11, 0.33)
                light.visible = False
    
    # =============================================================================
    # Updates opsin color to reflect the proportion in the closed (red) and open (blue) states
    # =============================================================================
    def opsinStateVisual(self):
        prob_Ostate = self.opsin_model.getModel().getO1() + self.opsin_model.getModel().getO2()
        if (self.prob_open == 0) and (prob_Ostate == 0):
            diff = 0
        else:
            diff = abs(self.prob_open - prob_Ostate)/((self.prob_open + prob_Ostate)/2)
        if diff > 0.001:
            for opsin in self.opsins:
                r = random.uniform(0,1.0)
                if r <= prob_Ostate:
                    opsin.color = color.blue
                else:
                    opsin.color = color.red
        self.prob_open = prob_Ostate

    # =============================================================================
    # Reset the model    
    # =============================================================================
    def reset(self, membrane_model, opsin_model):
        self.membrane_model = membrane_model
        self.opsin_model = opsin_model
        self.apVisual(False)
        self.opsinStateVisual()

    # =============================================================================
    # Calls on membrane model for update
    # =============================================================================
    def nextTimeStep(self, t, irradiance):
        self.opsin_model.nextTimeStep(self.membrane_model.getMembranePotential(), irradiance)
        
        self.membrane_model.updateMembranePotential(t, abs(self.opsin_model.getI()))
        self.apVisual(self.membrane_model.isSpiking())
        self.opsinStateVisual()
            
    def getMembraneModel(self):
        return self.membrane_model
    
    def getOpsinModel(self):
        return self.opsin_model
    
        
        
        
        
