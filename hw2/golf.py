#!/opt/anaconda3/bin/python 

import math 
import numpy as np
import matplotlib.pyplot as plt
import sys 
import pandas as pd

############                   TO DO:
############              -Optimize each function to shorten amount of lines 
############              -Document code 
############              -Add DEBUG=TRUE print statments throughout 
############              -Adjust sig figs in printing results statments 

class golfball:

    def __init__(self):
        self.m=46 #Mass of golf ball (grams)
        self.theta=None #Angle of the trajectory of the golf ball (degrees) [45,30,15,9]
        self.v_0=70 #Initial velocity of golf ball (m/s)
        self.rho=1.29*1000 #Density of air at sea level (g/m^3) 
        self.a=0.0014 #Frontal area of the golf ball (m^2)
        self.g=9.8 #Gravity (m/s^2)
        self.delta_t=0.1 #time step (seconds) 
        
    def parse_argv(self):
        for arg in sys.argv[1:]:
            self.theta=float(arg.split('=')[1])
        self.v_x_0=self.v_0*math.cos(math.radians(self.theta)) #Initial velocity in x-direction
        self.v_y_0=self.v_0*math.sin(math.radians(self.theta)) #Initial velocity in y-direction

    def ideal(self): #Case 1: Ideal
        #Parameters 
        x=0.0 #Initial postion of ball in x-direction
        y=0.0 #Initial postion of ball in y-direction
        v_x=self.v_x_0 #Initial velocity in x-direction
        v_y=self.v_y_0 #Initial velocity in y-direction
        #Lists to store positions 
        self.x_ideal=[]
        self.y_ideal=[]
        ideal_total_t=0 #Initializing time 
        while y>=0.0: #while the golf ball is in the air
            self.x_ideal.append(x)
            self.y_ideal.append(y)
            x=x+v_x*self.delta_t
            y=y+v_y*self.delta_t
            v_mag=math.sqrt(v_x**2+v_y**2) 
            v_x=v_x 
            v_y=v_y-self.g*self.delta_t
            ideal_total_t+=self.delta_t #Time of flight 
        vf_ideal=v_mag
        print(f'Ideal golf ball hit at {self.theta}: v_f (m/s)={vf_ideal} t_f (s)={ideal_total_t}')
 
    def smooth_drag(self): #Case 2: Smooth golf ball with drag
        #Parameters 
        x=0.0 #Initial postion of ball in x-direction
        y=0.0 #Initial postion of ball in y-direction
        v_x=self.v_x_0 #Initial velocity in x-direction
        v_y=self.v_y_0 #Initial velocity in y-direction
        self.x_smooth_drag=[]
        self.y_smooth_drag=[]
        C=0.5
        smooth_drag_t=0
        while y>=0.0:
            self.x_smooth_drag.append(x)
            self.y_smooth_drag.append(y)
            x=x+v_x*self.delta_t
            y=y+v_y*self.delta_t
            v_mag=math.sqrt(v_x**2+v_y**2)
            v_x=v_x-((C*self.rho*self.a*v_mag)/self.m)*v_x*self.delta_t
            v_y=v_y-(self.g*self.delta_t)-((C*self.rho*self.a*v_mag)/self.m)*v_y*self.delta_t
            smooth_drag_t+=self.delta_t
        vf_smooth_drag=v_mag
        print(f'Smooth golf ball with drag hit at {self.theta}: v_f (m/s)={vf_smooth_drag} t_f (s)={smooth_drag_t}')

    def dimpled_drag(self): #Case 3: Dimpled golf ball with drag
        #Parameters 
        x=0.0 #Initial postion of ball in x-direction
        y=0.0 #Initial postion of ball in y-direction
        v_x=self.v_x_0 #Initial velocity in x-direction
        v_y=self.v_y_0 #Initial velocity in y-direction
        self.x_dimple_drag=[]
        self.y_dimple_drag=[]
        dimple_drag_t=0
        while y>=0:
            self.x_dimple_drag.append(x)
            self.y_dimple_drag.append(y)
            x=x+v_x*self.delta_t
            y=y+v_y*self.delta_t
            v_mag=math.sqrt(v_x**2+v_y**2)
            if v_mag<=14:
                C=0.5
            elif v_mag>14:
                C=7/v_mag
            v_x=v_x-((C*self.rho*self.a*v_mag)/self.m)*v_x*self.delta_t
            v_y=v_y-(self.g*self.delta_t)-((C*self.rho*self.a*v_mag)/self.m)*v_y*self.delta_t
            dimple_drag_t+=self.delta_t
        vf_dimple_drag=v_mag
        print(f'Dimpled golf ball with drag hit at {self.theta}: v_f (m/s)={vf_dimple_drag} t_f (s)={dimple_drag_t}')

    def spin(self): #Case 4: Dimpled golf ball with drage and spin
        #Parameters 
        x=0.0 #Initial postion of ball in x-direction
        y=0.0 #Initial postion of ball in y-direction
        v_x=self.v_x_0 #Initial velocity in x-direction
        v_y=self.v_y_0 #Initial velocity in y-direction
        Sow=0.25 #s^-1
        self.x_spin=[]
        self.y_spin=[]
        spin_t=0
        while y>=0:
            self.x_spin.append(x)
            self.y_spin.append(y)
            x=x+v_x*self.delta_t
            y=y+v_y*self.delta_t
            v_mag=math.sqrt(v_x**2+v_y**2)
            if v_mag<=14:
                C=0.5
            elif v_mag>14:
                C=7/v_mag
            v_x=v_x+(-((C*self.rho*self.a*v_mag)/self.m)-(Sow/self.m)*v_y)*self.delta_t
            v_y=v_y+(-((C*self.rho*self.a*v_mag)/self.m)-(Sow/self.m)*v_x-self.g)*self.delta_t
            spin_t+=self.delta_t
        vf_spin=v_mag
        print(f'Dimpled golf ball with drag and spin hit at {self.theta}: v_f (m/s)={vf_spin} t_f (s)={spin_t}')

    #Plotting trajectories
    def plot(self): 
        plt.figure()
        plt.plot(self.x_ideal,self.y_ideal,label='Ideal')
        plt.plot(self.x_smooth_drag,self.y_smooth_drag,label='Smooth ball with drag')
        plt.plot(self.x_dimple_drag,self.y_dimple_drag,label='Dimpled ball with Drag')
        plt.plot(self.x_spin,self.y_spin,label='Dimpled ball with Drag and spin')
        plt.title(fr'Trajectory of a Golf Ball at $\theta$={self.theta}') 
        plt.xlabel('Distance in x (m)') 
        plt.ylabel('Height (m)')
        plt.legend(title='Cases',loc='upper left',fontsize='small')
        plt.show()

if __name__ == "__main__":
    data=golfball()
    data.parse_argv()
    data.ideal()
    data.smooth_drag()
    data.dimpled_drag()
    data.spin()
    data.plot() 