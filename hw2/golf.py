#!/opt/anaconda3/bin/python 

import math 
import numpy as np
import matplotlib.pyplot as plt
import sys 
from matplotlib.lines import Line2D


########### Things to do:   
########### -Add angle(s) plotted in title
########### -Keep track/print total time of flight 


class golfball:

    def __init__(self):
        self.m=46 #Mass of golf ball (grams)
        self.theta=None #Angle of the trajectory of the golf ball (degrees) [45,30,15,9]
        self.v_0=70 #Initial velocity of golf ball (m/s)
        self.rho=1.29 #Density of air at sea level (kg/m^3) 
        self.a=0.0014 #Frontal area of the golf ball (m^2)
        self.g=9.8 #Gravity (m/s^2)
        self.delta_t=0.5 #time step (seconds)
        
    def parse_argv(self):
        for arg in sys.argv[1:]:
            self.theta=[float(angle) for angle in arg.split('=')[1].split(',')]
        self.trajecs={angle:{} for angle in self.theta} #Dataframe to store a dictionary for each angle for each case 

    def ideal(self): #Case 1
        for angle in self.theta: 
            #Parameters 
            x=0.0 #Initial postion of ball in x-direction
            y=0.0 #Initial postion of ball in y-direction
            v_x=self.v_0*math.cos(math.radians(angle)) #Initial velocity in x-direction
            v_y=self.v_0*math.sin(math.radians(angle)) #Initial velocity in y-direction 

            #List to store positions 
            self.x_ideal=[]
            self.y_ideal=[]

            while y>=0.0: #while the golf ball is in the air
                self.x_ideal.append(x)
                self.y_ideal.append(y)

                x=x+v_x*self.delta_t
                y=y+v_y*self.delta_t
                v_mag=math.sqrt(v_x**2+v_y**2)
                v_x=v_x 
                v_y=v_y-self.g*self.delta_t

            self.trajecs[angle]['Ideal']=(self.x_ideal,self.y_ideal) 

    #Case 2: Smooth golf ball with drag 


    #Case 3: Dimpled golf ball with drag


    #Case 4: Dimpled golf ball with drage and spin 


    #Plotting trajectories
    def plot(self): 
        plt.figure()
        for angle, cases in self.trajecs.items():
            for case, (x_values, y_values) in cases.items():
                plt.plot(x_values,y_values,label=f'{angle} degrees - {case}')
        plt.title('Trajectory of a Golf Ball') 
        plt.xlabel('Distance in x (m)')
        plt.ylabel('Height (m)')
        plt.legend()
        plt.show()

    #Analysis (Final velocities, time of flight)

if __name__ == "__main__":
    data=golfball()
    data.parse_argv()
    data.ideal()
    data.plot() 