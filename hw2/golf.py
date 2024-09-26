#Golf Ball Code 

import math 
import numpy as np

#Variables 
m=46 #Mass of golf ball (grams)
theta=[45,30,15,9] #Angle of the trajectory of the golf ball (degrees)
v_0=70 #Initial velocity of golf ball (m/s)
rho=1.29 #Density of air at sea level (kg/m^3) 
a=0.0014 #Frontal area of the golf ball (m^2)
g=9.8 #Gravity (m/s^2)

delta_t=1 #time step (seconds)

#Parameters 
x_0=0 #Initial postion of ball in x-direction
y_0=0 #Initial postion of ball in y-direction
v_x0=v_0*math.cos(theta) #Initial velocity in x-direction
v_y0=v_0*math.sin(theta) #Initial velocity in y-direction
t=0 #Starting at 0 

#Case 1: No drag or spin 


#Case 2: Smooth golf ball with drag 


#Case 3: Dimpled golf ball with drag


#Case 4: Dimpled golf ball with drage and spin 


#Plotting trajectories 


#Analysis 