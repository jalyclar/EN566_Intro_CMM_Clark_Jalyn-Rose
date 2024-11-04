import numpy as np
import matplotlib.pyplot as plt 
import random

D=2 #diffusion constant 
L=100 #grid size 
nx=20 #gridpoints 
nt=100 #time steps
dx=L/(nx-1) #grid spacing 
x=np.arange(-L/2,L/2,dx)

t=0
for i in len(nt):
    t=t+1


