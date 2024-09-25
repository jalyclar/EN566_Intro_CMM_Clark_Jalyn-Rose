#Carbon Dating Code 

import numpy as np
import matplotlib.pyplot as plt

#Parameters 
n_a=6.022*10**23 #avagradros number (nuclei/mol)
mm=14 #molar mass of carbon (g/mol)
m_i=(10**-12)*1000 #intial mass of sample (g)
n0=m_i*(1/mm)*n_a #number of nuclei before decaying begins 

half_life=5700 #beta decay half life of sample (years)
tau=half_life/np.log(2) #the time constant for the decay 

total_t=20000 #total duration of decay (years)
delta_t=10 #time step (years) 
t=np.arange(0,total_t,delta_t) #time step array 

#Euler approx calculations
n=np.zeros(len(t)) 
n[0]=n0 #First boundary 
for i in range(0, len(t)-1):
    n[i+1]=n[i]-(np.divide(n[i],tau))*delta_t

#Activity caluclations 


#Plot 
