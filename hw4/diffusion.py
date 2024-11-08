#!/opt/anaconda3/bin/python 

import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

D=2 #diffusion constant 
L=10 #grid size 
t=1 #total time 
nx=100 #gridpoints 
nt=500 #time steps
dx=L/nx #grid spacing 
dt=t/nt 
x=np.linspace(-L/2,L/2,nx)

rho=np.zeros(nx)
center1=nx//2-2
center2=nx//2+3
rho[center1:center2]=1 #inital peak around x=0

snaps=[] #to store snapshots at specific times 
t_snaps=[]
rho_new=rho.copy()
for n in range(nt): #for snaps 
    for i in range(1,nx-1): #through grid points (x vlaues)
        rho_new[i]=rho[i]+D*dt/dx**2*(rho[i+1]-2*rho[i]+rho[i-1])
    rho=rho_new.copy()
    if n%(nt//5)==0: #Takes snaps at intervals of 5
        snaps.append(np.copy(rho_new))
        t_snaps.append(n*dt)

def gauss_fit(x,amplitude,sigma): #parameters to fit gaussian 
    return amplitude*np.exp(-x**2/(2*sigma**2)) #Approximating gaussian for snapshot

plt.figure(figsize=(9,6))
for time,snap in zip(t_snaps,snaps): 
    plt.plot(x,snap,label=f'Time at {time:.2f} seconds') #for each snapshot 
    if time>0.1: #moving to later times 
        popt,_=curve_fit(gauss_fit,x,snap,p0=[1,np.sqrt(2*D*time)]) #getting fitted values of parameters with inital guess
        plt.plot(x,gauss_fit(x,*popt),linestyle='--',label=f'Gaussian fit at {time:.2f} seconds') #passing args to gauss
plt.legend()
plt.xlabel('x')
plt.ylabel(r'Density $\rho(x,t)$')
plt.show()

