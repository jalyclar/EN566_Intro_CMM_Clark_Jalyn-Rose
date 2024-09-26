#Carbon Dating Code 
import numpy as np
import matplotlib.pyplot as plt
import math 

#Variables 
n_a=6.022*10**23 #avagradros number (nuclei/mol)
mm=14 #molar mass of carbon (g/mol)
m_i=(10**-12)*1000 #intial mass of sample (g)
n0=m_i*(1/mm)*n_a #number of nuclei before decaying begins 

half_life=5700 #beta decay half life of sample (years)
tau=half_life/np.log(2) #the time constant for the decay 

total_t=20000 #total duration of decay (years)
delta_t=10 #time step (years) 

#Parameters
t=np.arange(0,total_t,delta_t) #time step array
n=np.zeros(len(t))
a=np.zeros(len(t)) 
n[0]=n0 #First boundary
a[0]=np.divide(n0,tau) #First value 

#Euler approx and activity calculations
for i in range(0, len(t)-1):
    n[i+1]=n[i]-(np.divide(n[i],tau))*delta_t
    a[i+1]=np.divide(n[i+1],tau)

#Exact calculations 
a_exact=np.zeros(len(t))
a_exact=(np.divide(n0,tau))*np.exp(np.divide(-t,tau))

#Plot 
plt.figure()
plt.plot(t,a,'b',label='Euler Approx')
plt.plot(t,a_exact,'r',alpha=0.7,label='Exact') #Plot for analytical solution 
plt.title(r'Eluer and Analytical Solution for Activity of Carbon Isotope $^{14}_{2}\text{C}$')
plt.xlabel('Time (years)')
plt.ylabel('Activity (decay/year)')
plt.grid()
plt.legend()
plt.text(-2,1*10**9,f'Time step={delta_t} years',fontsize=10,bbox=dict(facecolor='grey',alpha=0.5,edgecolor='black'))
plt.show()

#Analysis 

#Percent deviation from exact after 2 half lives (For when delta t is 1000 years part c of question)
if delta_t >= 1000:
    two_half=2*half_life
    rounded=math.ceil(two_half/1000)*1000
    time_i=np.where(t==rounded)[0][0]
    a=a[time_i]
    a_exact=a_exact[time_i]
    perc_dev=np.divide((a-a_exact),a_exact)*100 #Should this be absolute value or no so negative shows its decrease
    print(f'After 2 half-lives the percent deviation from the exact result is {perc_dev:.2f}%')