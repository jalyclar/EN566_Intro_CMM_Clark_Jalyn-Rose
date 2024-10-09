import numpy as np
import matplotlib.pyplot as plt
import math

def partone(DEBUG=False): #Euler-cromer method 
   #Variables 
    g=9.8 #m/s^2
    l=9.8 #m
    gamma=0.25 #s^-1
    omega_values=[]
    theta_values=[]
    t_values=[]
    d_omega_values=[]
    d_theta_values=[]
    d_t_values=[]
    t=0 #seconds 
    omega=1 #m/s^2 initial velocity 
    theta=math.radians(15) #rad initial angle
    delta_t=0.01 #second
    q=2*gamma #Damp factor
    if DEBUG==True:
        print(f"q: {q}")
        print(f"theta: {theta}")

    #Condition for Damped 
    if q<(2*math.sqrt(g/l)): #Determing type of dampining based on damp factor 
        type_damp=r"Underdamped ($\gamma < \omega_0$)"
    elif q==(2*math.sqrt(g/l)):
        type_damp=r"Critically Damped ($\gamma = \omega_0$)"
    elif q>(2*math.sqrt(g/l)):
        type_damp=r"Overdamped ($\gamma > \omega_0$)"
    if DEBUG==True:
        print(f"type_damp: {type_damp}")
    
    #Numerical caluclations - Undamped
    for i in range(10000): #proving a signficant amount of steps 
        omega_values.append(omega)
        theta_values.append(theta)
        t_values.append(t)
        omega=omega-(g/l)*theta*delta_t 
        theta=theta+omega*delta_t
        t=t+delta_t
        if abs(omega)<=1e-8 and abs(theta)<=1e-8: #exit loop when close to 0 
            break
    if DEBUG==True:
        print(f"omega: {omega_values}")
        print(f"angle: {theta_values}")
        print(f"time values: {t_values}")

    #Resetting needed variables 
    t=0 #seconds 
    omega=1 #m/s^2 initial velocity 
    theta=math.radians(15) #rad initial angle

    #Numerical calculation - Damped
    for i in range(10000): #Over a signficant amount of steps 
        d_omega_values.append(omega)
        d_theta_values.append(theta)
        d_t_values.append(t)
        omega=omega-((g/l)*theta*delta_t+q*omega*delta_t)
        theta=theta+omega*delta_t
        t=t+delta_t
        if abs(omega)<=1e-8 and abs(theta)<=1e-8:
            break
    if DEBUG==True:
        print(f"omega: {d_omega_values}")
        print(f"angle: {d_theta_values}")
        print(f"time values: {d_t_values}")

    #Plot 
    fig,(ax1,ax2)=plt.subplots(2,1, figsize=(12,6))
    ax1.plot(t_values,theta_values,color='blue',label='Undamped')
    ax1.plot(d_t_values, d_theta_values,color='green',label=type_damp)
    ax1.set_ylabel("Angular Discplacment (rad)")
    ax1.legend(loc="upper left", bbox_to_anchor=(1, 1))
    ax2.plot(t_values,omega_values,color='blue',label='Undamped')
    ax2.plot(d_t_values,d_omega_values,color='green',label=type_damp)
    ax2.set_ylabel(r"Angular Velocity ($\frac{\text{rad}}{\text{s}}$)")
    ax2.legend(loc="upper left", bbox_to_anchor=(1, 1))
    fig.supxlabel("Time (s)")
    fig.suptitle(f"Oscillitory Motion - Euler-Cromer")
    for ax in [ax1,ax2]:
        ax.axhline(y=0,color='black',linewidth=1,zorder=1,alpha=0.3,linestyle='--') 
        y_min,y_max = ax.get_ylim()
        limit=max(abs(y_min),abs(y_max))
        ax.set_ylim(-limit,limit)
        ax.set_xlim([0, max(d_t_values)])
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)
    plt.show()

#calling functions 
partone(DEBUG=False)