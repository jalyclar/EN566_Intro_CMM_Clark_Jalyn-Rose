import numpy as np
import matplotlib.pyplot as plt
import math

def undamped(DEBUG=False): #Only using euler-cromer method 

    #Variables 
    g=9.8 #m/s^2
    l=9.8 #m
    omega_values=[]
    theta_values=[]
    t_values=[]
    delta_t=0.01 #second
    t=0 #seconds 
    omega=1 #m/s^2 initial velocity 
    theta=math.radians(15) #rad initial angle 

    #Numerical caluclations 
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

    #Plot
    fig,(ax1,ax2)=plt.subplots(2,1, figsize=(8,6))
    ax1.plot(t_values,theta_values,color='blue')
    ax1.set_ylabel("Angular Discplacment")
    ax2.plot(t_values,omega_values,color='pink')
    ax2.set_ylabel("Angular velocity")
    fig.supxlabel("Time (s)")
    fig.suptitle("Undamped Solution")
    for ax in [ax1,ax2]:
        ax.axhline(y=0, color='black', linewidth=1, zorder=1) 
        y_min,y_max = ax.get_ylim()
        limit=max(abs(y_min),abs(y_max))
        ax.set_ylim(-limit,limit)
    plt.tight_layout()
    plt.show()

def damped(DEBUG=False): #Only using euler-cromer method 
   
   #Variables 
    g=9.8 #m/s^2
    l=9.8 #m
    gamma=0.25 #s^-1
    alpha_d=0.2 #rad/s^2 
    d_omega_values=[]
    d_theta_values=[]
    d_t_values=[]
    delta_t=0.01 #second
    t=0 #seconds 
    omega=1 #m/s^2 initial velocity 
    theta=math.radians(15) #rad initial angle
    q=2*gamma #Damp factor

    #Condition
    if q<(2*math.sqrt(g/l)): #Determing type of dampining based on damp factor 
        type_damp="Underdamped"
    elif q==(2*math.sqrt(g/l)):
        type_damp="Critically Damped"
    elif q>(2*math.sqrt(g/l)):
        type_damp="Overdamped" 

    #Numerical calculation 
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
    fig,(ax1,ax2)=plt.subplots(2,1, figsize=(8,6))
    ax1.plot(d_t_values, d_theta_values,color='pink')
    ax1.set_ylabel("Angular Discplacment")
    ax2.plot(d_t_values,d_omega_values,color='blue')
    ax2.set_ylabel("Angular velocity")
    fig.supxlabel("Time (s)")
    fig.suptitle(f"{type_damp} Solution")
    for ax in [ax1,ax2]:
        ax.axhline(y=0, color='black', linewidth=1, zorder=1) 
        y_min,y_max = ax.get_ylim()
        limit=max(abs(y_min),abs(y_max))
        ax.set_ylim(-limit,limit)
    plt.tight_layout()
    plt.show()

#calling functions 
undamped(DEBUG=False)
damped(DEBUG=False)