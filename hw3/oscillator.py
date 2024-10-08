import numpy as np
import matplotlib.pyplot as plt
import math

def linear(DEBUG=False):
    g=9.8 #m/s^2
    l=9.8 #m

    #Linear motion nondamped 
    omega_values=[]
    theta_values=[]
    t_values=[]

    delta_t=0.01 #second
    t=0 #seconds 
    omega=1 #m/s^2 initial velocity 
    theta=math.radians(15) #rad initial angle 
    for i in range(1000):
        omega_values.append(omega)
        theta_values.append(theta)
        t_values.append(t)
        omega=omega-(g/l)*theta*delta_t
        theta=theta+omega*delta_t
        t=t+delta_t
        if abs(omega) < 1e-5 and abs(theta) < 1e-5:
            break
    if DEBUG==True:
        print(f"omega: {omega_values}")
        print(f"angle: {theta_values}")
        print(f"time values: {t_values}")

    plt.figure()
    plt.plot(t_values, theta_values,label="Theta(t)")
    plt.plot(t_values,omega_values,label="Omega(t)")
    plt.xlabel("Time (s)")
    plt.ylabel("Angular velocity and Angular Discplacment")
    plt.title("Angular velocity and Angular Displacment vs Time")
    plt.legend()
    plt.show()

def damped(DEBUG=False):
    g=9.8 #m/s^2
    l=9.8 #m
    gamma=0.25 #s^-1
    alpha_d=0.2 #rad/s^2

    #Linear motion nondamped 
    d_omega_values=[]
    d_theta_values=[]
    d_t_values=[]

    delta_t=0.01 #second
    t=0 #seconds 
    omega=1 #m/s^2 initial velocity 
    theta=math.radians(15) #rad initial angle
    q=2*gamma #Damp factor

    for i in range(10000):
        if q<(2*math.sqrt(g/l)): #underdamped
            d_omega_values.append(omega)
            d_theta_values.append(theta)
            d_t_values.append(t)

            omega=omega-((g/l)*theta*delta_t+q*omega*delta_t)
            theta=theta+omega*delta_t
            t=t+delta_t
            if abs(omega) < 1e-5 and abs(theta) < 1e-5:
                break
    if DEBUG==True:
        print(f"omega: {d_omega_values}")
        print(f"angle: {d_theta_values}")
        print(f"time values: {d_t_values}")

        #elif q=3*omega: #critically damped
        #elif q>2*omega: #overdamped 
    
    plt.figure()
    plt.plot(d_t_values, d_theta_values,label="Theta(t)")
    plt.plot(d_t_values,d_omega_values,label="Omega(t)")
    plt.xlabel("Time (s)")
    plt.ylabel("Angular velocity and Angular Discplacment")
    plt.title("Angular velocity and Angular Displacment vs Time - Damped")
    plt.legend()
    plt.show()
    


#calling functions 
#linear(DEBUG=False)
damped(DEBUG=False)