import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.signal import find_peaks 
from scipy.signal import peak_widths

def partone(DEBUG=False):
    def pendulum_simulation(omega_D):
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
        rk_omega_values=[]
        rk_theta_values=[]
        rk_t_values=[]
        rk_drive_values=[]
        t=0 #seconds 
        omega=1 #m/s^2 initial velocity 
        theta=math.radians(15) #rad initial angle
        delta_t=0.001 #second
        q=2*gamma #Damp factor
        if DEBUG==True:
            print(f"q: {q}\n\n")
            print(f"theta: {theta}\n\n") 
        t_max=100
        t_array=np.arange(t,t_max,delta_t)

        #Condition for Damped 
        if q<(2*math.sqrt(g/l)): #Determing type of dampining based on damp factor 
            type_damp=r"Underdamped ($\gamma < \omega_0$)"
        elif q==(2*math.sqrt(g/l)):
            type_damp=r"Critically Damped ($\gamma = \omega_0$)"
        elif q>(2*math.sqrt(g/l)):
            type_damp=r"Overdamped ($\gamma > \omega_0$)"
        if DEBUG==True:
            print(f"type_damp: {type_damp}\n\n")
        
        #Numerical caluclations - Undamped - Euler-cromer method 
        for i in range(1, len(t_array)): #proving a signficant amount of steps 
            omega_values.append(omega)
            theta_values.append(theta)
            t_values.append(t)
            omega=omega-(g/l)*theta*delta_t 
            theta=theta+omega*delta_t
            t=t+delta_t
            #if abs(omega)<=1e-8 and abs(theta)<=1e-8: #exit loop when close to 0 
                #break
        if DEBUG==True:
            print(f"undamped omega: {omega_values[1:30]}\n\n")
            print(f"undamped angle: {theta_values[1:30]}\n\n")
            print(f"undamped time values: {t_values[1:30]}\n\n")

        #Resetting needed variables 
        t=0 #seconds 
        omega=1 #m/s^2 initial velocity 
        theta=math.radians(15) #rad initial angle

        #Numerical calculation - Damped - Euler-cromer method 
        for i in range(1, len(t_array)): #Over a signficant amount of steps 
            d_omega_values.append(omega)
            d_theta_values.append(theta)
            d_t_values.append(t)
            omega=omega-((g/l)*theta*delta_t+q*omega*delta_t)
            theta=theta+omega*delta_t
            t=t+delta_t
            #if abs(omega)<=1e-8 and abs(theta)<=1e-8:
                #break
        if DEBUG==True:
            print(f"damped omega: {d_omega_values[1:30]}\n\n")
            print(f"damped angle: {d_theta_values[1:30]}\n\n")
            print(f"damped time values: {d_t_values[1:30]}\n\n")
        
        #Resetting needed variables 
        t=0 #seconds 
        omega=1 #m/s^2 initial velocity 
        theta=math.radians(15) #rad initial angle
        alpha=0.2 #rad/s^2
        #omega_D=0.935 #rad/s

        #Numerical calculation - Driven - Runge-Kutta 4th Order
        #for i in range(10000): 
        for i in range(1, len(t_array)):
            k1_omega=delta_t*((-g/l)*theta-2*gamma*omega+alpha*math.sin(omega_D*t)) 
            k1_theta=delta_t*omega
            k2_omega=delta_t*((-g/l)*(theta+(k1_theta/2))-2*gamma*(omega+(k1_omega/2))+alpha*math.sin(omega_D*(t+(delta_t/2)))) 
            k2_theta=delta_t*(omega+(k1_omega/2))
            k3_omega=delta_t*((-g/l)*(theta+(k2_theta/2))-2*gamma*(omega+(k2_omega/2))+alpha*math.sin(omega_D*(t+(delta_t/2))))
            k3_theta=delta_t*(omega+(k2_omega/2))
            k4_omega=delta_t*((-g/l)*(theta+k3_theta)-2*gamma*(omega+k3_omega)+alpha*math.sin(omega_D*(t+delta_t))) 
            k4_theta=delta_t*(omega+k3_omega)
            omega_update=omega+(1/6)*(k1_omega+2*k2_omega+2*k3_omega+k4_omega)
            theta_update=theta+(1/6)*(k1_theta+2*k2_theta+2*k3_theta+k4_theta)
            rk_omega_values.append(omega_update)
            rk_theta_values.append(theta_update)
            rk_t_values.append(t)
            rk_drive_values.append(math.cos(omega_D * t))
            t=t+delta_t
            #if abs(omega_update-omega)<=1e-5 and abs(theta_update-theta)<=1e-5:
                #break
            omega=omega_update
            theta=theta_update
        if DEBUG==True:
            print(f"driven omega: {rk_omega_values[1:30]}]\n\n")
            print(f"driven angle: {rk_theta_values[1:30]}\n\n")
            print(f"driven time values: {rk_t_values[1:30]}\n\n")
        return t_values,theta_values, d_t_values, d_theta_values, type_damp, rk_t_values,rk_theta_values,omega_values,d_omega_values,rk_omega_values,rk_drive_values

    #Calculating pendulum simulation using single omega_D value
    omega_D_fix=0.935
    t_values,theta_values, d_t_values, d_theta_values, type_damp, rk_t_values,rk_theta_values,omega_values,d_omega_values,rk_omega_values,rk_drive_values=pendulum_simulation(omega_D_fix)

    #Plots - Pendulum Simulation 
    fig,(ax1,ax2)=plt.subplots(2,1, figsize=(12,6))
    ax1.plot(t_values,theta_values,color='blue',label='Undamped')
    ax1.plot(d_t_values, d_theta_values,color='green',label=type_damp)
    ax1.plot(rk_t_values,rk_theta_values,color='red',label='Driven')
    ax1.set_ylabel("Angular Discplacment (rad)")
    ax1.legend(loc="upper left", bbox_to_anchor=(1, 1))
    ax2.plot(t_values,omega_values,color='blue',label='Undamped')
    ax2.plot(d_t_values,d_omega_values,color='green',label=type_damp)
    ax2.plot(rk_t_values,rk_omega_values,color='red',label='Driven')
    ax2.set_ylabel(r"Angular Velocity ($\frac{\text{rad}}{\text{s}}$)")
    ax2.legend(loc="upper left", bbox_to_anchor=(1, 1))
    fig.supxlabel("Time (s)")
    fig.suptitle(r"Oscillitory Motion - $\Omega_D$=0.935 $\frac{rad}{sec}$")
    for ax in [ax1,ax2]:
        ax.axhline(y=0,color='black',linewidth=1,zorder=1,alpha=0.3,linestyle='--') 
        ax.set_xlim([0, max(d_t_values)])
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)
    plt.show()

    #Different driving frequencies analysis 
    theta_0_values=[]
    phi_values=[]
    g=9.8 #m/s^2
    l=9.8 #m
    omega_0=math.sqrt(g/l)
    omega_D_range=np.arange(omega_0*0.1,omega_0*3,0.1)
    if DEBUG==True:
        print(omega_D_range)

    #Calculating driven pendulum simulation for range of omega_D values 
    for omega_D in omega_D_range:
        t_values,theta_values, d_t_values, d_theta_values, type_damp, rk_t_values,rk_theta_values,omega_values,d_omega_values,rk_omega_values,rk_drive_values=pendulum_simulation(omega_D)
        #Extracting around steady-state 
        ss_theta=rk_theta_values[-len(rk_theta_values)//4:] #Last 1/4 of values (around steady portion)
        ss_time=rk_t_values[-len(rk_t_values)//4:]
        ss_drive=rk_drive_values[-len(rk_drive_values)//4:]
        #Finding amplitude as the max in the ssteady state 
        theta_0=max(np.abs(ss_theta))
        theta_0_values.append(theta_0)
        #Calculating phase shift comparing peak of theta and peak of driving force
        t_theta_max = ss_time[np.argmax(ss_theta)]
        t_drive_max = ss_time[np.argmax(ss_drive)]
        phi = omega_D * (t_theta_max - t_drive_max)
        phi = (phi + 2*np.pi) % (2*np.pi) #Normalizes phi to between 0-2pi
        phi_values.append(phi)
    if DEBUG==True:
        print(f"phi: {phi_values}\n\n")
        print(f"omega_D range: {omega_D_range}\n\n")
        print(f"theta_0 values: {theta_0_values}\n\n")

    #Plots - Different driving frequencies resonance structure mapping analysis  
    fig,(ax3,ax4)=plt.subplots(2,1, figsize=(12,6))
    ax3.plot(omega_D_range, theta_0_values)
    ax3.set_xlabel(r'Driving Frequency $\Omega_\text{D} (\frac{\text{rad}}{\text{sec}})$')
    ax3.set_ylabel(r'Amplitude $\theta_0(\Omega_\text{D})$(rad)')
    ax4.plot(omega_D_range, phi_values)
    ax4.set_xlabel(r'Driving Frequency $\Omega_\text{D} (\frac{\text{rad}}{\text{sec}})$')
    ax4.set_ylabel(r'Phase Shift $\phi(\Omega_\text{D})$(rad)')
    ax.set_ylim([0, max(theta_0_values)])
    fig.suptitle('Resonance Structure Mapping')
    plt.tight_layout()
    plt.show()

    #FWHM of resonnace curve calculation
    res_values=np.array(theta_0_values)
    res_peaks,_=find_peaks(res_values,distance=2)
    final_peak=res_values[res_peaks]
    #plt.plot(final_peak)
    #plt.plot(omega_D_range, theta_0_values)
    plt.show()
    if DEBUG==True:
        print(res_values)
        print(final_peak)


#calling functions 
partone(DEBUG=False)