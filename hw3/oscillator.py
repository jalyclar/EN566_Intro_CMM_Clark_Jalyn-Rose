import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.signal import find_peaks 
from scipy.signal import peak_widths
from scipy.stats import linregress

################       Add error analysis on calculated gamma, compare to error given in notes     ##############################
################      Add optional arguments sys.argv for the different parts of problem    ##############################
################      Adjust self. varibales assinged in init    ##############################
################       Add more DEBUG statments   ##############################
################       Adjust text and spacial formatting on plots     ##############################

class oscillatory:
    def __init__(self):
        #Global variables 
        self.g=9.8 #m/s^2
        self.l=9.8 #m
        self.gamma=0.25 #s^-1
        self.omega_values=[]
        self.theta_values=[]
        self.t_values=[]
        self.d_omega_values=[]
        self.d_theta_values=[]
        self.d_t_values=[]
        self.rk_omega_values=[]
        self.rk_theta_values=[]
        self.rk_t_values=[]
        self.rk_drive_values=[]
        self.nl_omega=[]
        self.nl_theta=[]
        self.nl_t=[]
        self.nl_drive=[] 
        self.t_0=0 #seconds 
        self.omega_i=1 #m/s^2 initial velocity 
        self.theta_0=math.radians(15) #rad initial angle
        self.delta_t=0.01 #second
        self.q=2*self.gamma #Damp factor
        self.DEBUG=False 
    
    def linear_simulation(self,omega_D=1,t_array=np.arange(1,10,0.1),delta_t=1):
        self.omega_values=[]
        self.theta_values=[]
        self.t_values=[]
        self.d_omega_values=[]
        self.d_theta_values=[]
        self.d_t_values=[]
        self.rk_omega_values=[]
        self.rk_theta_values=[]
        self.rk_t_values=[]
        self.rk_drive_values=[]

        t=self.t_0 
        omega=self.omega_i
        theta=self.theta_0
        q=self.q
        if self.DEBUG==True:
            print(f"q: {q}\n\n")
            print(f"theta: {theta}\n\n") 

        #Condition for Damped 
        if q<(2*math.sqrt(self.g/self.l)): #Determing type of dampining based on damp factor 
            self.type_damp=r"Underdamped ($\gamma < \omega_0$)"
        elif q==(2*math.sqrt(self.g/self.l)):
            self.type_damp=r"Critically Damped ($\gamma = \omega_0$)"
        elif q>(2*math.sqrt(self.g/self.l)):
            self.type_damp=r"Overdamped ($\gamma > \omega_0$)"
        if self.DEBUG==True:
            print(f"type_damp: {self.type_damp}\n\n")
        
        #Numerical caluclations - Undamped - Euler-cromer method 
        for i in range(1, len(t_array)): #providing a signficant amount of steps 
            self.omega_values.append(omega)
            self.theta_values.append(theta)
            self.t_values.append(t)
            omega=omega-(self.g/self.l)*theta*delta_t 
            theta=theta+omega*delta_t
            t=t+delta_t
            #if abs(omega)<=1e-8 and abs(theta)<=1e-8: #exit loop when close to 0 
                #break
        if self.DEBUG==True:
            print(f"undamped omega: {self.omega_values[1:30]}\n\n")
            print(f"undamped angle: {self.theta_values[1:30]}\n\n")
            print(f"undamped time values: {self.t_values[1:30]}\n\n")

        #Resetting needed variables 
        t=self.t_0
        omega=self.omega_i 
        theta=self.theta_0 

        #Numerical calculation - Damped - Euler-cromer method 
        for i in range(1, len(t_array)): #Over a signficant amount of steps 
            self.d_omega_values.append(omega)
            self.d_theta_values.append(theta)
            self.d_t_values.append(t)
            omega=omega-((self.g/self.l)*theta*delta_t+q*omega*delta_t)
            theta=theta+omega*delta_t
            t=t+delta_t
            #if abs(omega)<=1e-8 and abs(theta)<=1e-8:
                #break
        if self.DEBUG==True:
            print(f"damped omega: {self.d_omega_values[1:30]}\n\n")
            print(f"damped angle: {self.d_theta_values[1:30]}\n\n")
            print(f"damped time values: {self.d_t_values[1:30]}\n\n")
        
        #Resetting needed variables 
        t=self.t_0
        omega=self.omega_i
        theta=self.theta_0 
        alpha=0.2 #rad/s^2

        #Numerical calculation - Driven - Runge-Kutta 4th Order
        #for i in range(10000): 
        for i in range(1, len(t_array)):
            k1_omega=delta_t*((-self.g/self.l)*theta-2*self.gamma*omega+alpha*math.sin(omega_D*t))
            k1_theta=delta_t*omega
            k2_omega=delta_t*((-self.g/self.l)*(theta+(k1_theta/2))-2*self.gamma*(omega+(k1_omega/2))+alpha*math.sin(omega_D*(t+(delta_t/2)))) 
            k2_theta=delta_t*(omega+(k1_omega/2))
            k3_omega=delta_t*((-self.g/self.l)*(theta+(k2_theta/2))-2*self.gamma*(omega+(k2_omega/2))+alpha*math.sin(omega_D*(t+(delta_t/2))))
            k3_theta=delta_t*(omega+(k2_omega/2))
            k4_omega=delta_t*((-self.g/self.l)*(theta+k3_theta)-2*self.gamma*(omega+k3_omega)+alpha*math.sin(omega_D*(t+delta_t))) 
            k4_theta=delta_t*(omega+k3_omega)
            omega_update=omega+(1/6)*(k1_omega+2*k2_omega+2*k3_omega+k4_omega)
            theta_update=theta+(1/6)*(k1_theta+2*k2_theta+2*k3_theta+k4_theta)
            self.rk_omega_values.append(omega_update)
            self.rk_theta_values.append(theta_update)
            self.rk_t_values.append(t)
            self.rk_drive_values.append(math.cos(omega_D * t))
            t=t+delta_t
            #if abs(omega_update-omega)<=1e-5 and abs(theta_update-theta)<=1e-5:
                #break
            omega=omega_update
            theta=theta_update
        if self.DEBUG==True:
            print(f"driven omega: {self.rk_omega_values[1:30]}]\n\n")
            print(f"driven angle: {self.rk_theta_values[1:30]}\n\n")
            print(f"driven time values: {self.rk_t_values[1:30]}\n\n")

    def nonlinear_simulation(self,omega_D=1,t_array=np.arange(1,10,0.1),delta_t=1):
        self.nl_omega=[]
        self.nl_theta=[]
        self.nl_t=[]
        self.nl_drive=[]  

        delta_t=self.delta_t
        alpha=0.2 #rad/s^2
        t_max=100
        t_array=np.arange(self.t_0,t_max,delta_t)
        theta=self.theta_0
        omega=self.omega_i
        t=self.t_0
        for i in range(1, len(t_array)):
            k1_omega=delta_t*((-self.g/self.l)*math.sin(theta)-2*self.gamma*omega+alpha*math.sin(omega_D*t))
            k1_theta=delta_t*omega
            k2_omega=delta_t*((-self.g/self.l)*(math.sin(theta+(k1_theta/2)))-2*self.gamma*(omega+(k1_omega/2))+alpha*math.sin(omega_D*(t+(delta_t/2)))) 
            k2_theta=delta_t*(omega+(k1_omega/2))
            k3_omega=delta_t*((-self.g/self.l)*(math.sin(theta+(k2_theta/2)))-2*self.gamma*(omega+(k2_omega/2))+alpha*math.sin(omega_D*(t+(delta_t/2))))
            k3_theta=delta_t*(omega+(k2_omega/2))
            k4_omega=delta_t*((-self.g/self.l)*(math.sin(theta+k3_theta))-2*self.gamma*(omega+k3_omega)+alpha*math.sin(omega_D*(t+delta_t))) 
            k4_theta=delta_t*(omega+k3_omega)
            omega_update=omega+(1/6)*(k1_omega+2*k2_omega+2*k3_omega+k4_omega)
            theta_update=theta+(1/6)*(k1_theta+2*k2_theta+2*k3_theta+k4_theta)
            self.nl_omega.append(omega_update)
            self.nl_theta.append(theta_update)
            self.nl_t.append(t)
            self.nl_drive.append(math.cos(omega_D*t))
            t=t+delta_t
            omega=omega_update
            theta=theta_update
    
    def partone(self):
        #Calculating pendulum simulation using single omega_D value
        omega_D_fix=0.935
        t_max=100
        t_array=np.arange(self.t_0,t_max,self.delta_t)
        self.linear_simulation(omega_D_fix,t_array,self.delta_t)
        
        #Plots - Pendulum Simulation 
        fig,(ax1,ax2)=plt.subplots(2,1, figsize=(10,6))
        ax1.plot(self.t_values,self.theta_values,color='blue',label='Undamped')
        ax1.plot(self.d_t_values, self.d_theta_values,color='green',label=self.type_damp)
        ax1.plot(self.rk_t_values,self.rk_theta_values,color='red',label='Driven')
        ax1.set_ylabel("Angular Discplacment (rad)")
        ax1.legend(loc="upper right")
        ax2.plot(self.t_values,self.omega_values,color='blue',label='Undamped')
        ax2.plot(self.d_t_values,self.d_omega_values,color='green',label=self.type_damp)
        ax2.plot(self.rk_t_values,self.rk_omega_values,color='red',label='Driven')
        ax2.set_ylabel(r"Angular Velocity ($\frac{\text{rad}}{\text{s}}$)")
        ax2.legend(loc="upper right")
        fig.supxlabel("Time (s)")
        fig.suptitle(r"Oscillitory Motion (Linear) Pendulum When $\Omega_D$=0.935 $\frac{rad}{sec}$")
        for ax in [ax1,ax2]:
            ax.axhline(y=0,color='black',linewidth=1,zorder=1,alpha=0.3,linestyle='--') 
        plt.tight_layout()
        plt.subplots_adjust(right=0.8)
        plt.show()

        #Different driving frequencies analysis 
        theta_values=[]
        phi_values=[]
        omega_natu=math.sqrt(self.g/self.l)
        omega_D_range=np.arange(omega_natu*0.3,omega_natu*2.5,0.01)
        if self.DEBUG==True:
            print(omega_D_range)

        #Calculating driven pendulum simulation for range of omega_D values 
        for omega_D in omega_D_range:
            self.rk_theta_values=[]#Ressting lists before each run and from before 
            self.rk_t_values=[]
            self.rk_drive_values=[]
            t_max=100
            t_array=np.arange(self.t_0,t_max,self.delta_t)
            self.linear_simulation(omega_D,t_array,self.delta_t)
            #self.t_values,self.theta_values, self.d_t_values, self.d_theta_values, self.type_damp, self.rk_t_values,self.rk_theta_values,self.omega_values,self.d_omega_values,self.rk_omega_values,self.rk_drive_values=
            #Extracting around steady-state 
            ss_theta=self.rk_theta_values[-len(self.rk_theta_values)//4:] #Last 1/4 of values (around steady portion)
            ss_time=self.rk_t_values[-len(self.rk_t_values)//4:]
            ss_drive=self.rk_drive_values[-len(self.rk_drive_values)//4:]
            #Finding amplitude as the max in the ssteady state 
            theta=max(np.abs(ss_theta))
            theta_values.append(theta)
            #Calculating phase shift comparing peak of theta and peak of driving force
            t_theta_max = ss_time[np.argmax(ss_theta)]
            t_drive_max = ss_time[np.argmax(ss_drive)]
            phi = omega_D * (t_theta_max - t_drive_max)
            phi = (phi + 2*np.pi) % (2*np.pi) #Normalizes phi to between 0-2pi
            phi_values.append(phi)
        if self.DEBUG==True:
            print(f"phi: {phi_values}\n\n")
            print(f"omega_D range: {omega_D_range}\n\n")
            print(f"theta_0 values: {theta_values}\n\n")

        #Plots - Different driving frequencies resonance structure mapping analysis  
        fig,(ax3,ax4)=plt.subplots(2,1, figsize=(10,6))
        ax3.plot(omega_D_range, theta_values)
        ax3.set_xlabel(r'Driving Frequency $\Omega_\text{D} (\frac{\text{rad}}{\text{sec}})$')
        ax3.set_ylabel(r'Amplitude $\theta_0(\Omega_\text{D})$(rad)')
        ax4.plot(omega_D_range, phi_values)
        ax4.set_xlabel(r'Driving Frequency $\Omega_\text{D} (\frac{\text{rad}}{\text{sec}})$')
        ax4.set_ylabel(r'Phase Shift $\phi(\Omega_\text{D})$(rad)')
        ax.set_ylim([0, max(theta_values)])
        fig.suptitle('Resonance Structure Mapping (Linear Pnedulum)')
        plt.tight_layout()
        plt.show()

        #FWHM of resonnace curve calculation
        res_values=np.array(theta_values) #Peaks function only takes 1D arrays
        omeg_D_values=np.array(omega_D_range)
        res_peaks,_=find_peaks(res_values,distance=2) #Calculating peak of resonance curve 
        final_peak=res_values[res_peaks] #Assinging indicie of peak to value from array 
        result=peak_widths(res_values,res_peaks,rel_height=0.5) #Calculates FWHM
        fwhm=result[0]*omeg_D_values[1]-omeg_D_values[0] #Assings just the width 
        print(f"fwhm of resoannce: {fwhm}")
        plt.figure()
        plt.plot(res_values)
        plt.plot(res_peaks,res_values[res_peaks],'x')
        plt.hlines(*result[1:],color='C2')
        plt.show()
        ################       Add error analysis between fwhm and gamma    ##############################
        if self.DEBUG==True:
            print(f"resoannce values type: {type(res_values)}\n\n")
            print(f"peak of resonance: {final_peak}\n\n")
            plt.figure()
            plt.plot(res_values)
            plt.plot(res_peaks,res_values[res_peaks],'x')
            plt.hlines(*result[1:],color='C2')
            plt.show()
        
    def parttwo(self): #Simulating near resoannce and calculating energies
        ke_values=[]
        pe_values=[]
        te_values=[]
        periods=[]
        period=(2*math.pi)/0.935
        t_array=np.arange(self.t_0,10*period,self.delta_t)
        self.rk_omega_values=[]
        self.rk_theta_values=[]#Ressting lists before each run and from before 
        self.rk_t_values=[]
        self.rk_drive_values=[]
        self.linear_simulation(0.935,t_array,self.delta_t) #At a driving frequency close to resonance 
        m=0.5 #kg 
        for omega in self.rk_omega_values:
        #KE calculation
            ke=(1/2)*m*(self.l**2)*omega
            ke_values.append(ke)
        for theta in self.rk_theta_values:
            #PE calculation 
            pe=m*self.g*self.l*(1-math.cos(theta))
            pe_values.append(pe)
        #TE calculation
        for ke,pe in zip(ke_values,pe_values): 
            te=pe+ke
            te_values.append(te)
        if self.DEBUG==True:
            print(f"Length of ke_values: {len(ke_values)}\n\n")
            print(f"Length of pe_values: {len(pe_values)}\n\n")
            print(f"Length of te_values: {len(te_values)}\n\n")
            print(f"Length of time values: {len(self.rk_t_values)}\n\n")
        for t in self.rk_t_values:
            p=t/period
            periods.append(p)
        #Plotting energies 
        plt.figure()
        plt.plot(periods,ke_values,label='Kinetic Energy')
        plt.plot(periods,pe_values,label='Potential Energy')
        plt.plot(periods,te_values,label='Total Energy')
        plt.xlabel("Periods (s)")
        plt.ylabel('Energy (J)')
        plt.title('Energy at Driving Frequency Near Resonnace')
        plt.legend()
        plt.show()

    def partthree(self): #Non-linear effects
        omega_D=0.935 #Omega_D clost to resonance
        t_max=100
        t_array=np.arange(self.t_0,t_max,self.delta_t)
        self.nonlinear_simulation(omega_D,t_array,self.delta_t)
        self.linear_simulation(omega_D,t_array,self.delta_t)

        fig,(ax1,ax2)=plt.subplots(2,1, figsize=(10,6))
        ax1.plot(self.nl_t,self.nl_theta,label='Non-Linear')
        ax1.plot(self.rk_t_values,self.rk_theta_values,color='red',label='Linear')
        ax1.set_ylabel("Angular Discplacment (rad)")
        ax1.legend(loc="best")
        ax2.plot(self.nl_t,self.nl_omega,label='Non-Linear')
        ax2.plot(self.rk_t_values,self.rk_omega_values,color='red',label='Linear')
        ax2.set_ylabel(r"Angular Velocity ($\frac{\text{rad}}{\text{s}}$)")
        ax2.legend(loc="best")
        fig.supxlabel('Time (s)')
        fig.suptitle(r"Oscillitory Motion of Pendulum When $\Omega_D$=0.935 $\frac{rad}{sec}$")
        plt.show()

        omega_D=1.2 #Omega_D higher than resonance 
        self.nonlinear_simulation(omega_D,t_array,self.delta_t)
        self.linear_simulation(omega_D,t_array,self.delta_t)
        
        fig,(ax1,ax2)=plt.subplots(2,1, figsize=(10,6))
        ax1.plot(self.nl_t,self.nl_theta,label='Non-Linear')
        ax1.plot(self.rk_t_values,self.rk_theta_values,color='red',label='Linear')
        ax1.set_ylabel("Angular Discplacment (rad)")
        ax1.legend(loc="best")
        ax2.plot(self.nl_t,self.nl_omega,label='Non-Linear')
        ax2.plot(self.rk_t_values,self.rk_omega_values,color='red',label='Linear')
        ax2.set_ylabel(r"Angular Velocity ($\frac{\text{rad}}{\text{s}}$)")
        ax2.legend(loc="best")
        fig.supxlabel('Time (s)')
        fig.suptitle(r"Oscillitory Motion Pendulum When $\Omega_D$=1.2 $\frac{rad}{sec}$")
        plt.show()

    def partfour(self):
        delta_t=self.delta_t
        alpha=[0.2,0.5,1.2]
        omega_D=0.666
        t_max=100
        t_array=np.arange(self.t_0,t_max,delta_t)
        theta=[15.000,15.001]
        stability={}
        for traj in alpha:
            delta_theta_values=[[],[]]
            time_values=[[],[]]
            for idx, angle in enumerate(theta):
                omega=self.omega_i
                t=self.t_0
                for i in range(1, len(t_array)):
                    k1_omega=delta_t*((-self.g/self.l)*math.sin(angle)-2*self.gamma*omega+traj*math.sin(omega_D*t))
                    k1_theta=delta_t*omega
                    k2_omega=delta_t*((-self.g/self.l)*(math.sin(angle+(k1_theta/2)))-2*self.gamma*(omega+(k1_omega/2))+traj*math.sin(omega_D*(t+(delta_t/2)))) 
                    k2_theta=delta_t*(omega+(k1_omega/2))
                    k3_omega=delta_t*((-self.g/self.l)*(math.sin(angle+(k2_theta/2)))-2*self.gamma*(omega+(k2_omega/2))+traj*math.sin(omega_D*(t+(delta_t/2))))
                    k3_theta=delta_t*(omega+(k2_omega/2))
                    k4_omega=delta_t*((-self.g/self.l)*(math.sin(angle+k3_theta))-2*self.gamma*(omega+k3_omega)+traj*math.sin(omega_D*(t+delta_t))) 
                    k4_theta=delta_t*(omega+k3_omega)
                    omega_update=omega+(1/6)*(k1_omega+2*k2_omega+2*k3_omega+k4_omega)
                    theta_update=angle+(1/6)*(k1_theta+2*k2_theta+2*k3_theta+k4_theta)
                    theta_update=(theta_update+np.pi)%(2*np.pi)-np.pi
                    time_values[idx].append(t)
                    delta_theta_values[idx].append(theta_update)
                    t=t+delta_t
                    omega=omega_update
                    angle=theta_update
            differences_theta=[]
            for i in range(len(delta_theta_values[0])):
                difference_theta=np.abs(delta_theta_values[0][i]-delta_theta_values[1][i])
                differences_theta.append(difference_theta)
            stability[traj]=[time_values[0],differences_theta]
        plt.figure()
        for traj,values in stability.items():
            time_values,differences_theta=values
            y=np.log(np.array(differences_theta)+1e-10)
            plt.plot(time_values,y,label=f'{traj}')
        plt.xlabel('Time (s)')
        plt.ylabel(r'log ($ \Delta \theta $)')
        plt.title('Pendulum Stability for Lyapunov Estimation')
        plt.legend()
        plt.show()

        start_i=1000
        end_i=8000
        lyap_values={}
        plt.figure()
        for traj,values in stability.items():
            time_values,differences_theta=values
            log_differences_theta=np.log(np.array(differences_theta)+1e-10)
            lin_reg_x=np.array(time_values[start_i:end_i])
            lin_reg_y=log_differences_theta[start_i:end_i]
            if self.DEBUG==True:
                print(f"lin reg x values: {lin_reg_x}\n\n")
                print(f"lin reg y values: {lin_reg_y}\n\n")
            slope,intercept,r,p,e=linregress(lin_reg_x,lin_reg_y)
            lyap_values[traj]=slope
            plt.plot(time_values,log_differences_theta,label=f'{traj}')
            plt.plot(lin_reg_x,slope*lin_reg_x+intercept,'r--',label=fr'Fit for {traj}: $\lambda$={slope:.2f}')
        plt.xlabel('Time (s)')
        plt.ylabel(r'log ($ \Delta \theta $)')
        plt.title('Lyapunov Exponent Estimation')
        plt.legend()
        plt.show()
        print(lyap_values)

if __name__=="__main__":
    data=oscillatory()
    data.nonlinear_simulation()
    data.linear_simulation()
    data.partone()
    data.parttwo()
    data.partthree()
    data.partfour()
    