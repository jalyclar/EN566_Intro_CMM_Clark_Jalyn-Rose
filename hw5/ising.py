#!/opt/anaconda3/bin/python 

import numpy as np 
import random 
import matplotlib.pyplot as plt 
import sys 

class ising_simulation():

    def __init__(self):
        self.J=1.5 #nn interaction stregnth 
        self.kb=1.0 #boltzman constant 
        self.iterations=5000 #iterations to reach equlibrium 
        self.iterations2=5000 #iterations for magnetism calulation
        self.temps=np.linspace(1,6,100) #Detailed temperature sweep
        self.parts=None

    def parse_arg(self):
        for arg in sys.argv[1:]:
            if "--part=" in arg:
                self.parts=arg.split('=')[1].split(',') 

    def get_lattice(self,n): #Generating the lattice with spin up or spin down 
        lattice=np.empty((n,n),dtype=int)
        for i in range(n):
            for j in range(n):
                lattice[i,j]=1 if random.random() <0.5 else -1
        return lattice

    def montecarlo_sim(self,spins,a,N,J): #monte carlo simulation for equlibrium to be reached 
        for space in range(N):
            i=random.randint(0,spins.shape[0]-1) #getting random sites 
            j=random.randint(0,spins.shape[0]-1)
            spin=spins[i,j]
            lattice_config=(spins[(i-1)%spins.shape[0],j]+spins[(i+1)%spins.shape[0],j]+spins[i,(j-1)%spins.shape[0]]+spins[i,(j+1)%spins.shape[0]])
            deltaE=2*J*spin*lattice_config
            if deltaE<=0 or random.random()<np.exp(-deltaE*a):
                spins[i,j]*=-1

    def energy_calc(self,spins,J): #calculating energy configuration of lattice 
        n=spins.shape[0]
        E=0
        for i in range(n):
            for j in range(n):
                nbs=(spins[(i-1)%n,j]+spins[(i+1)%n,j]+spins[i,(j-1)%n]+spins[i,(j+1)%n])
                E-=spins[i,j]*nbs
        return J*E/2

    def partone(self): #metropolis algorithim to calculate magnetism of system for a given T
        if '1' in self.parts: 
            #Initilize Values 
            n=50 
            N=n*n #total number of spins 
            spins=self.get_lattice(n)
            mag=[] #to hold magnetism values 

            for t in self.temps:
                a=1/(self.kb*t) #to use with probaiblity calulation of spin in deltaE
                for iter in range(self.iterations):
                    self.montecarlo_sim(spins,a,N,self.J)  #allowing equlibrium to be reached 
                m=0
                for it in range(self.iterations2):
                    self.montecarlo_sim(spins,a,N,self.J) #calulating magnetism values 
                    m+=np.abs(np.sum(spins)) #calculating magnetism of system for T
                mag.append(m/(self.iterations2*N)) 

            plt.figure() #Plotting Mag vs T for the system 
            plt.plot(self.temps,mag,label='Magnetization')
            plt.xlabel('Temperature (K)')
            plt.ylabel('Magnetization')
            plt.title('Magnetization as a Function fo Temperature - 2D Ising Model')
            plt.legend()
            plt.savefig('part1plot.png',dpi=300)

    def parttwo(self): #metropolis algortithim for specific heat/spin and verify approx finite scaling relation
        if '2' in self.parts:
            #Initilize Values
            lattices=[5,10,20,30,40,50,75,100,200,500] #lattice sizes to loop through 
            cmax_n=[] #to store specific heat/spin 
            logn=[] #to hold log(n) total spins 

            for n in lattices:
                N=n*n
                spins=lattices(n)
                c_n=[] #to store specific heat 

                for t in self.temps:
                    a=1/(self.kb*t) #to use with probaiblity calulation of spin in deltaE
                    for iter in range(self.iterations):
                        self.montecarlo_sim(spins,a,N,self.J) #allowing equlibrium to be reached 
                    Ener=[]
                    for it in range(self.iterations2):
                        self.montecarlo_sim(spins,a,N,self.J) 
                        en=self.energy_calc(spins,self.J) #for energy calculation using lattice from monte carlo
                        Ener.append(en)
                    av_E=np.mean(Ener) #average energy of system 
                    av_E2=np.mean(np.square(Ener)) #average enrgy squared 
                    dE2=av_E2-av_E**2 #change in energy delta E
                    c=dE2/(self.kb*t**2) #specific heat 
                    c_n.append(c/N) 
                logn.append(np.log(n))
                cmax_n.append(max(c_n))

                plt.figure() #plot of specific heat vs T 
                plt.plot(self.temps,c_n,label=f'Lattice size {n}')
                plt.xlabel('Temperature')
                plt.ylabel('Specific Heat per Spin (C/N)')
                plt.title(f'Specific Heat for Spin vs. Temperature - Lattice Size={n}')
                plt.legend()
                plt.savefig(f'C_T_lattice_{n}.png',dpi=300)

            plt.figure() #plot of lattice size vs max specific heat from temperatures looped over 
            plt.plot(lattices,cmax_n,label='C_max/N')
            plt.xlabel('Lattice Size (n)')
            plt.ylabel('C_max/N')
            plt.title('C_max/N vs. Lattice Size (n)')
            plt.legend()
            plt.savefig('cmaxn_vs_n.png',dpi=300)

            plt.figure() #plot of log(n) vs max specific heat per spin from temperatures looped over 
            plt.plot(logn,cmax_n,label='C_max/N')
            plt.xlabel('log(n)')
            plt.ylabel('C_max/N')
            plt.title('Finite-Size Scaling Relation')
            plt.legend()
            plt.savefig('cmaxn_vs_logn.png',dpi=300)

if __name__=="__main__":
    data=ising_simulation()
    data.parse_arg()
    data.partone()
    data.parttwo()


