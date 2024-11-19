import numpy as np 
import random 
import matplotlib.pyplot as plt 

def partone():
    #Variables 
    n=50
    N=n*n
    j=1.5 #nn interaction stregnth 
    h=0 #no external magnetic field 
    kb=1.0 #boltzman constant 
    iterations=2000
    eqt=2000

    #Initialization 
    lats=np.empty((n,n),dtype=int)
    for i in range(n):
        for j in range(n):
            lats[i,j]=1 if random.random() <0.5 else -1
    temps=np.arange(1,5,0.01)
    mag=[] #to hold magnetism 

    #equlibrium first 
    for t in temps:
        for iter in range(iterations):
            for space in range(N):
                i=random.randint(0,n-1) #getting random sites 
                j=random.randint(0,n-1)
                spin=lats[i,j]
                nbs=(lats[(i-1)%n,j]+lats[(i+1)%n,j]+lats[i,(j-1)%n]+lats[i,(j+1)%n])
                deltaE=2*j*spin*nbs
                if deltaE<=0 or random.random()<np.exp(-deltaE/(kb*t)):
                    lats[i,j]*=-1
    #magnetism calculations
        m=0
        for it in range(eqt):
            for s in range(N):
                i=random.randint(0,n-1)
                j=random.randint(0,n-1)
                spin=lats[i,j]
                nbs=(lats[(i-1)%n,j]+lats[(i+1)%n,j]+lats[i,(j-1)%n]+lats[i,(j+1)%n])
                deltaE=2*j*spin*nbs
                if deltaE<=0 or random.random()<np.exp(-deltaE/(kb*t)):
                    lats[i,j]*=-1
            m+=np.abs(np.sum(lats))
        mag.append(m/(eqt*N))

    plt.figure()
    plt.plot(temps,mag,label='Magnetization')
    plt.axvline(x=2.27, color='red', linestyle='--', label='Critical Temperature')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Magnetization')
    plt.title('Magnetization as a Function fo Temperature - 2D Ising Model')
    plt.legend()
    plt.savefig('part1plot.png',dpi=300)

def parttwo():
    #Variables 
    lats=[5,10,20,30,40,50,75,100,200,500]
    j=1.5 #nn interaction stregnth 
    h=0 #no external magnetic field 
    kb=1.0 #boltzman constant 
    iterations=2000
    eqt=2000
    temps=np.arange(1,5,0.01)

    #Initilization 
    cmax_n=[]
    logn=[]
    nval=[]

    for n in lats:
        N=n*n 
        lats=np.empty((n,n),dtype=int)
        for i in range(n):
            for j in range(n):
                lats[i,j]=1 if random.random() <0.5 else -1
        c_n=[]

        for t in temps:
            #Eq
            for iter in range(iterations):
                for space in range(N):
                    i=random.randint(0,n-1) #getting random sites 
                    j=random.randint(0,n-1)
                    spin=lats[i,j]
                    nbs=(lats[(i-1)%n,j]+lats[(i+1)%n,j]+lats[i,(j-1)%n]+lats[i,(j+1)%n])
                    deltaE=2*j*spin*nbs
                    if deltaE<=0 or random.random()<np.exp(-deltaE/(kb*t)):
                        lats[i,j]*=-1
            #En
            energy=[]
            for it in range(eqt):
                for s in range(N):
                    i=random.randint(0,n-1)
                    j=random.randint(0,n-1)
                    spin=lats[i,j]
                    nbs=(lats[(i-1)%n,j]+lats[(i+1)%n,j]+lats[i,(j-1)%n]+lats[i,(j+1)%n])
                    deltaE=2*j*spin*nbs
                    if deltaE<=0 or random.random()<np.exp(-deltaE/(kb*t)):
                        lats[i,j]*=-1
                E=-j*np.sum(lats*(np.roll(lats,1,axis=0)+np.roll(lats,-1,axis=0)+np.roll(lats,1,axis=1)+np.roll(lats,-1,axis=1)))/2
                energy.append(E)

            av_E=np.mean(energy)
            av_E2=np.mean(np.square(energy))
            dE2=av_E2-av_E**2
            c=dE2/(kb*t**2)
            c_n.append(c/N)

        logn.append(np.log(n))
        cmax_n.append(max(c_n))
        nval.append(n)

        plt.figure()
        plt.plot(temps,c_n,label=f'Lattice size {n}')
        plt.xlabel('Temperature')
        plt.ylabel('Specific Heat per Spin (C/N)')
        plt.title(f'Specific Heat for Spin vs. Temperature - Lattice Size={n}')
        plt.legend()
        plt.savefig(f'C_T_lattice_{n}.png',dpi=300)

    plt.figure()
    plt.plot(nval,cmax_n,label='C_max/N')
    plt.xlabel('Lattice Size (n)')
    plt.ylabel('C_max/N')
    plt.title('C_max/N vs. Lattice Size (n)')
    plt.legend()
    plt.savefig('cmaxn_vs_n.png',dpi=300)

    plt.figure()
    plt.plot(logn,cmax_n,label='C_max/N')
    plt.xlabel('log(n)')
    plt.ylabel('C_max/N')
    plt.title('Finite-Size Scaling Relation')
    plt.legend()
    plt.savefig('cmaxn_vs_logn.png',dpi=300)

partone()
parttwo()




