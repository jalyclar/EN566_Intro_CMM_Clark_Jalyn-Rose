#rnumbers 

import random 
import matplotlib.pyplot as plt 
import numpy as np

def partone():
    N_values=[1000,1000000]
    bin_values=[10,20,50,100]

    for N in N_values:
        rnum_values=[]
        for i in range(N):
            rnum=random.random()
            rnum_values.append(rnum)
        
        for bin_size in bin_values:
            plt.hist(rnum_values,bin_size)
            plt.xlabel('Random Number')
            plt.ylabel('Frequency')
            plt.title(f'Random Numbers - Bin Size={bin_size} - Random Numbers={N}')
    plt.show()

def parttwo():
    #Box - Mueller Method
    sigma=1.0
    N=1000
    seeds=np.linspace(10,30,10)
    num_plots=len(seeds)
    for j, seed in enumerate(seeds):
        seed=int(seed)
        random.seed(seed)
        p_values=[]
        for i in range(N//2):
            u1=random.uniform(0,1) #2 randomly distrubted sets between 0 and 1 
            u2=random.uniform(0,1)
            X=np.sqrt(-2*np.log(u2))*np.cos(2*np.pi*u1)*sigma #Gauss distubution 
            Y=np.sqrt(-2*np.log(u2))*np.sin(2*np.pi*u1)*sigma
            p_values.extend([X,Y])
        x=np.linspace(-3*sigma,3*sigma,100)  # Full Gaussian
        gauss=(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(x**2)/(2*sigma**2))
        plt.figure()
        plt.hist(p_values,30,alpha=0.7,density=True)
        plt.plot(x,gauss)
        plt.xlabel('P(x)')
        plt.ylabel('Frequency')
        plt.title(f'Gaussian Distributed Random Numbers - Run {j}')
    plt.show()
    
if __name__=="__main__":
    partone()
    parttwo()