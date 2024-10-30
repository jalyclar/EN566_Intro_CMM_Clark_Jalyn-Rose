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
            plt.title(f'Random Numbers - Bin Size={bin_size}')
            plt.show()

def parttwo():
    sigma=1.0
    N=1000
    p_values=[]
    for i in range(N):
        #Use Rejection method  
        #Use LGC to generate linear combination of random numbers 
        rnum=random.random()
        p=(1/sigma*np.sqrt(2*np.pi))**(-(x**2)/2*sigma**2)
        p_values.append(p)
    



if __name__=="__main__":
    partone()