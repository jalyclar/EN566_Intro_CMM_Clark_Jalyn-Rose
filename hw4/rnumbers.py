#rnumbers 

import random 
import matplotlib.pyplot as plt 
import numpy as np

def partone(N=1000,bins=10):
    rnum_values=[]

    for i in range(N):
        rnum=random.random()
        rnum_values.append(rnum)

    plt.hist(rnum_values,bins)
    plt.xlabel('Random Number')
    plt.ylabel('Frequency')
    plt.show()

if __name__=="__main__":
    partone()