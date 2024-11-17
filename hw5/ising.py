import numpy as np 
import random 
import matplotlib.pyplot as plt 

#Variables 
n=50
N=n*n
j=1.5 #nn interaction stregnth 
h=0 #no external magnetic field 
kb=1.0 #boltzman constant 
iterations=2000

#Initialization 
lats=np.empty((n,n),dtype=int)
for i in range(n):
    for j in range(n):
        lats[i,j]=1 if random.random() <0.5 else -1
temps=np.arange(1,5,0.05)
mag=[] #to hold magnetism 

for t in temps:
    m=0
    for iter in range(iterations):
        for space in range(N):
            i=random.randint(0,n-1) #getting random sites 
            j=random.randint(0,n-1)
            spin=lats[i,j]
            nbs=(lats[(i-1)%n,j]+lats[(i+1)%n,j]+lats[i,(j-1)%n]+lats[i,(j+1)%n])
            deltaE=2*j*spin*nbs
            if deltaE<=0 or random.random()<np.exp(-deltaE/(kb*t)):
                lats[i,j]*=-1
        m+=np.abs(np.sum(lats))
    mag.append(m/(iterations*N))

plt.figure()
plt.plot(temps,mag)
plt.show()
    







