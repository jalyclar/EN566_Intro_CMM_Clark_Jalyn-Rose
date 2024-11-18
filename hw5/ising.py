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
    with open('magdata.txt','w') as file:
        for item in mag:
            file.write(item+'\n')

plt.figure()
plt.plot(temps,mag,label='Magnetization')
plt.axvline(x=2.27, color='red', linestyle='--', label='Critical Temperature')
plt.xlabel('Temperature (K)')
plt.ylabel('Magnetization')
plt.title('Magnetization as a Function fo Temperature - 2D Ising Model')
plt.legend()
plt.savefig('part1plot.png')
plt.show()








