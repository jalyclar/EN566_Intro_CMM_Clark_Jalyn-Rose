import random 
import numpy as np 
import matplotlib.pyplot as plt

###Double check code follows instructions??
###Add class
###Add statment at top to run 
###Add sys arg 


n=100
m=10**4
rwy=np.zeros(n)
rwx=np.zeros(n)
xn=np.zeros(n)
msd=np.zeros(n)

for i in range(m):
    x=0
    y=0
    for j in range(2,n):
        r1=random.random()
        r2=random.random()
        if r1<0.5:
            x=x+1
        else:
            x=x-1
        if r2<0.5:
            y=y+1
        else:
            y=y-1
        xn[j]=xn[j]+x
        rwx[j]=rwx[j]+x**2
        rwy[j]=rwy[j]+y**2
for i in range(n):
    xn[j]=xn[j]/m
    rwy[i]=rwy[i]/m
    rwx[i]=rwx[i]/m
    msd[i]=rwx[i]+rwy[i]

plt.figure()
plt.plot(range(n),xn,label=r'$\langle x_n \rangle$')
plt.xlabel('Timestep (i)')
plt.ylabel('Displacement')
plt.title('2D Random Walk')
plt.legend()
plt.show()

plt.figure()
plt.plot(range(n),rwx,label=r'$\langle (x_n^2) \rangle$')
plt.xlabel('Timestep (i)')
plt.ylabel('Displacement')
plt.title('2D Random Walk')
plt.legend()
plt.show()

slope,intercept=np.polyfit(range(n),msd,1)
D=slope/2

plt.figure()
plt.plot(range(n),msd,label=r'$\langle r^2 \rangle$')
plt.plot(range(n),slope*np.array(range(n))+intercept,'--', color='pink',label=f'Diffusion Fit - D={slope}')
plt.plot()
plt.xlabel('Timestep (i)')
plt.ylabel('MSD')
plt.title('2D Random Walk - Diffusion')
plt.legend()
plt.show()