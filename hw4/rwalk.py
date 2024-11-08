#!/opt/anaconda3/bin/python 

import random 
import numpy as np 
import matplotlib.pyplot as plt
import sys 

class randomwalk():

    def __init__(self):
        self.part=None
        self.n=100

        m=10**4
        rwy=np.zeros(self.n)
        self.rwx=np.zeros(self.n)
        self.xn=np.zeros(self.n)
        self.msd=np.zeros(self.n)

        for i in range(m):
            x=0
            y=0
            for j in range(2,self.n):
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
                self.xn[j]=self.xn[j]+x
                self.rwx[j]=self.rwx[j]+x**2
                rwy[j]=rwy[j]+y**2
        for i in range(self.n):
            self.xn[j]=self.xn[j]/m
            rwy[i]=rwy[i]/m
            self.rwx[i]=self.rwx[i]/m
            self.msd[i]=self.rwx[i]+rwy[i]

    def parse_arg(self):
        for arg in sys.argv[1:]:
            if "--part=" in arg:
                self.parts=arg.split('=')[1].split(',')

    def partone(self):
        if '1' in self.parts:
            plt.figure()
            plt.plot(range(self.n),self.xn,label=r'$\langle x_n \rangle$')
            plt.xlabel('Timestep (i)')
            plt.ylabel('Displacement')
            plt.title('2D Random Walk')
            plt.legend()
            plt.show()

            plt.figure()
            plt.plot(range(self.n),self.rwx,label=r'$\langle (x_n^2) \rangle$',color='red')
            plt.xlabel('Timestep (i)')
            plt.ylabel('Displacement')
            plt.title('2D Random Walk')
            plt.legend()
            plt.show()
    
    def parttwo(self):
        if '2' in self.parts:
            slope,intercept=np.polyfit(range(self.n),self.msd,1)
            D=slope/2

            plt.figure()
            plt.plot(range(self.n),self.msd,label=r'$\langle r^2 \rangle$',color='black')
            plt.plot(range(self.n),slope*np.array(range(self.n))+intercept,'--', color='pink',label=f'Diffusion Fit - D={slope}')
            plt.plot()
            plt.xlabel('Timestep (i)')
            plt.ylabel('MSD')
            plt.title('2D Random Walk - Diffusion')
            plt.legend()
            plt.show()

if __name__=="__main__":
    data=randomwalk()
    data.parse_arg()
    data.partone()
    data.parttwo()
