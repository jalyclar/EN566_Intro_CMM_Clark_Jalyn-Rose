import numpy as np
import matplotlib.pyplot as plt 
import random

###Add class
###Add statment at top to run 
###Add sys arg 

h=60
w=40
iterations=[1,1000,100000,1000000,10000000]

for iter in iterations:
    grid=np.zeros((w,h),dtype=int)
    grid[:,:h//3]=-1 #species A left 
    grid[:,2*h//3:]=1 #species B right 

    movement=[(-1,0),(0,1),(1,0),(0,-1)] #movement in x and y 

    trials=100
    na_tot=np.zeros(h)
    nb_tot=np.zeros(h)

    for i in range(iter):
        while True:
            x=random.randint(0,w-1)
            y=random.randint(0,h-1)
            if grid[x,y]!=0:
                gas_par=grid[x,y]
                break
        dy,dx=random.choice(movement)
        ny=y+dy
        nx=x+dx
        if 0<=ny<h and 0<=nx<w and grid[nx,ny]==0:
            grid[nx,ny]=gas_par
            grid[x,y]=0
    plt.figure(figsize=(9,6)) #Part one 
    plt.imshow(grid,cmap='RdBu',origin='upper')
    plt.colorbar(label='Species')
    plt.title(f'Diffusion of Gas Particles A and B at iter:{iter}')
    plt.show()

    na=np.zeros(h) #Part two
    nb=np.zeros(h)
    for x in range(h):
        col=grid[:,x]
        A=np.sum(col==-1)
        B=np.sum(col==1)
        na[x]=A/w
        nb[x]=B/w
    plt.figure(figsize=(9,6))
    plt.plot(np.arange(h),na,color='red',label=r'$n_A(x)$')
    plt.plot(np.arange(h),nb,color='blue',label=r'$n_B(x)$')
    plt.xlabel('(x)')
    plt.ylabel('Population Density')
    plt.title(f'Population Densities at iter:{iter}')
    plt.legend()
    plt.show()

    for t in range(trials): #Part three
        na_trial=np.zeros(h)
        nb_trial=np.zeros(h)
        for x in range(h):
            col=grid[:,x]
            A=np.sum(col==-1)
            B=np.sum(col==1)
            na_trial[x]=A/w
            nb_trial[x]=B/w
        na_tot=na_tot+na_trial
        nb_tot=nb_tot+nb_trial
    na_av=na_tot/trials
    nb_av=nb_tot/trials
    plt.figure(figsize=(9,6))
    plt.plot(np.arange(h),na_av,color='red',label=r'Average Density $n_A(x)$')
    plt.plot(np.arange(h),nb_av,color='blue',label=r'Average Density $n_B(x)$')
    plt.xlabel('(x)')
    plt.ylabel('Average Population Density')
    plt.title(f'Average Population Densities at iter:{iter}')
    plt.legend()
    plt.show()