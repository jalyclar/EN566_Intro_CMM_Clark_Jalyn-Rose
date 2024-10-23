import numpy as np 
import matplotlib.pyplot as plt

#In 2d
def gauss_seidel(DEBUG=True):
    #Variables 
    a=0.6 #inital distance seperation
    h=0.1 #grid spacing 
    grid_size=10 #setting grid size
    tolerance=10**-3
    e=8.854*10**-12 #F/m

    #Gauss Seidel 
    x=np.linspace(-10,10,grid_size) #R boundaries 
    y=np.linspace(-10,10,grid_size)
    X,Y=np.meshgrid(x,y) #for contour plot 
    V=np.zeros((grid_size,grid_size)) #Initializing grid 
    rho=np.zeros((grid_size,grid_size)) #For modified caluclation 
    center_grid=grid_size//2 #finding the center of the grid
    if DEBUG==True:
        print(f"center of grid: {center_grid}\n\n")
    dipole_distance=int(a/h) #converting initial distance seperation into grid points
    if DEBUG==True:
        print(f"dipole distance: {dipole_distance}\n\n")
    q_pos_x=center_grid
    q_pos_y=center_grid-dipole_distance//2
    q_neg_x=center_grid
    q_neg_y=center_grid+dipole_distance//2
    if DEBUG==True:
        print(f"q_pos_x: {q_pos_x}\n\n")
        print(f"q_pos_y: {q_pos_y}\n\n")
        print(f"q_neg_x: {q_neg_x}\n\n")  
        print(f"q_neg_y: {q_neg_y}\n\n")     
    V[q_pos_x,q_pos_y]=1 
    V[q_neg_x,q_neg_y]=-1
    rho[q_pos_x,q_pos_y]=1 
    rho[q_neg_x,q_neg_y]=-1
    if DEBUG==True:
        print(f"Initial V(r) grid: {V}\n\n")
    iterations=int((grid_size**2)/2)
    for iter in range(iterations):
        max_change=0 #for convergence criteria 
        for i in range(1, grid_size-1): #rows
            for j in range(1, grid_size-1): #columns 
                if (i,j)==(q_pos_x,q_pos_y) or (i,j)==(q_neg_x,q_neg_y): #skipping past intial dipole point charges 
                    continue 
                r=np.sqrt((i-center_grid)**2+(j-center_grid)**2)
                if r<10: #boundary condition 
                    if rho[i,j]==0: #using laplace where there is no cahrge density 
                        V_n=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])
                    else: #using poisson where there is charge density 
                        V_n=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])+(rho[i,j]*h**2/e)
                    del_V=np.absolute(V_n-V[i,j]) #max chnage in V
                    V[i,j]=V_n
                    if del_V>max_change: 
                        max_change=del_V #assinging max chnage in V to variable 
        if max_change<tolerance:
            print(f'Converged after {iter} iterations.\n\n')
            break #when max chnage is less than totelerence stops while loop 
    if DEBUG==True:
        print(f"Final V(r) grid: {V}")
    plt.contour(X,Y,V, levels=np.linspace(-1,1,20),colors='black',linewidths=0.3)
    plt.contourf(X,Y,V,levels=np.linspace(-1,1,20),cmap='RdBu',vmin=-1,vmax=1)
    plt.title('Electric Potential of a Static Electric Dipole: Gauss-Seidel')
    plt.xlabel('(m)')
    plt.ylabel('(m)')
    plt.colorbar() 
    plt.show()

def jacobi_relax(DEBUG=True):
    #Variables 
    a=0.6 #inital distance seperation
    h=0.1 #grid spacing 
    grid_size=10 #setting grid size
    e=8.854*10**-12 #F/m

    #jacobi relaxation 
    x=np.linspace(-10,10,grid_size) #R boundaries 
    y=np.linspace(-10,10,grid_size)
    X,Y=np.meshgrid(x,y) #for contour plot 
    V=np.zeros((grid_size,grid_size)) #Initializing grid 
    new_V=np.zeros((grid_size,grid_size)) #For updated grid
    rho=np.zeros((grid_size,grid_size)) #For modified caluclation 
    center_grid=grid_size//2 #finding the center of the grid
    if DEBUG==True:
        print(f"center of grid: {center_grid}\n\n")
    dipole_distance=int(a/h) #converting initial distance seperation into grid points
    if DEBUG==True:
        print(f"dipole distance: {dipole_distance}\n\n")
    q_pos_x=center_grid
    q_pos_y=center_grid-dipole_distance//2
    q_neg_x=center_grid
    q_neg_y=center_grid+dipole_distance//2
    if DEBUG==True:
        print(f"q_pos_x: {q_pos_x}\n\n")
        print(f"q_pos_y: {q_pos_y}\n\n")
        print(f"q_neg_x: {q_neg_x}\n\n")  
        print(f"q_neg_y: {q_neg_y}\n\n")     
    V[q_pos_x,q_pos_y]=1 
    V[q_neg_x,q_neg_y]=-1
    rho[q_pos_x,q_pos_y]=1 
    rho[q_neg_x,q_neg_y]=-1
    if DEBUG==True:
        print(f"Initial V(r) grid: {V}\n\n")
    iterations=grid_size**2
    print(f'Number of iterations for convergence: {iterations}\n\n')
    tolerance=10**-3

    for iter in range(iterations):
        converge=0
        for i in range(1,grid_size-1):
            for j in range(1,grid_size-1):
                if (i,j)==(q_pos_x,q_pos_y) or (i,j)==(q_neg_x,q_neg_y): #skipping past intial dipole point charges 
                    continue 
                r=np.sqrt((i-center_grid)**2+(j-center_grid)**2)
                if r<10: #boundary condition
                    if rho[i,j]==0: #using laplace where there is no cahrge density 
                        V[i,j]=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])
                    else: #using poisson where there is charge density 
                        V[i,j]=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])+(rho[i,j]*h**2/e)
                    del_V=np.absolute(new_V[i,j]-V[i,j]) #max chnage in V
                    if del_V>converge:
                        converge=del_V
        new_V[:,:]=V[:,:]
        if converge<tolerance:
            print(f'Converged after {iter} iterations.\n\n')
            break                           
    if DEBUG==True:
        print(f"Final V(r) grid: {new_V}")
    plt.contour(X,Y,new_V, levels=np.linspace(-1,1,20),colors='black',linewidths=0.3)
    plt.contourf(X,Y,new_V,levels=np.linspace(-1,1,20),cmap='RdBu',vmin=-1,vmax=1)
    plt.title('Electric Potential of a Static Electric Dipole: Jacobi Relaxation')
    plt.xlabel('(m)')
    plt.ylabel('(m)')
    plt.colorbar() 
    plt.show()

if __name__=="__main__":
    gauss_seidel()
    jacobi_relax()




