import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

#In 2d

################     WHats wrong with jacobiiiiiiiiii !!!!!!!   ##########################
################     when WHILE statment instead for for range iterations ???????    ##########################
################     Adjust structure of code     ##########################
################    Allow command line input      ##########################

def gauss_seidel(DEBUG=False):
    #Variables 
    a=0.6 #inital distance seperation
    h=0.1 #grid spacing 
    grid_size=10 #setting grid size
    tolerance=10**-6
    e=8.854*10**-12 #F/m
    Q=1

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
    V[q_pos_x,q_pos_y]=Q
    V[q_neg_x,q_neg_y]=-Q
    rho[q_pos_x,q_pos_y]=Q
    rho[q_neg_x,q_neg_y]=-Q
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
                    max_change+=del_V #assinging max chnage in V to variable 
        if max_change<tolerance:
            print(f'Converged after {iter} iterations.\n\n')
            break #when max chnage is less than totelerence stops while loop 
    if DEBUG==True:
        print(f"Final V(r) grid: {V}")
        
    #Equipotential Lines
    plt.contour(X,Y,V, levels=np.linspace(np.min(V),np.max(V),20),cmap='inferno')
    #plt.contourf(X,Y,V,levels=np.linspace(-1,1,20),cmap='RdBu',vmin=-1,vmax=1)
    plt.title('Electric Potential of a Static Electric Dipole: Gauss-Seidel')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar() 
    plt.grid()
    plt.show()

    #V(r)
    r = np.sqrt(X**2 + Y**2)
    V_r = V[int(grid_size / 2), :]  # Potential along y=0 line
    plt.figure(figsize=(10, 6))
    plt.plot(x, V_r)
    plt.title('Electric Potential V(r) as a function of distance from the origin')
    plt.xlabel('Distance r')
    plt.ylabel('Electric Potential V(r)')
    plt.grid()
    plt.xlim(-10, 10)
    plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
    plt.show()

def jacobi_relax(DEBUG=False,grid_size=10,tolerance=10**-3):
    #Variables 
    a=0.6 #inital distance seperation
    h=0.1 #grid spacing 
    e=8.854*10**-12 #F/m
    Q=1

    #jacobi relaxation 
    x=np.linspace(-10,10,grid_size) #R boundaries 
    y=np.linspace(-10,10,grid_size)
    X,Y=np.meshgrid(x,y) #for contour plot 
    #V=np.zeros((grid_size,grid_size)) #Initializing grid 
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
    new_V[q_pos_x,q_pos_y]=Q
    new_V[q_neg_x,q_neg_y]=-Q
    rho[q_pos_x,q_pos_y]=Q
    rho[q_neg_x,q_neg_y]=-Q
    if DEBUG==True:
        print(f"Initial V(r) grid: {new_V}\n\n")
    #iterations=grid_size**2
    #if DEBUG==True:
        #print(f'Number of iterations for convergence: {iterations}\n\n')
    iterations=1000000
    N_iter=0
    for iter in range(iterations):
        old_V=new_V.copy()
        del_V=0
        for i in range(1,grid_size-1):
            for j in range(1,grid_size-1):
                if (i,j)==(q_pos_x,q_pos_y) or (i,j)==(q_neg_x,q_neg_y): #skipping past intial dipole point charges 
                    continue
                r=np.sqrt((i-center_grid)**2+(j-center_grid)**2) #Checking to see where point is relative to center
                if r<10: #boundary condition
                    if rho[i,j]==0: #using laplace where there is no cahrge density 
                        new_V[i,j]=1/4*(old_V[i+1,j]+old_V[i-1,j]+old_V[i,j+1]+old_V[i,j-1])
                    else: #using poisson where there is charge density 
                        new_V[i,j]=1/4*(old_V[i+1,j]+old_V[i-1,j]+old_V[i,j+1]+old_V[i,j-1])+(rho[i,j]*h**2/e)
                    #del_V=np.absolute(new_V[i,j]-V[i,j]) #max chnage in V
                    #change=max(change,del_V)
                    #if DEBUG==True:
                        #print(f'Current delta(V): {change:.8f}\n\n')
                    del_V+=np.absolute(new_V[i,j]-old_V[i,j])
        N_iter=iter
        #del_V=np.max(np.absolute(new_V[i,j]-old_V[i,j]))
        #print(f'{new_V}\n\n')
        if DEBUG==True:
            print(f'Current iteration: {N_iter}\n\n')
        if iter>5 and del_V<tolerance:
            break
    #print(f'Converged after {iter} iterations.\n\n')
    #print(f"Final V(r) grid: {new_V}\n\n")
    return X,Y,new_V,N_iter

def sor(DEBUG=False,grid_size=10,tolerance=10**-3):
    #Variables 
    a=0.6 #inital distance seperation
    h=0.1 #grid spacing 
    e=8.854*10**-12 #F/m
    alpha=2/(1+(np.pi/grid_size))
    Q=1
    if DEBUG==True:
        print(f'alpha value: {alpha}\n\n')
    x=np.linspace(-10,10,grid_size) #R boundaries 
    y=np.linspace(-10,10,grid_size)
    X,Y=np.meshgrid(x,y) #for contour plot 
    V=np.zeros((grid_size,grid_size)) #Initializing grid gauss
    V_new=np.zeros((grid_size,grid_size)) #Initializing grid for SOR
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
    V[q_pos_x,q_pos_y]=Q #For updates using gauss
    V[q_neg_x,q_neg_y]=-Q
    V_new[q_pos_x,q_pos_y]=Q #For final update for sor 
    V_new[q_neg_x,q_neg_y]=-Q
    rho[q_pos_x,q_pos_y]=Q 
    rho[q_neg_x,q_neg_y]=-Q
    if DEBUG==True:
        print(f"Initial V(r) grid: {V}\n\n")
    #iterations=grid_size
    #if DEBUG==True:
        #print(f'Number of iterations for convergence: {iterations}\n\n')
    N_iter=0
    #for iter in range(iterations):
    iterations=1000000
    for iter in range(iterations):
        change=0
        for i in range(1,grid_size-1):
            for j in range(1,grid_size-1):
                if (i,j)==(q_pos_x,q_pos_y) or (i,j)==(q_neg_x,q_neg_y): #skipping past intial dipole point charges 
                    continue 
                r=np.sqrt((i-center_grid)**2+(j-center_grid)**2)
                if r<10: #boundary condition
                    if rho[i,j]==0: #using laplace where there is no cahrge density 
                        V_new[i,j]=(1-alpha)*V_new[i,j]+alpha*(1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1]))
                        #V_n=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])
                    else: #using poisson where there is charge density 
                        V_new[i,j]=(1-alpha)*V_new[i,j]+alpha*(1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])+(rho[i,j]*h**2/e))
                        #V_n=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])+(rho[i,j]*h**2/e)
                    #V_new[i,j]=(1-alpha)*V_new[i,j]+alpha*V_n
                    #print(f"{V_new}\n\n")
                    del_V=np.absolute(V_new[i,j]-V[i,j]) #max chnage in V
                    #print(del_V)
                   #V[i,j]=V_new[i,j]
                    change=max(change,del_V)
                V[i,j]=V_new[i,j]
        N_iter=iter
        #print(f'Change and iter step: {change:.8f} and {N_iter}')
        if change<tolerance:
            break
    if DEBUG==True:
        print(f"Final V(r) grid: {V_new}\n\n")
    #print(f'Converged after {N_iter} iterations.\n\n')
    #plt.contour(X,Y,V_new, levels=np.linspace(-1,1,50),cmap='inferno')
    #plt.title('Electric Potential of a Static Electric Dipole: Jacobi Relaxation')
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.colorbar() 
    #plt.show()                        
    return X,Y,V,N_iter

def partone(DEBUG=True):
    grid_size=10
    X,Y,new_V,N_iter=jacobi_relax(DEBUG=False,grid_size=grid_size,tolerance=10**-3)
    if DEBUG==True:
        print(f"Final V(r) grid: {new_V}\n\n")
    print(f'Converged after {N_iter} iterations.\n\n')
    plt.contour(X,Y,new_V, levels=np.linspace(-1,1,50),cmap='inferno') #RdYlBu
    #plt.contourf(X,Y,new_V,levels=np.linspace(-1,1,50),cmap='inferno',vmin=-1,vmax=1)
    plt.title('Electric Potential of a Static Electric Dipole: Jacobi Relaxation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar()
    plt.grid()
    plt.show()

    #V(r)
    center_grid=grid_size//2
    h=0.1
    r = np.linspace(-grid_size/2,grid_size/2,grid_size)*h #r distances along y=0 line 
    V_r = new_V[center_grid,:]  # Potential along y = 0 line

    plt.figure(figsize=(10, 6))
    plt.plot(r, V_r)
    plt.title('Electric Potential V(r) as a function of distance from the origin')
    plt.xlabel('Distance r')
    plt.ylabel('Electric Potential V(r)')
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
    plt.show()

def partwo(DEBUG=False):
    error_limits=np.logspace(-11,-1,num=20)
    N_iter_values=[]
    errors=[]
    for tolerance in error_limits:
        errors.append(tolerance)
        _,_,_,N_iter=jacobi_relax(tolerance=tolerance)
        N_iter_values.append(N_iter)
    if DEBUG==True:
        print(f'N iter values: {N_iter_values}\n\n')
        print(f'Error values: {errors}\n\n')
    plt.figure(figsize=(8,5))
    plt.plot(errors,N_iter_values,label='Jacobi Relaxation Model')
    plt.xscale('log')
    plt.title(r'Relationship Between The Number of Iterations ($N_{\text{iter}}$) and Tolerance Limit ($\epsilon$)')
    plt.ylabel(r'$N_{\text{iter}}$')
    plt.xlabel(r'$\epsilon$ (log-scale)')
    plt.legend()
    plt.show()

def partthree(DEBUG=True):
    grids=np.linspace(10,60,5,dtype=int) 
    iter_jacob=[]
    iter_sor=[]
    if DEBUG==True:
        print(f'grid sizes to use: {grids}')
    for L in grids:
        _,_,_,N_iter=jacobi_relax(False,L,10**-3)
        _,_,_,N_iter_s=sor(False,L,10**-3)
        iter_jacob.append(N_iter)
        iter_sor.append(N_iter_s)
    if DEBUG==True:
        print(f'N iter values Jacobi: {iter_jacob}\n\n')
        print(f'N iter values SOR : {iter_sor}\n\n')
    plt.figure()
    plt.plot(grids,iter_jacob,label='Jacobi Relaxation Method')
    plt.plot(grids,iter_sor,label='Simultaneous Over-Relaxation Method')
    plt.xlabel('Grid Size')
    plt.ylabel(r'$N_{\text{iter}}$')
    plt.title(r'Iterations ($N_{\text{iter}}$) dependence on Grid Size')
    plt.legend()
    plt.show()


if __name__=="__main__":
    #gauss_seidel()
    #jacobi_relax()
    #sor()
    #partone()
    partwo()
    #partthree()




