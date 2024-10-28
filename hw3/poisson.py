#!/opt/anaconda3/bin/python 

import numpy as np 
import matplotlib.pyplot as plt
import sys

#In 2d

################    Allow command line input      ##########################

class poisson:
    def __init__(self):
        self.a=0.6 #inital distance seperation
        self.h=0.1 #grid spacing 
        self.e=8.854*10**-12 #F/m
        self.Q=1
        self.grid_size=11 #setting grid size

        self.x=np.linspace(-10,10,self.grid_size) #R boundaries
        self.y=np.linspace(-10,10,self.grid_size)
        self.X,self.Y=np.meshgrid(self.x,self.y) #for contour plot
        self.V=np.zeros((self.grid_size,self.grid_size)) #Initializing grid 
        self.rho=np.zeros((self.grid_size,self.grid_size)) #For modified caluclation
        self.new_V=np.zeros((self.grid_size,self.grid_size)) #For updated grid
        self.center_grid=self.grid_size//2 #finding the center of the grid
        self.dipole_distance=int(self.a/self.h) #converting initial distance seperation into grid points

        self.q_pos_x=self.center_grid
        self.q_pos_y=self.center_grid-self.dipole_distance//2
        self.q_neg_x=self.center_grid
        self.q_neg_y=self.center_grid+self.dipole_distance//2

        self.V[self.q_pos_x,self.q_pos_y]=self.Q
        self.V[self.q_neg_x,self.q_neg_y]=self.Q*-1
        self.rho[self.q_pos_x,self.q_pos_y]=self.Q
        self.rho[self.q_neg_x,self.q_neg_y]=self.Q*-1
        self.part=None

    def parse_argv(self):
        for arg in sys.argv[1:]:
            if "part=" in arg:
                self.part = arg.split('=')[1]

    def jacobi_relax(self,DEBUG=False,tolerance=10**-20):
        if self.part=='partone' or 'parttwo':
            new_V=self.V.copy()
            rho=self.rho
            if DEBUG==True:
                print(f"center of grid: {self.center_grid}\n\n")
            if DEBUG==True:
                print(f"dipole distance: {self.dipole_distance}\n\n")
            if DEBUG==True:
                print(f"q_pos_x: {self.q_pos_x}\n\n")
                print(f"q_pos_y: {self.q_pos_y}\n\n")
                print(f"q_neg_x: {self.q_neg_x}\n\n")  
                print(f"q_neg_y: {self.q_neg_y}\n\n")     
            if DEBUG==True:
                print(f"Initial V(r) grid: {new_V}\n\n")
            iterations=int(self.grid_size**2)
            if DEBUG==True:
                print(f'Number of iterations for convergence: {iterations}\n\n')
            for iter in range(iterations):
                old_V=new_V.copy()
                del_V=0
                for i in range(1,self.grid_size-1):
                    for j in range(1,self.grid_size-1):
                        if (i,j)==(self.q_pos_x,self.q_pos_y) or (i,j)==(self.q_neg_x,self.q_neg_y): #skipping past intial dipole point charges 
                            continue
                        r=np.sqrt((i-self.center_grid)**2+(j-self.center_grid)**2) #Checking to see where point is relative to center
                        if r<10: #boundary condition
                            if rho[i,j]==0: #using laplace where there is no cahrge density 
                                new_V[i,j]=(1/4)*(old_V[i+1,j]+old_V[i-1,j]+old_V[i,j+1]+old_V[i,j-1])
                            else: #using poisson where there is charge density 
                                new_V[i,j]=(1/4)*(old_V[i+1,j]+old_V[i-1,j]+old_V[i,j+1]+old_V[i,j-1])+(rho[i,j]*self.h**2/self.e)
                            #del_V=np.absolute(new_V[i,j]-V[i,j]) #max chnage in V
                            #change=max(change,del_V)
                            #if DEBUG==True:
                                #print(f'Current delta(V): {change:.8f}\n\n')
                            del_V+=np.absolute(new_V[i,j]-old_V[i,j])
                N_iter_j=iter
                #del_V=np.max(np.absolute(new_V[i,j]-old_V[i,j]))
                #print(f'{new_V}\n\n')
                if DEBUG==True:
                    print(f'Current iteration: {N_iter_j}\n\n')
                if iter>5 and del_V<tolerance:
                    if DEBUG==True:
                        print(f'Converged after {N_iter_j} iterations.\n\n')
                    break
            jacobi_V=new_V
            return jacobi_V,N_iter_j

    def partone(self,DEBUG=False):
        if self.part=='partone':
            jacobi_V,_=self.jacobi_relax()
            if DEBUG==True:
                print(f"Final V(r) grid: {jacobi_V}\n\n")
            plt.figure()
            plt.contour(self.X*self.h,self.Y*self.h,jacobi_V,levels=np.linspace(-1,1,50),cmap='inferno') #RdYlBu
            #plt.contourf(X,Y,new_V,levels=np.linspace(-1,1,50),cmap='inferno',vmin=-1,vmax=1)
            plt.title('Electric Potential of a Static Electric Dipole: Jacobi Relaxation')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.colorbar()
            plt.grid()
            plt.show()

            #V(r)
            r=np.linspace(-self.grid_size,self.grid_size,self.grid_size)*self.h #r distances along y=0 line 
            V_r=jacobi_V[self.center_grid,:]  # Potential along y = 0 line

            plt.figure(figsize=(10, 6))
            plt.plot(r, V_r)
            plt.title('Electric Potential V(r) as a function of distance from the origin')
            plt.xlabel('Distance r')
            plt.ylabel('Electric Potential V(r)')
            plt.grid()
            plt.xlim(np.min(r),np.max(r))
            plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
            plt.show()

    def partwo(self,DEBUG=False):
        if self.part=='parttwo':
            error_limits=np.logspace(-9,-1,num=20)
            N_iter_values=[]
            errors=[]
            for tolerance in error_limits:
                errors.append(tolerance)
                _,N_iter_j=self.jacobi_relax(tolerance=tolerance)
                N_iter_values.append(N_iter_j)
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

    def sor(self,DEBUG=False,tolerance=10**-20):
        if self.part=='partthree':
            alpha=2/(1+(np.pi/self.grid_size))
            if DEBUG==True:
                print(f'alpha value: {alpha}\n\n')
            if DEBUG==True:
                print(f"center of grid: {self.center_grid}\n\n")
            if DEBUG==True:
                print(f"dipole distance: {self.dipole_distance}\n\n")
            if DEBUG==True:
                print(f"q_pos_x: {self.q_pos_x}\n\n")
                print(f"q_pos_y: {self.q_pos_y}\n\n")
                print(f"q_neg_x: {self.q_neg_x}\n\n")  
                print(f"q_neg_y: {self.q_neg_y}\n\n")     
            V=self.V.copy() #For updates using gauss
            V_new=V.copy() #For final update for sor 
            rho=self.rho
            if DEBUG==True:
                print(f"Initial V(r) grid: {V}\n\n")
            iterations=self.grid_size
            if DEBUG==True:
                print(f'Number of iterations for convergence: {iterations}\n\n')
            for iter in range(iterations):
                change=0
                for i in range(1,self.grid_size-1):
                    for j in range(1,self.grid_size-1):
                        if (i,j)==(self.q_pos_x,self.q_pos_y) or (i,j)==(self.q_neg_x,self.q_neg_y): #skipping past intial dipole point charges 
                            continue 
                        r=np.sqrt((i-self.center_grid)**2+(j-self.center_grid)**2)
                        if r<10: #boundary condition
                            if rho[i,j]==0: #using laplace where there is no cahrge density 
                                V_new[i,j]=(1-alpha)*V_new[i,j]+alpha*(1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1]))
                                #V_n=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])
                            else: #using poisson where there is charge density 
                                V_new[i,j]=(1-alpha)*V_new[i,j]+alpha*(1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])+(rho[i,j]*self.h**2/self.e))
                                #V_n=1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])+(rho[i,j]*h**2/e)
                            #V_new[i,j]=(1-alpha)*V_new[i,j]+alpha*V_n
                            #print(f"{V_new}\n\n")
                            del_V=np.absolute(V_new[i,j]-V[i,j]) #max chnage in V
                            #print(del_V)
                        #V[i,j]=V_new[i,j]
                            change=max(change,del_V)
                        V[i,j]=V_new[i,j]
                #print(f'Change and iter step: {change:.8f} and {N_iter}')
                N_iter_s=iter
                if change<tolerance:
                    break
            self.sor_V=V_new
            if DEBUG==True:
                print(f"Final V(r) grid: {V_new}\n\n")
            #print(f'Converged after {N_iter} iterations.\n\n')
            #plt.contour(X,Y,V_new, levels=np.linspace(-1,1,50),cmap='inferno')
            #plt.title('Electric Potential of a Static Electric Dipole: Jacobi Relaxation')
            #plt.xlabel('x')
            #plt.ylabel('y')
            #plt.colorbar() 
            #plt.show()
            return N_iter_s

    def partthree(self,DEBUG=False):
        if self.part=='partthree':
            grids=np.linspace(10,30,3,dtype=int) 
            iter_jacob=[]
            iter_sor=[]
            if DEBUG==True:
                print(f'grid sizes to use: {grids}')
            for L in grids:
                self.grid_size = L
                self.x = np.linspace(-10, 10, self.grid_size)
                self.y = np.linspace(-10, 10, self.grid_size)
                self.X, self.Y = np.meshgrid(self.x, self.y)
                self.V = np.zeros((self.grid_size, self.grid_size))
                self.rho = np.zeros((self.grid_size, self.grid_size))
                # Reinitialize dipole positions and charges
                self.center_grid = self.grid_size // 2
                self.dipole_distance = int(self.a / self.h)
                self.q_pos_x = self.center_grid
                self.q_pos_y = self.center_grid - self.dipole_distance // 2
                self.q_neg_x = self.center_grid
                self.q_neg_y = self.center_grid + self.dipole_distance // 2
                self.V[self.q_pos_x, self.q_pos_y] = self.Q
                self.V[self.q_neg_x, self.q_neg_y] = -self.Q
                self.rho[self.q_pos_x, self.q_pos_y] = self.Q
                self.rho[self.q_neg_x, self.q_neg_y] = -self.Q

                _,N_iter_j=self.jacobi_relax()
                N_iter_s=self.sor()
                iter_jacob.append(N_iter_j)
                iter_sor.append(N_iter_s)
            if DEBUG==True:
                print(f'N iter values Jacobi: {iter_jacob}\n\n')
                print(f'N iter values SOR : {iter_sor}\n\n')
            fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
            ax1.plot(grids,iter_jacob,label='Numerical Data',color='pink')
            ax1.plot(grids,np.square(grids),label=r'Expected Data $N_{\text{iter}} \propto n^2 $',linestyle='--',color='black',alpha=0.3)
            ax1.set_xlabel('Grid Size (L)')
            ax1.set_ylabel(r'$N_{\text{iter}}$')
            ax1.set_title('Jacobi Relaxation Method ')
            ax2.plot(grids,iter_sor,label='Simultaneous Over-Relaxation Method',color='lightblue')
            ax2.plot(grids,grids,label=r'Expected Data $N_{\text{iter}} \propto n$',linestyle='--',color='black',alpha=0.3)
            ax2.set_xlabel('Grid Size (L)')
            ax2.set_ylabel(r'$N_{\text{iter}}$')
            ax2.set_title('Simultaneous Over-Relaxation Method')
            ax1.legend()
            ax2.legend()
            plt.suptitle(r'Iterations ($N_{\text{iter}}$) dependence on Grid Size')
            plt.show()                   

if __name__=="__main__":
    data=poisson()
    data.parse_argv()
    data.jacobi_relax()
    data.partone()
    data.partwo()
    data.sor()
    data.partthree()




