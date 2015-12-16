"""Project part 3, solve the 1d advection equation in fortran"""
import numpy as np
import matplotlib.pyplot as plt
import adv

#gfortran -c fdmoduleB.f90 ode.f90 advmodule.f90
#The above step should generate three .mod files
#f2py -llapack -c fdmoduleB.f90 ode.f90 advmodule.f90 -m adv --f90flags='-fopenmp' -lgomp

def advection1f(nt,tf,n,dx,c=1.0,S=0.0,display=False,numthreads=4):
    """solve advection equation, df/dt + c df/dx = S
    for x=0,dx,...,(n-1)*dx, and returns f(x,tf),fp(x,tf),
    and f4(x,tf) which are solutions obtained using the fortran
    routines ode_euler, ode_euler_omp, ode_rk4
    -    f(x,t=0) = sin(2 pi x/L), L = n*dx
    -    nt time steps are taken from 0 to tf
    -    The solutions are plotted if display is true
    """
    #adv.ode.euler(t0 float, y0 rank-1 array, dt, float, nt int)
    
    #construct grid and initial condition, set time span for odeint
    x = np.arange(0.0,n*dx,dx)
    f0 = np.sin(2.0*np.pi*x/(n*dx))
    
    #initalize values in the fortran modules
    adv.fdmodule.dx = dx
    adv.fdmodule.n = n
    adv.advmodule.s_adv = S
    adv.advmodule.c_adv = c
    
    f = adv.ode.euler(0,f0,tf/float(nt),nt)
    fp = adv.ode.euler_omp(0,f0,tf/float(nt),nt,numthreads)
    f4 = adv.ode.rk4(0,f0,tf/float(nt),nt)
    
    #display results
    if display:
        plt.figure()
        plt.plot(x,f0,'b')
        plt.plot(x,f,'g')
        plt.plot(x,fp,'r')
        plt.plot(x,f4,'k')
        plt.xlabel('x')
        plt.ylabel('f')
        plt.title('Yadu Bhageria, advection1f \n advection eqn. solution for n,dx,c,S=%d,%2.3f,%2.1f,%2.1f' %(n,dx,c,S))
        plt.grid()
        plt.axis('tight')
        plt.legend(('f0','f','fp','f4'),loc='best') 
            
    return f,fp, f4
    
 
def test_advection1f(n):
    """compute scaled L1 errors for solutions with n points 
    produced by advection1f when c=1,S=1,tf=1,nt=16000, L=1"""
    #set conditions
    c = 1.0
    S = 1.0
    tf = 1.0
    nt = 16000
    L = 1.0
    dx = L/n

    #analytical f
    x = np.arange(0.0,n*dx,dx)
    fa = np.sin(2.0*np.pi*np.mod(x-c*tf,L)/L) + S*tf
    
    f,fp,f4 = advection1f(nt,tf,n,dx,c,S,False)
    
    e = np.mean(np.abs(fa-f))
    ep = np.mean(np.abs(fa-fp))
    e4 = np.mean(np.abs(fa-f4))
    
    return e,ep,e4  
    #return e,ep,e4 #errors from euler, euler_omp, rk4
        
"""This section is included for assessment and must be included as is in the final
file that you submit,"""
if __name__ == '__main__':
    n = 200
    nt = 320000
    tf = 1.21
    dx = 1.0/n
    f,fp,f4 = advection1f(nt,tf,n,dx,1,0,True)
    e,ep,e4 = test_advection1f(200)
    eb,epb,e4b = test_advection1f(400)
    print e,ep,e4
    print eb,epb,e4b
    print e/eb,ep/epb,e4/e4b
    plt.show()