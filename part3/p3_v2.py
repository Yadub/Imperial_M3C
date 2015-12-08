"""Project part 3, solve the 1d advection equation in fortran"""
import numpy as np
import matplotlib.pyplot as plt
import adv

#gfortran -c fdmoduleB.f90 ode.f90 advmodule.f90
#The above step should generate three .mod files
#f2py -llapack -c fdmoduleB.f90 ode.f90 advmodule.f90 -m adv --f90flags='-fopenmp' -lgomp

def advection1f(nt,tf,n,dx,c=1.0,S=0.0,display=False,numthreads=1):
    """solve advection equation, df/dt + c df/dx = S
    for x=0,dx,...,(n-1)*dx, and returns f(x,tf),fp(x,tf),
    and f4(x,tf) which are solutions obtained using the fortran
    routines ode_euler, ode_euler_omp, ode_rk4
    -    f(x,t=0) = sin(2 pi x/L), L = n*dx
    -    nt time steps are taken from 0 to tf
    -    The solutions are plotted if display is true
    """
    

 
def test_advection1f(n):
    """compute scaled L1 errors for solutions with n points 
    produced by advection1f when c=1,S=1,tf=1,nt=16000, L=1"""
    
          
    return e,ep,e4 #errors from euler, euler_omp, rk4
        
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