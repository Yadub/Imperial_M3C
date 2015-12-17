# -*- coding: utf-8 -*-
#Yadu Bhageria
#CID:00733164
#-*- coding: utf-8 -*-
#Project part 4
import numpy as np
import matplotlib.pyplot as plt
import adv2d
#gfortran -c advmodule2d.f90
#f2py -llapack -c advmodule2d.f90 fdmoduleB.f90 fdmodule2dp.f90 ode2d.f90 -m adv2d --f90flags='-fopenmp' -lgomp


#add functions as needed
def advection2d(nt,tf,nx,ny,c1=1.0,c2=1.0,S=0.0,display=False,par=0,numthreads=1):
    """solve advection equation, df/dt + c1 df/dx + c2 df/dy = S
    for x=0,dx,...,(nx-1)*dx, and for y=0,dy,...,(ny-1)*dy,
    returns f2(x,y,tf) which is a solution obtained using the fortran
    routine ode_rk4 from ode2d
    -    f(x,y,t=0) = exp(−100[(x−x0)2+(y−y0)2]), L = n*dx
    -    nt time steps are taken from 0 to tf
    -    The solutions are plotted if display is true
    """
    dx = 1.0/(nx)
    dy = 1.0/(ny)
    
    #set fortran module variables    
    adv2d.fdmodule2d.n1 = nx
    adv2d.fdmodule2d.n2 = ny
    adv2d.fdmodule2d.dx1 = dx
    adv2d.fdmodule2d.dx2 = dy
    adv2d.advmodule.c1_adv = c1
    adv2d.advmodule.c2_adv = c2
    adv2d.advmodule.s_adv = S
    adv2d.ode2d.numthreads = numthreads
    adv2d.ode2d.par = par
    
    x = np.arange(0.0,nx*dx,dx)
    y = np.arange(0.0,ny*dy,dy)

    X, Y = np.meshgrid(x, y)
    x0 = 0.5
    y0 = 0.5
    dt = tf/nt
    
    f0 = np.exp(-100.0*(np.power(X-x0,2)+np.power(Y-y0,2)))
    
    f,time = adv2d.ode2d.rk4(0,f0,dt,nt)
    #f = adv2d.ode2d.euler(0,f0,dt,nt)
    
    fa = np.exp(-100.0*(np.power((np.mod(X-c1*tf,1.0)-x0),2)+np.power((np.mod(Y-c2*tf,1.0)-y0),2))) + S*tf
    
    error = np.sum(np.abs(fa-f))/(nx*ny)
    
    #display results
    if display:
        plt.figure()
        #plt.contourf(X,Y,f0)
        CP = plt.contourf(X,Y,fa)
        plt.clabel(CP, inline=1, fontsize=10)
        plt.contourf(X,Y,f)
        plt.title('Yadu Bhageria, advection2d with tf = %s'%(tf))
    
    return error, time
    
if __name__=='__main__':
    nx = 50
    ny = 50
    nt = 800
    tf = 1.0
    e,t = advection2d(nt,tf,nx,ny,display=True)
    ep,tp = advection2d(nt,tf,nx,ny,par=1,numthreads=8)
    print(e, ep)
    print('Serial code took %s seconds. Speedup from serial to parallel is %s times'%(t,t/tp))
    plt.show()