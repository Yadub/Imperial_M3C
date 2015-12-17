# -*- coding: utf-8 -*-
#Yadu Bhageria
#CID:00733164
#Project part 4
import numpy as np
import matplotlib.pyplot as plt
import adv2d
#gfortran -c advmodule2d.f90
#f2py -llapack -c advmodule2d.f90 fdmoduleB.f90 fdmodule2dp.f90 ode2d.f90 -m adv2d --f90flags='-fopenmp' -lgomp

#add functions as needed
def advection2d(tf,nt,nx,ny,c1=1.0,c2=1.0,S=0.0,display=False,par=0,numthreads=1):
    """solve advection equation, df/dt + c1 df/dx + c2 df/dy = S
    for x=0,dx,...,(nx-1)*dx, and for y=0,dy,...,(ny-1)*dy,
    returns f2(x,y,tf) which is a solution obtained using the fortran
    routine ode_rk4 from ode2d
    -    f(x,y,t=0) = exp(−100[(x−x0)2+(y−y0)2]), L = n*dx
    -    nt time steps are taken from 0 to tf
    -    The solutions are plotted if display is true
    -    User can choose to run the code in parallel by setting par = 1 and choosing numthreads
    """
    #initalize dx,dy,dt by the given nx,ny,nt,tf
    dx = 1.0/(nx)
    dy = 1.0/(ny)
    dt = tf/nt
    #initialize array of x and y values: [0,1)
    x = np.arange(0.0,nx*dx,dx)
    y = np.arange(0.0,ny*dy,dy)
    #create grids X and Y with varying x and y respectively
    X, Y = np.meshgrid(x, y)
    #set initial conditions
    x0 = 0.5
    y0 = 0.5
    f0 = np.exp(-100.0*(np.power(X-x0,2)+np.power(Y-y0,2)))
    
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
    
    f,time = adv2d.ode2d.rk4(0,f0,dt,nt)
    #f = adv2d.ode2d.euler(0,f0,dt,nt) #for using the Euler method instead of RK4 method
    
    #display resulting estimate with contour slices
    if display:
        plt.figure()
        plt.contourf(X,Y,f)
        plt.title('Yadu Bhageria, advection2d, with tf=%s, nt=%s, \n (nx,ny)=(%s,%s) and (c1,c2,S)=(%s,%s,%s)'%(tf,nt,nx,ny,c1,c2,S))
        plt.xlabel('X')
        plt.ylabel('Y')
        #plt.savefig('p4_fig.png') #uncomment if you want to save file
    return f, time
    
def test_advection2d(num_threads):
    """input: number of threads used by my OMP
    Calls advection2d with a range of values of nx, ny, and
    parameters, tf=1,c1=1,c2=1,S=1,dx=1/nx, dy=1/ny
    -    Computes and plots error, estimates convergence rate of error 
    -    Computes and plots speedup dependent on numthreads
    -    Computes and plots the effect of variation of dt on accuracy 
        for some n values with grid = (n,n)
    
    Note: There is no difference between accuracy in the parallel code and the serial code
    
    Trends:
    -    The rate of convergence of the error is 1 with respect to grid size N 
        where N is total grid points. This is obviously scewed somewhat when 
        we use grids that are not symmetric in size. ie nx > ny or vice versa
    -    As grid size increases speedup increases and a square grid increases in speedup faster
        than a non square grid. Also although speedup increases it does not tend to
        the value of numthreads due to a large part of the code - rk4 - not being parallized
    -   for small dt, the code looses accuracy. After dt>(~400) the code reaches
        optimum accuracy and stops gaining any more accuracy
    """
    #set problem parameters
    nxvalues = np.array([50,100,150,200,300,400])
    nyvalues = nxvalues
    nn = np.size(nxvalues)
    
    tf = 1.0
    c1 = 1.0
    c2 = 1.0
    S = 1.0
    nt = 500

    e = np.empty((nn,nn))
    ep = np.empty((nn,nn))
    speedup = np.empty((nn,nn))
    ntot = np.empty((nn,nn))
    
    #loop through 50<n1,n2<1600 storing speedup
    for i in enumerate(nxvalues):
        for j in enumerate(nyvalues):
            
            print "running with n1,n2=",i[1],j[1] #to keep the user informed of the codes progress
            
            #for each n, construct grid and exact solution            
            dx = 1.0/(i[1])
            dy = 1.0/(j[1])
            x = np.arange(0.0,i[1]*dx,dx)
            y = np.arange(0.0,j[1]*dy,dy)
            x0 = 0.5
            y0 = 0.5
            X, Y = np.meshgrid(x, y)
            #analytical solution for any value of tf, not just integer.
            fa = np.exp(-100.0*(np.power((np.mod(X-c1*tf,1.0)-x0),2)+np.power((np.mod(Y-c2*tf,1.0)-y0),2))) + S*tf
            
            #solve advection eqn and compute error and time
            f,time = advection2d(tf,nt,i[1],j[1],c1,c2,S)  
            fp,timep = advection2d(tf,nt,i[1],j[1],c1,c2,S,par=1,numthreads=num_threads)                                 
            e[i[0],j[0]] = np.mean(np.abs(f-fa))
            ep[i[0],j[0]] = np.mean(np.abs(fp-fa))
            speedup[i[0],j[0]] = time/timep
            ntot[i[0],j[0]] = i[1]*j[1]
    
    m,p = np.polyfit(np.log(ntot.reshape(nn**2)),np.log(e.reshape(nn**2)),1)
    
    #display results
    plt.figure()
    plt.loglog(ntot.reshape(nn**2),e.reshape(nn**2),'D')
    plt.loglog(ntot.reshape(nn**2),np.exp(p)*(ntot.reshape(nn**2))**m,'--')
    plt.legend(('computed','lsq linear fit'),loc='best')
    plt.xlabel('n')
    plt.ylabel('error')
    plt.title('Yadu Bhageria, test_advection2d, fig1, variation of error with n \n for c1=c2=S=1 and m=%.3g'%(-m))
    plt.axis('tight')
    plt.savefig('p4_fig1.png')
        
    #display results
    plt.figure()
    plt.plot(nxvalues,speedup[:,::1])
    plt.legend(('n2=50','100','150','200','300','400'),loc='best')
    plt.xlabel('nx')
    plt.ylabel('speedup')
    plt.title('Yadu Bhageria, test_advection2d, fig2, Speedup vs GridSize\n with numthreads=%s'%(num_threads))
    plt.axis('tight')
    plt.savefig('p4_fig2.png')    
    
    ntvalues = np.array([100,200,400,800]) 
    ent = np.empty((np.size(ntvalues),nn))
    
    
    for ntval in enumerate(ntvalues):
        print('Running for nt = %s'%(ntval[1]))        
        for n in enumerate(nxvalues):
                #for each n, construct grid and exact solution            
                dx = 1.0/(n[1])
                dy = 1.0/(n[1])
                x = np.arange(0.0,n[1]*dx,dx)
            	y = np.arange(0.0,n[1]*dy,dy)
                x0 = 0.5
                y0 = 0.5
                X, Y = np.meshgrid(x, y)
                #analytical solution for any value of tf, not just integer.
                fa = np.exp(-100.0*(np.power((np.mod(X-c1*tf,1.0)-x0),2)+np.power((np.mod(Y-c2*tf,1.0)-y0),2))) + S*tf
                
                print 'running for nx,ny=',n[1],n[1]
                #solve advection eqn and compute error and time for just the parallel equation  
                fp,timep = advection2d(tf,ntval[1],n[1],n[1],c1,c2,S,par=1,numthreads=num_threads)                                 
            	ent[ntval[0],n[0]] = np.mean(np.abs(fp-fa))      

    #display results
    plt.figure()
    plt.plot(ntvalues,ent[:,::1],'x-')
    plt.legend(('nx,ny=50','100','150','200','300','400'),loc='best')
    plt.xlabel('nt')
    plt.ylabel('error')
    plt.title('Yadu Bhageria, test_advection2d, fig3 \n variation of error with nt for c=S=1, dx=1/n along with varying grid sizes')
    plt.axis('tight')
    plt.savefig('p4_fig3.png')
    
    plt.figure()
    plt.plot(ntvalues[2:4],ent[2:4,::1],'x-')
    plt.legend(('nx,ny=50','100','150','200','300','400'),loc='best')
    plt.xlabel('nt')
    plt.ylabel('error')
    plt.title('Yadu Bhageria, test_advection2d, fig4 - zoomed in fig3 \n zoomed in version')
    plt.savefig('p4_fig4.png')
    
    return -m,speedup 
    
if __name__=='__main__':
    nx = 400
    ny = 200
    nt = 500
    tf = 1.0
    f,t = advection2d(tf,nt,nx,ny,S=1.0,display=True)
    fp,tp = advection2d(tf,nt,nx,ny,S=1.0,par=1,numthreads=8)
    print('Is there a difference between serial and parallel code accuracy?: %r'%((np.sum(f-fp)!=0)))
    print('speedup from serial to parallel is %s times with numthreads=8'%(t/tp))
    m,speedup = test_advection2d(2)
    print('m = %.5g'%(m))
    plt.show()
    #Values above can be editted to provide varying results. 
    #For example changing numthreads gives significantly different values for speedup