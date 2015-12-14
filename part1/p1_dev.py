"""
Yadu Bhageria
CID: 00733164
Test gradient routines in fdmodule2d
To build the .so module:
rm *.mod
gfortran -c fdmodule.f90 fdmodule2d.f90
f2py  --f90flags='-fopenmp' -lgomp -c fdmodule.f90 fdmodule2d.f90 -m p1_3 -llapack
"""
from p1_3 import fdmodule as f1
from p1_3 import fdmodule2d as f2d
import numpy as np
import matplotlib.pyplot as plt

def test_grad1(n1,n2,numthreads):
    """call fortran test routines, test_grad,test_grad_omp,teat_grad2d
    with n2 x n1 matrices, run f2d_test_grad_omp with numthreads threads
    return errors and computation times
    """
    #set fortran module variables    
    f2d.n1=n1
    f2d.n2=n2

    f2d.dx1=1.0/(f2d.n1-1)
    f2d.dx2=1.0/(f2d.n2-1)

    f2d.numthreads = numthreads
    #-----------

    tarray = []
    tparray = []    
    t2darray = []

    #call fortran test routines ten times and return error and average times    
    for i in range(10):
        e,t = f2d.test_grad()
        ep,tp=f2d.test_grad_omp()
        e2d,t2d = f2d.test_grad2d()

        tarray += [t]
        tparray += [tp]
        t2darray += [t2d]
    

    return e,ep,e2d,np.mean(tarray),np.mean(tparray),np.mean(t2darray)

def test_gradN(numthreads):
    """input: number of threads used in test_grad_omp
    call test_grad1 with a range of grid sizes and
    assess speedup
    """

    n1values = np.array([50, 100, 200, 400, 800, 1600, 3200])
    n2values = n1values
    nn = np.size(n1values)
    speedup = np.empty((nn,nn))
    speedup2d = np.empty((nn,nn))

    #loop through 50<n1,n2<1600 storing speedup
    for i in enumerate(n1values):
        for j in enumerate(n2values):
            print "running with n1,n2=",i[1],j[1]              
            f2d.n1=i[1]
            f2d.n2=j[1]
    
            e,ep,e2d,t,tp,t2d = test_grad1(f2d.n1,f2d.n2,2)
            speedup[i[0],j[0]] = t/tp #compute and store speedup
            speedup2d[i[0],j[0]] = t/t2d
            
    #display results
    plt.figure()
    plt.plot(n1values,speedup[:,::2])
    plt.legend(('n2=50','200','800','3200'),loc='best')
    plt.xlabel('n1')
    plt.ylabel('speedup')
    plt.title('Yadu Bhageria, test_gradN')
    plt.axis('tight')
    plt.savefig("p1_fig1.png")
    
    plt.figure()
    plt.plot(n1values,speedup2d[:,::2])
    plt.legend(('n2=50','200','800','3200'),loc='best')
    plt.xlabel('n1')
    plt.ylabel('speedup2d')
    plt.title('Yadu Bhageria, test_gradN')
    plt.axis('tight')
    plt.savefig("p1_fig2.png")   
    
    return speedup,speedup2d
    
if __name__ == "__main__":
    
    e,ep,e2d,t,tp,t2d = test_grad1(400,200,2)
#    print(e,ep,e2d,t,tp,t2d)
    s,s2d = test_gradN(2) 
    plt.show()
