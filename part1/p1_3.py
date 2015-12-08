"""
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

    
    
def test_gradN(numthreads):
    """input: number of threads used in test_grad_omp
    call test_grad1 with a range of grid sizes and
    assess speedup
    """
 
        
   
    


if __name__ == "__main__":
    
    e,ep,e2d,t,tp,t2d = test_grad1(400,200,2)
    s,s2d = test_gradN(2) 
    plt.show()
