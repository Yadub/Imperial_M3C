"""Project part 2
Use gradient thresholding to process image.
Gradient is computed using grad_omp routine in fdmodule2d
f2py  --f90flags='-fopenmp' -lgomp -c fdmodule.f90 fdmodule2d.f90 -m p2 -llapack
"""
from p2 import fdmodule2d as f2d
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt

def threshold(M,fac,display,savefig):
    """set all elements of M to zero where grad_amp < fac*max(grad_amp)
    """
 
    


if __name__ == '__main__':
    M=misc.face()
    print "shape(M):",np.shape(M)
    plt.figure()
    plt.imshow(M)
    plt.show() #may need to close figure for code to continue
#    N,dfmax=threshold(M,0.2,True,True) #Uncomment this line 
    plt.show()
    