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
    
    Inial code checks the input and thus makes the function more robust and user-friendly
    It makes sure the correct type of data inputted and if not it tells the user what is wrong
    """
    
    #check input
    assert np.size(M.shape) == 3, "Inputted matrix, M, must have 3 dimensions"
    assert M.shape[2] == 3, "Inputted matrix, M, must have 3 values in the third dimension for the RGB colours"
    assert type(fac) is int or type(fac) is float, "fac, a must be a number"
    assert fac>=0, "fac, a must be non-negative"
    assert type(display) == bool, "Display flag must be a boolean. i.e. Either True or False"
    assert type(savefig) == bool, "Savefig flag must be a boolean. i.e. Either True or False"
    
    m1 = np.shape(M)[0]
    m2 = np.shape(M)[1]
    
    #set fortran module variables    
    f2d.n1 = m2
    f2d.n2 = m1
    
    f2d.dx1=1.0
    f2d.dx2=1.0
    
    #initalize arrays for computation
    df_amp = np.zeros((m1,m2))
    df_max = [0,0,0]
    N = M
        
    #grad_omp takes (f) as input and outputs (df1,df2,df_amp,df_max)
    for i in range(3):   
        N_temp = M[:,:,i]
        (dum1,dum2,df_amp,df_max[i]) = f2d.grad_omp(M[:,:,i])
        N_temp[df_amp < fac*df_max[i]] = 0
        N[:,:,i] = N_temp
        
    #displays the figure if the flag if true    
    if display == True:
        plt.figure()
        plt.imshow(N)
        plt.title('Yadu Bhageria, threshold, thresholded matrix')
    
    #saves the figure if the flag if true
    if savefig == True:
        plt.figure()
        plt.imshow(N)
        plt.title('Yadu Bhageria, threshold, thresholded matrix')
        plt.savefig("p2.png")
        plt.close()
    
    return N,df_max

if __name__ == '__main__':
    M=misc.face()
    print "shape(M):",np.shape(M)
    plt.figure()
    plt.imshow(M)
    plt.show() #may need to close figure for code to continue
    N,dfmax=threshold(M,0.2,True,True) #Uncomment this line 
    plt.show()
    