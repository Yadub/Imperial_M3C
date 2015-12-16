#Project part 4
import numpy as np
import matplotlib.pyplot as plt
import adv2d
#gfortran -c advmodule2d.f90
#f2py -llapack -c advmodule2d.f90 fdmoduleB.f90 fdmodule2dp.f90 ode2d.f90 -m adv2d --f90flags='-fopenmp' -lgomp


#add functions as needed



if _name_=='__main__':
    #add code as needed