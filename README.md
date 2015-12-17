# README #

This README details all files in the project repo.

## What is this repository for? ##

* Final Project for the [M3C course](http://imperialhpsc.bitbucket.org) at Imperial College London 2015-2016

## Contact Details ##

* Email: [yrb13@ic.ac.uk](mailto:yrb13@ic.ac.uk)

## Files and Brief Descriptions ##

Files that have been provided for the project but have not been editted in any way whatsoever have also been listed here. They state unchanged next to them.

Files that have been provided for the project but needed to be editted have originals copied with a "_template" suffix. They also state unchanged next to them

Each folder contains some *.mod, *.o files that have no been listed.

* contributors.txt: contains my name, CID, and a statement about this is project being my own unaided work.

* part1/: sub-directory containing the files relevant to "Part1: Improving your serial gradient code".
* part2/: sub-directory containing the files relevant to "Part2: Applying your parallel gradient code".
* part3/: sub-directory containing the files relevant to "Part3: Parallel advection equation solver".
* part4/: sub-directory containing the files relevant to "Part4: Solving the two-dimensional advection equation".

### part1/ ###

* data.in: used as input for testing the fortran code
* fdmodule.f90: final version for marking. Completed the subroutine, cfd4_2d
* fdmodule_dev.f90: copy of fdmodule.f90 that I was making working changes to
* fdmodule_template.f90: unchanged
* fdmodule2d.f90: final version for marking. Added grad2d and test_grad2d subroutines 
* fdmodule2d_dev.f90: copy of fdmodule2d from hw4soln that I am making working changes to
* p1_3: final version for marking. completed the script and functions
* p1_3.so: shared ojbect file created for python using f2py --f90flags='-fopenmp' -lgomp -c fdmodule_dev.f90 fdmodule2d_dev.f90 -m p1_3 -llapack
* p1_dev.py: copy of p1_3.py that I have made working changes to
* p1_fig1.png: plot produced by p1_dev.py for speedup (serial vs OMP parallel code)
* p1_fig2.png: plot produced by p1_dev.py for speedup2d (serial vs 2D-serial code)
* p1_3_template.py: unchanged
* testfortran: bash file with commands to compile and run the fortran code for the 2D-serial code with input from data.in
* testing.exe: executable created using the fortran modules and test function in fdmodule2d_dev.f90
* some *.mod, *.o files that have not been listed

### part2/ ###

* fdmodule.f90: unchanged
* fdmodule2d.f90: from solution form hw4. We use the grad_omp function along with fdmodule.f90 for image processing
* p2.png: image created by threshold function in p2_dev.py of the thresholded matrix
* p2.py: final version for marking. made changes so it outputs thresholded matrix, dfmax and possibly produces an image
* p2.so: shared object file created for python using f2py  --f90flags='-fopenmp' -lgomp -c fdmodule.f90 fdmodule2d.f90 -m p2 -llapack
* p2_dev.py: copy of p2_template.py that I am making working changes to
* p2_template.py: unchanged

### part3/ ###

* adv.so: shared object file created for python
* adv_mpi.exe: execulable produced from ode_pi_dev.f90
* advection_mpi.f90: final version for marking. developed euler_mpi to solve a 1D advection with a parallel distributed memory code
* advmodule.f90: unchanged
* data_mpi.in: file for reading in values when running adv_mpi.exe
* fdmoduleB.f90: unchanged
* fmpi.dat: python friendly data file produced upon running adv_mpi.exe
* ode.f90: final version for marking. Developed RHS, RHS_omp, euler_omp
* ode_dev.f90: copy of ode_template.f90 that I am making working changes to
* ode_mpi_dev.f90: copy of ode_mpi_template.f90 that I am making working changes to
* ode_mpi_template.f90: unchanged
* ode_template.f90: unchanged
* p3.py: final version for marking. Completed advection1f and testadvection1f 
* p3_dev: copy of p3_v2.py that I am making working changes to
* p3_template.py: unchanged
* p3_v2_template.py: unchanged
* some *.mod, *.o files that have not been listed

### part4/ ###

* adv2d.so: shared object file created for python
* advmodule2d.f90: unchanged
* fdmodule2dp: copied from solutions for hw4 part 1. Editted to call fd2 instead of cfd4 and removed variables df_amp, df_max
* fdmoduleB.f90: unchanged
* ode2d.f90: final version for marking. modified RHS
* ode2d_dev.f90: copy of ode_template.f90 that I am making working changes to
* ode2d_template.f90: unchanged
* p4.py: final version for marking. developed functions that compute and analyze the 2D advection equation. 
* p4_fig.png: produced from p4.py and the advection2d function
* p4_fig1.png: produced from p4.py and the test_advection2d function: grid size vs error
* p4_fig2.png: produced from p4.py and the test_advection2d function: grid size vs speedup
* p4_fig3.png: produced from p4.py and the test_advection2d function: dt vs error
* p4_fig4.png: produced from p4.py and the test_advection2d function: zoomed in version of fig3
* p4_dev.py: copy of p4_template that I am making working changes to
* p4_template.py: unchanged
* some *.mod, *.o files that have not been listed

## Helpful Links for the admin's benefit ##
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)