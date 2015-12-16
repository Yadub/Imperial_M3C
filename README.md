# README #

This README details all files in the project repo.

## What is this repository for? ##

* Final Project for the [M3C course](http://imperialhpsc.bitbucket.org) at Imperial College London 2015-2016

## Contact Details ##

* Email: [yrb13@ic.ac.uk](mailto:yrb13@ic.ac.uk)

## Files and Brief Descriptions ##

Files that have been provided for the project but have not been editted in any way whatsoever have also been listed here. They state unchanged next to them.
Files that have been provided for the project but needed to be editted have originals copied with a "_template" suffix. They also state unchanged next to them

* contributors.txt: contains my name, CID, and a statement about this is project being my own unaided work.

* part1/: sub-directory containing the files relevant to "Part1: Improving your serial gradient code".
* part2/: sub-directory containing the files relevant to "Part2: Applying your parallel gradient code".
* part3/: sub-directory containing the files relevant to "Part3: Parallel advection equation solver".
* part4/: sub-directory containing the files relevant to "Part4: Solving the two-dimensional advection equation".

### part1/ ###

* data.in: used as input for testing the fortran code
* fdmodule.f90: unchaged
* fdmodule_dev.f90: copy of fdmodule.f90 that I am making working changes to
* fdmodule2d_dev.f90: copy of fdmodule2d from hw4soln that I am making working changes to
* p1_3.py: unchanged
* p1_3.so: shared ojbect file created for python using f2py --f90flags='-fopenmp' -lgomp -c fdmodule_dev.f90 fdmodule2d_dev.f90 -m p1_3 -llapack
* p1_dev.py: copy of p1_3.py that I have made working changes to
* p1_fig1.png: plot produced by p1_dev.py for speedup (serial vs OMP parallel code)
* p1_fig2.png: plot produced by p1_dev.py for speedup2d (serial vs 2D-serial code)
* testfortran: bash file with commands to compile and run the fortran code for the 2D-serial code with input from data.in
* testing.exe: executable created using the fortran modules and test function in fdmodule2d_dev.f90

### part2/ ###

* fdmodule.f90: unchanged
* fdmodule2d.f90: from solution form hw4. We use the grad_omp function along with fdmodule.f90 for image processing
* p2.png: image created by threshold function in p2_dev.py of the thresholded matrix
* p2.so: shared object file created for python using f2py  --f90flags='-fopenmp' -lgomp -c fdmodule.f90 fdmodule2d.f90 -m p2 -llapack
* p2_dev.py: copy of p2_template.py that I am making working changes to
* p2_template.py: unchanged

### part3/ ###

* adv_mpi.exe: execulable produced from ode_mpi_dev.f90
* advmodule.f90: unchanged
* data_mpi.in: file for reading in values when running adv_mpi.exe
* fdmoduleB.f90: unchanged
* fmpi.dat: python friendly data file produced upon running adv_mpi.exe
* ode_dev.f90: copy of ode_template.f90 that I am making working changes to
* ode_mpi_dev.f90: copy of ode_mpi_template.f90 that I am making working changes to
* ode_mpi_template.f90: unchanged
* ode_template.f90: unchanged
* p3.py: unchanged
* p3_dev: copy of p3_v2.py that I am making working changes to
* p3_v2.py: unchanged

### part4/ ###

* advmodule2d.f90: unchanged
* fdmoduleB.f90: unchanged
* ode2d.f90: unchanged
* p4.py: unchanged

## Helpful Links for the admin's benefit ##
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)
