# A-Monte-Carlo-method-for-tracking-dust-properties-during-coagulation-in-protoplanetary-disks
This code is provided to reproduce the plots in the paper "A Monte Carlo method for tracking dust properties during coagulation in protoplanetary disks".
This code is a modified version of the open source code mcdust https://github.com/vicky1997/mcdust 


This readme file was generated on 2026-02-16 by Nerea Gurrutxaga

## GENERAL INFORMATION

Title of code: Monte Carlo code to track dust properties
Description: This code calculates the dust evolution of particles, conserving the global inventory of dust properties

Corresponding Author Information
- Name: Nerea Gurrutxaga
- ORCID: 0009-0008-3256-9564
- Institution: Max Planck Institute for Solar System Research
- Email: gurrutxaga@mps.mpg.de


## SHARING/ACCESS Code

Data of the paper is available in Zenodo. 
doi: 10.5281/zenodo.18673340
link: https://zenodo.org/uploads/18673340

## FILE OVERVIEW

Main files:
- src: all .F90 programs for this
- src_mcdust: all .F90 programs for original mcdust, adapted from https://github.com/vicky1997/mcdust
- setup: compilation and parameter information
- outputs: for data storage
- scripts: generate all figures from the main manuscript 
- obj: storage of temporary ".mod" and ".o" files

## Prerequisites

`gfortran` `hdf5-serial` `python`

To install the required software in Ubuntu (this requires root permissions):
`sudo apt-get install gfortran`
`sudo apt-get install libhdf5-serial-dev`

Python is not required to run the code. But if you want to use the routines to read/write data from the simulation, you will need a Python installation.

## To reproduce tests: 
All setup files are already ready
### Kernel tests: 
These simulations test the validity of the algorithm against three analytical solutions of the Smoluchowski equation. There are three different coagulation kernels: constant, linear, and product. We run each kernel test with 200 particles and 10000 particles (hr: high resolution). The setup files of the six different simulations:
kernel_constant/
kernel_constant_hr/
kernel_linear/
kernel_linear_hr/
kernel_product/
kernel_product_hr/

Example of how to run it:
make SETUP_FILE=kernel_constant
run kernel_constant setups/kernel_constant/setup.par

Run script/Fig3-4.py to reproduce Figures 3 and 4

### Global simulation test: 
These simulations test the validity of the algorithm in a global protoplanetary disk. We test the original mcdust versus the algorithm presented here with three different values of the merging parameter X=0.1, 0.01, 0.001. The setup files of the four different simulations:
global_disk_mcdust/
global_disk_X1/
global_disk_X2/
global_disk_X3/

Example of how to run it:

- make mcdust SETUP_FILE=global_disk_mcdust
run global_disk_mcdust setups/global_disk_mcdust/setup.par

- make SETUP_FILE=global_disk_X1
run global_disk_X1 setups/global_disk_X1/setup.par

Run script/Fig5.py to reproduce Figure 5

### Conservation of the global budget test: 
This simulation tests the validity of the algorithm for conserving the global budget of dust properties. The setup file of the simulation:
water_mcdust/
water/

Example of how to run it:

- make mcdust SETUP_FILE=water_mcdust
run water_mcdust setups/water_mcdust/setup.par

- make SETUP_FILE=water
run water setups/water/setup.par

Run script/Fig6.py to reproduce Figure 6

