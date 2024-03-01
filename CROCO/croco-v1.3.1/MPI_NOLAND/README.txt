=======================================================================
              USING MPP_PREP to suppress LAND-ONLY processors
=======================================================================

MPP_PREP is a fortran utility to determine the optimal MPI decomposition to supress land computation, when one cpu contains only land points. It was originaly developped by Maurice Imbard and Jean-Marc Molines for NEMO model.

Content:
- Makefile : an example of makefile
- mpp_optimiz.f90 : fortran code
- namelist
- mpp_plot.py : python script for diagnostic

Compilation
===========

- Edit Makefile to put your own compiler, options and NetCDF location
- gmake

An executable called mpp_optimiz is generated

Settings
========

Edit the namelist. Note :
- &NAMPARAM section is not used for now
- Npts, the number of ghostcells used by CROCO depends on the numerical choices. Set to 2 in general, to 3 with advection schemes of order 5 or 6
- jprocx is the maximum number of CPU you will be able to use
Copy (or link) your grid file in MPI_NOLAND directory

Execution
=========
Run ./mpp_optimiz
It will output :
- on the screen : the optimum decomposition to suppress the maximum of land points, with values for NNODES, NP_XI and NP_ETA. NNODES is the number of cpu to use, NP_XI and NP_ETA (NP_XIxNP_ETA <=NNODES) the cartesian decomposition
- a file called processor.layout with the description of all combinations tested
- a COVDATA file called YOURCONFIG-NXxNY_NTOT : use for diagnostic only

Visualisation
=============
Launch mpp_plot.py
Usage : mpp_plot.py GRID_FILE COVDTA_FILE
It will display the decomposition with sea processors. Note that the figure diplayed may not be exactly accurate for graphical reasons.

How to use in CROCO
===================
- edit param.h to change values of NP_XI, NP_ETA and NNODES to the one provided by mpp_optimiz
- if your grid file is not called croco_grd.nc (for exampel located in CROCO_FILES/croco_grd.nc) 
  - edit MPI_Setup.F and change line 136 the name of the grid file ( from croco_grd.nc to CROCO_FILES/croco_grd.nc for exemple) 
  - or link your grid file in the launching directory as ln -sf CROCO_FILES/mycroco_frd.nc croco_grd.nc
- activate  the cpp key MPI_NOLAND
- compile as usual
- launch on a number of CPU equal to NNODES


