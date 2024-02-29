Scripts to run all test cases defined in cppdefs.h and make some plots :

- test_real.sh : 
main script calling run_real.sh and plot_real.sh
variables : list of test cases
            root repository for croco (default is the same thant CVTK)
            Run directory for jobcomp, modified files
            MPI/OPENMP
            MAX_PROCS
Usage      : test_real.sh  [-h] [-n EXAMPLE] [-d ROOT_DIR] [-r RUNDIR] [-p PARALLEL] [-m MAX_PROC]
 -h               : help
 -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all
 -d ROOTDIR       : Root of the git repository, default : same as CVTK
 -r RUNDIR        : Run repository for cppdefs.hand modified .F files, default : Run
 -p PARALLEL      : Type of parallelism (MPI or OPENMP), default : no
 -m MAX_PROC      : Max number of cpus available, default : 1
 -s MPI_RUN       : mpirun command, default : mpirun

- run_real.sh : 
compile and run a given test_case, called from test_real or alone
variables : list of test cases
            root repository for croco (default is the same thant CVTK)
            Run directory for jobcomp, modified files
            MPI/OPENMP
            MAX_PROCS
Usage      : run_real.sh  [-h] [-n EXAMPLE] [-d ROOT_DIR] [-r RUNDIR] [-p PARALLEL] [-m MAX_PROC]
 -h               : help
 -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all
 -d ROOTDIR       : Root of the git repository, default : same as CVTK
 -r RUNDIR        : Run repository for cppdefs.h and modified .F files, default : Run
 -p PARALLEL      : Type of parallelism (MPI or OPENMP), default : no
 -m MAX_PROC      : Max number of cpus available, default : 1
 -s MPI_RUN       : mpirun command, default : mpirun
 
- plot_real.sh : 
compile and run a given test_case, called from test_real or alone
requirement : ghostscript and matlab
variables : list of test cases
            root repository for croco (default is the same thant CVTK)
            MPI/OPENMP
            MAX_PROCS
Usage      : plot_real.sh  [-h] [-n EXAMPLE] [-d ROOT_DIR]
 -h               : help
 -n EXAMPLE       : TEST name, as listed in cppdefs.h, default : all
 -d ROOTDIR       : Root of the git repository, default : same as CVTK
 
 Method :
 - copy jobcomp, TEST_CASES and cppdefs.h from the RUNDIR repository 
 - copy .h .F files from ROOTDIR/OCEAN and RUNDUR
 - Compile and Run
 - for test cases with input .nc files, if not there, they are created (attempt)     
 - LOG files stored in LOG
 - final pdf file is merged.pdf