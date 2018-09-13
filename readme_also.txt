From: Per Jönsson
Sent: Friday, October 2, 2015 10:02 AM
To: GaigalasG
Subject: GRASP with new diagonalization routines
 
Hello Gediminas,

In 

../../scratch/tspejo/GRASP_FEAST/till_PerJ

on cluster 1 you will find the two directories FEAST and grasp2Kdev_intel
This is the version of grasp that now utilizes the fast diagonalization routines.

You have to start to compile FEAST with the intel compiler. Start by removing all
files in the  FEASTS/lib/ifort and FEAST/lib/lx64 directories

Then go to the src directory and issue the command
make ARCH=ifort LIB=feast all

you should now get the libfeast.a  libfeast_banded.a  libfeast_dense.a  libfeast_sparse.a
libraries. They reside in FEASTS/lib/ifort but copy the also to FEASTS/lib/ifort

Now start to compile the grasp package using intel.
Maybe you have to edit the path to the ifort compiler in the file
make-environment_ifort as to make things work on your system.
Then do

source  ./make-environment_ifort 

and compile as usual. We have a problem on our cluster and I don't know how you handled this on your cluster. In the make-environment_ifort file there is

FC_MPI="mpif90"  

By default this command calls gfortran (see note from intel) 



However, instead of calling gfortran we want to call ifortran. According to intel 

https://software.intel.com/en-us/forums/intel-clusters-and-hpc-technology/topic/288354

we should write

FC_MPI="mpiifort" 

instead. However, mpiifort this is not installed on our cluster (hopefully it is installed on yours) and we are working on it. This is as far as I have come. 

 If you are able to compile everything properly with ifort you should then go into rscf_test and compile this separately. It is here Per Andersson has done changes so that FEAST is called instead of dvdson (llok at the Makefile and in maneig where Per has done changes).

I hope to get this working soon and then we can together evaluate it to see if it speeds up calculations.

Per Andersson will come to Malmö in two weeks and then I can ask him more but I wanted to get started already now so that I know if there are any problems.

Best wishes
Per


 
#/bin/bash
# -------------------------------------------------------------------------------------------------------------------------------------
# GRASP2K ENVIRONMENT FLAGS
# -------------------------------------------------------------------------------------------------------------------------------------
#
# Define the following global variables according to your environment and
# source this script or add these definitions to your terminal configuration
# file, eg. ~/.cshrc, ~/.bashrc or ~/.profile.
#
# Current version: Linux, ifort
#
# Note: jjgen can't be compiled by ifort and is therefore compiled with gcc/gfortran
#
# -------------------------------------------------------------------------------------------------------------------------------------
# Set up main flags
# -------------------------------------------------------------------------------------------------------------------------------------
export   FC="ifort"                                            # Fortran compiler
export   FCP="mpif90"                                          # Fortran compiler
export   FC_FLAGS=" -g -traceback -C -O3 -xHost -parallel -save -zero -align all -W0 -mcmodel=large" # Serial code compiler flag
export   FC_LD=" "                                             # Serial linker flags
export   GRASP="${PWD}"                                        # Location of the grasp2k root directory
export   LAPACK_LIB="-llapackd -llapacku -lblas"               # Library to be searched
export   LAPACK_DIR="${GRASP}/lib"                             # Location of LAPACK library
export   GRASPLIBS="-l92 -lnjgraf -ldvdson"                    # Libraries to be searched by v1 and v2 code
export   NEWGRASPLIBS="-lrang -l92 -ldvdson"                   # Libraries to be searched by v3 code
export   FC_MALLOC="LINUX"                                     # Memory allocation routine (available for Linux and other systems)
export   GRASP_INCLUDES="${GRASP}/src/lib/def"                 # Location of parameter definition file
export   F90="ifort"                                           # Fortran compiler
export   F90_FLAGS=" -mcmodel=large"                           # Serial code compiler flag
export   CPP="g++"                                             # C++ compiler
export   CPP_FLAGS="-O3"                                       # C++ compiler flags
export   CPP_LD="-static"                                      # C++ linker
# -------------------------------------------------------------------------------------------------------------------------------------
# For running MPI
# -------------------------------------------------------------------------------------------------------------------------------------
export   FC_MPI="mpiifort"                                     # MPICH/mpif90 compiler
#export   FC_MPI="mpif90"                                      # MPICH/mpif90 compiler
export   FC_MPIFLAGS="-O2 -save -zero -align all -W0 -mcmodel=large"          # Parallel code compiler flags
export   FC_MPILD="-O -mcmodel=large"                          # MPI Loader
export   MPI_TMP="/tmp/$USER"                                  # Temporary directory for storage of files written to disk
