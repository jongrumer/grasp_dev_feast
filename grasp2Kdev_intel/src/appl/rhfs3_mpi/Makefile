# FC = mpif77

EXE = rhfs_mpi

BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

# MPI_INC = ${MPI_INCLUDES}
# MPI_FFLAGS = -L$(MPI_LIBDIR) $(MPI_LIBLST)

APP_LIBS = ${GRASPBLOCKD} -L ${GRASPLIB} -lmpiu ${NEWGRASPLIBS}

APP_OBJ= \
        chkplt.o gethfdmpi.o hfs_gg.o hfs92mpi.o matelt.o  rinthf.o setdbg.o \
        readmixmpi.o lodcslmpijb.o setcslajb.o setNCFjb.o lodcsljb.o \
        setsum.o strsum.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
	-L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o: 
	$(FC_MPI) -c $(FC_MPIFLAGS) -I $(GRASP_INCLUDES) $*.f

clean:
	-rm -f *.o core *.mod

