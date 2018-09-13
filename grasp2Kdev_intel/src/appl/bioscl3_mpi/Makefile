EXE = rtransition_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = $(GRASPLIB)/graspBD.o
BINFILE = $(BINDIR)/$(EXE)          

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} -lmpiu ${NEWGRASPLIBS} 

APP_OBJ= biosclmpi.o fname.o ldcsl1mpi.o ldcsl2mpi.o lodcslm.o  \
         merg12mpi.o mrgcslmpi.o ldlbl1.o ldlbl2.o printaLS.o \
         readmixmpi.o iqr.o isparr.o itjpor.o jcupr.o jqsr.o lodrwfimpi.o \
         lodrwffmpi.o getosdmpi.o brkt.o getrmpmpi.o strsum.o osclmpi.o \
         connect.o ichkq1.o cpmix.o mctoutmpi_gg.o\
         alcnsa.o alcnta.o  mctinmpi.o csfm.o printa.o spme.o trsortmpi.o \
         bessj.o alclla.o alcnma.o setcslm.o angdatampi.o engout1.o testmix.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
              -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

