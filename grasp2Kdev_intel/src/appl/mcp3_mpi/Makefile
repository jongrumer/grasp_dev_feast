EXE = rangular_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

 
APP_LIBS = ${GRASPBLOCKD}  -L ${GRASPLIB}  -lmpiu ${NEWGRASPLIBS} 

APP_OBJ= \
         cxk.o genmcpmpi_${FC_MALLOC}.o mcpmpi_gg.o outsdampi.o \
         setdbgmpi.o setmcpmpi.o \
         fndbeg.o getinf.o setdbg.o setmcp.o setsda.o setsum.o settmp.o \
         snrc.o sort.o strsum.o talk.o
 
$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
	-L ${LAPACK_DIR} ${LAPACK_LIB} 

.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

