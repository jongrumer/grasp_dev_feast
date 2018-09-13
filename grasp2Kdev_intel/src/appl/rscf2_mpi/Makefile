EXE = rmcdhf_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L ${GRASPLIB}  -lmpiu ${GRASPLIBS} 

APP_OBJ= \
     alcsca.o consis.o csfwgt.o dacon.o dampck.o dampor.o defcor.o \
     dsubrs.o eigen.o endsum.o engout.o estim.o fco.o gco.o getaldwt.o \
     getoldwt.o in.o iq.o ispar.o itjpo.o lagcon.o maxarr.o newe.o orbout.o \
     orthor.o orthy.o out.o outbnd.o prtrsl.o prwf.o setcof.o setdbg.o \
     setham.o setmcp.o setmix.o setsum.o setxuv.o setxv.o setxz.o solve.o \
start.o     strsum.o xpot.o ypot.o \
     cofpotmpi.o getaldmpi.o getoldmpi.o getscdmpi.o hmoutmpi.o  \
     improvmpi.o maneigmpi.o matrixmpi.o newcompi.o  rscfmpivu_${FC_MALLOC}.o scfmpi.o \
     setcslmpi.o setdbgmpi.o setlagmpi.o

$(EXE): $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
	-L ${LAPACK_DIR} ${LAPACK_LIB} 

.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

