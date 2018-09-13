EXE = rmcdhf_pitest
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS} -L/up/perand/kod/FEAST/3.0/lib/dx64 -lpfeast_sparse -lpfeast -L/opt/intel/composerxe/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm

APP_OBJ= \
     alcsca.o consis.o csfwgt.o dacon.o dampck.o dampor.o defcor.o \
     dsubrs.o eigen.o endsum.o engout.o estim.o fco.o gco.o getaldwt.o \
     getoldwt.o in.o iq.o ispar.o itjpo.o lagcon.o maxarr.o newe.o orbout.o \
     orthor.o orthy.o out.o outbnd.o prtrsl.o prwf.o setcof.o setdbg.o \
     setham.o setmcp.o setmix.o setsum.o setxuv.o setxv.o setxz.o solve.o \
     strsum.o xpot.o ypot.o \
     cofpot.o getaldmpi.o getoldmpi.o getscdpmpi.o hmout.o improv.o maneig.o \
     matrix.o newco.o rscfvu.o scf.o setcsl.o setlag.o mpiu.o setisompi.o \
     setrwfmpi.o lodrwfmpi.o lodcshpmpi.o lodcsh2pmpi.o iniest3p.o

$(EXE): $(APP_OBJ)
	$(FCP) -openmp -o $(BINFILE) -I/opt/intel/composerxe/mkl/include  $(APP_OBJ) $(APP_LIBS)  
#$(FC) -o $(BINFILE) $(FC_LD) -openmp -I/opt/intel/composerxe/mkl/include  $(APP_OBJ) $(APP_LIBS)  \
#                  -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FCP) -c $(FC_FLAGS) -openmp -I $(GRASP_INCLUDES) -I/opt/intel/composerxe/mkl/include $< -o $@

clean:
	-rm -f *.o core 

