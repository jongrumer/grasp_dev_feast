#EXE = erwf_ORG
EXE = rwfnestimate
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
   chkplt.o     erwf.o       frmhyd.o     frmrwf.o     frmtfp.o     genrwf.o \
   getinf.o     prtrem.o     sbstep.o     screenpar.o  setdbg.o     setsum.o \
   solvh.o      strsum.o     summry.o     tail.o       tfpot.o      wrtrwf.o \
   canonical.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
             -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

