EXE = sms
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASP}/lib/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
         chkplt.o densmcp.o densnew.o gco.o getsmd.o polint.o rintdens.o \
         rintiso.o setdbg.o setmcp.o setsum.o sms.o sms92.o smsmcp.o \
         smsnew.o strsum.o vinti.o getmixblock.o
$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD)  $(APP_OBJ) $(APP_LIBS) \
              -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@


clean:
	-rm -f *.o core 

