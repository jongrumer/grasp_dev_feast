EXE = rcsfdelete
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= chkplt.o iqr.o isparr.o itjpor.o jcupr.o \
         jqsr.o ldcsl1.o ldcsl2.o xcld12.o xcsl.o \
         setdbg.o length.o rcsfblock.f


$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD)  $(APP_OBJ) $(APP_LIBS) \
                -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

