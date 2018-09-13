EXE = rcsfmerge
#BINDIR = ../../../bin
#GRASPLIB = ../../../lib
#GRASPBLOCKD = ../../../lib/graspBD.o
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${GRASPLIBS}

APP_OBJ= \
       chkplt.o iqr.o isparr.o itjpor.o jcupr.o jqsr.o ldcsl1.o ldcsl2.o \
       merg12.o mrgcsl.o setdbg.o rcsfblock.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

