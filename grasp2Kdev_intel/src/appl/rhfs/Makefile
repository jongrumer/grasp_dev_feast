EXE = rhfs
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASP}/lib/graspBD.o
BINFILE = $(BINDIR)/$(EXE)
# GRASP_INCLUDES = ${GRASP}/src/lib/def

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${NEWGRASPLIBS}

APP_OBJ= \
        chkplt.o engouth.o gethfd.o hfs_gg.o hfs92.o matelt.o  rinthf.o setdbg.o \
        setsum.o strsum.o getmixblock.o angdata.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
               -L ${LAPACK_DIR} ${LAPACK_LIB} 

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@


clean:
	-rm -f *.o core 


