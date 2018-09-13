EXE = rhfs_lsj
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASP}/lib/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${NEWGRASPLIBS}

APP_OBJ= rhfs_lsj.o 
$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD)  $(APP_OBJ) $(APP_LIBS) \
              -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@


clean:
	-rm -f *.o core 

