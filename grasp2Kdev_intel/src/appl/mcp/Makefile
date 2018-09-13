EXE = rangular
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)


APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${NEWGRASPLIBS}

APP_OBJ= \
          cxk.o genmcp.o  mcp_gg.o outsda.o setmcp2.o \
          fndbeg.o getinf.o setdbg.o setmcp.o setsda.o setsum.o settmp.o \
          snrc.o sort.o strsum.o talk.o
$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

