EXE = bioscl_energy
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = $(GRASPLIB)/graspBD.o
BINFILE = $(BINDIR)/$(EXE)


APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${NEWGRASPLIBS}

APP_OBJ= bioscl.o fname.o ldcsl1.o ldcsl2.o lodcslm.o merg12.o mrgcsl.o \
         readmix.o iqr.o isparr.o itjpor.o jcupr.o jqsr.o lodrwfi.o \
         lodrwff.o getosd.o brkt.o getrmp.o strsum.o oscl.o connect.o \
         alcnsa.o alcnta.o  mctin.o csfm.o printa.o spme.o trsort.o \
         bessj.o alclla.o alcnma.o setcslm.o angdata.o engout1.o mctout_gg.o \
         cpmix.o ldlbl1.o ldlbl2.o printaLS.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(FC_LD) $(APP_OBJ) $(APP_LIBS) \
                   -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

