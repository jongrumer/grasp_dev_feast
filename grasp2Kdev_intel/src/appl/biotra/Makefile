EXE = rbiotransform
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} ${NEWGRASPLIBS}

APP_OBJ= biotr.o fname.o intrpqf.o intrpqi.o kapdata.o lodrwff.o lodrwfi.o \
        radpar.o rintff.o rintii.o tcsl.o brkt.o gets.o wrtmat.o  \
        biotr1.o radfile.o ifnmnx.o ielsum.o copvec.o invmat.o \
        ulla.o trpmat.o matml4.o scalve.o setvec.o pamtmt.o bndinv.o lulu.o \
        inprod.o genmcp.o mcpin.o getmix.o qqsort.o citrag.o tiinig.o \
        ti1tv.o vecsum.o prsym.o setcslb.o angdata.o mcpout_gg.o

$(EXE) : $(APP_OBJ)
	$(FC) $(FC_FLAGS) -o $(BINFILE) $(FC_LD)  $(APP_OBJ) $(APP_LIBS) \
             -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

