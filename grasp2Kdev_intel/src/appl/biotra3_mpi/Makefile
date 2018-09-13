EXE = rbiotransform_mpi
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)

APP_LIBS = ${GRASPBLOCKD} -L${GRASPLIB} -lmpiu ${NEWGRASPLIBS}

APP_OBJ= biotrmpi.o fname.o intrpqf.o intrpqi.o kapdata.o \
        lodrwffmpi.o lodrwfimpi.o mcpoutmpi_gg.o mcpinmpi.o\
        radparmpi.o rintff.o rintii.o tcsl.o brkt.o getsmpi.o wrtmat.o  \
        biotr1.o radfilempi.o ifnmnx.o ielsum.o copvec.o invmat.o \
        ulla.o trpmat.o matml4.o scalve.o setvec.o pamtmt.o bndinv.o lulu.o \
        inprod.o getmixmpi.o qqsortmpi.o citragmpi.o tiinigmpi.o \
        ti1tv.o vecsum.o prsym.o setcslbmpi.o setcslampi.o angdatampi.o

$(EXE) : $(APP_OBJ)
	$(FC_MPI) -o $(BINFILE) $(FC_MPILD) $(APP_OBJ) $(APP_LIBS) \
                -L ${LAPACK_DIR} ${LAPACK_LIB}

.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

