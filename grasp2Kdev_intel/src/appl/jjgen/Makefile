EXE = rcsfgenerate
BINDIR = ${GRASP}/bin
GRASPLIB = ${GRASP}/lib
GRASPBLOCKD = ${GRASPLIB}/graspBD.o
BINFILE = $(BINDIR)/$(EXE)
FC = gfortran
FC_FLAGS = -fno-automatic

APP_OBJ= \
       adder.o blanda.o blandc.o blandb.o fivefirst.o fivelines.o \
       genb.o jjgen15b.o kopp1.o kopp2.o lasa1.o lasa2.o lika.o \
       lockad.o matain.o matcin.o matbin.o merge.o mergeb.o \
       reada.o reffa.o slug.o sluggo.o test.o copy7t9.o open79.o \
       rcsfblock.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(APP_OBJ) -L$(GRASPLIB) $(GRASPLIBS) 

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@
clean:
	-rm -f *.o core 
