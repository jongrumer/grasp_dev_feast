.SUFFIXES: .f90 .mod

EXE = rcsfzerofirst
BINDIR = ${GRASP}/bin
BINFILE = $(BINDIR)/$(EXE)
FC = gfortran
FC_FLAGS = -O2

APP_OBJ= \
       rcsfzerofirst.o

$(EXE): $(APP_OBJ)
	$(FC) -o $(BINFILE) $(APP_OBJ) 

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -o $@
clean:
	-rm -f *.o core 
