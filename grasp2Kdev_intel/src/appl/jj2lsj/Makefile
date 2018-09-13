.SUFFIXES: .f90 .mod

# executable :: jj2lsj
EXE = jj2lsj
BINDIR = ${GRASP}/bin
BINFILE = $(BINDIR)/$(EXE)

#   Define data types
VASTO = vast_kind_param_M.o

#   Define Commons
Commons = common.o jj2lsj_data_1_C.o jj2lsj_data_2_C.o jj2lsj_data_3_C.o

#  Define memory management module
Memory = memory_man.o

#  Define interface to routines from the library
Interface = lib_interface.o

#  Define routines from the library
Routines = \
   convrt.o convrt_double.o dracah.o factt.o  \
   getmixblock.o getyn.o ichop.o idigit.o iq.o \
   ispar.o itjpo.o ittk.o ixjtik.o jcup.o jqs.o  \
   lodcsl.o lodiso.o lval.o nine.o nine0.o openfl.o pack.o \
   packLS.o parsjl.o prsrcn.o prsrsl.o setcsla.o setiso.o \
   sixj.o sixj1.o sixj2.o sixj3.o sixj35.o sixj4.o sixj5.o \
   starttime.o stoptime.o


JJ2LSJ_OBJ = \
   jj2lsj_code.o jj2lsj2K.o


$(EXE): ${VASTO} ${Commons} ${Memory} ${Interface} ${Routines} $(JJ2LSJ_OBJ)
	$(F90) -o $(BINFILE) $(F90_LD) $(JJ2LSJ_OBJ) $(Commons)  $(Memory) \
	${Routines}

.f90.o:
	 $(F90) -c $(F90_FLAGS) $< -I . -o $@

.f90.mod:
	 $(F90) -c $(F90_FLAGS) $< -I . -o $@

clean:
	-rm -f *.o core *.mod

JJ2LSJ_SRC =  \
      jj2lsj_data_1.f90 jj2lsj_data_2.f90 jj2lsj_data_3.f90 \
      jj2lsj_code.f90 jj2lsj2K.f90


