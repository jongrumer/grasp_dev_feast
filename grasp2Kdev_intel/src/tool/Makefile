.SUFFIXES: .f90 .mod

LIBDIR = $(GRASP)/lib
BIN = $(GRASP)/bin
LIB = $(GRASPLIBS)
BD = $(LIBDIR)/graspBD.o
UTIL = rmixextract  rcsfblock rwfnmchfmcdf readrwf \
       rlevels rtabtrans1 rotate rwfnrotate \
       ris_lsj risosupport rcsfexcitation HF rwfnplot rcsfsplit rtablevels \
       rtabtransE1 rmixreduce rtabtrans2 rlevelseV

install: EXE 
	cp rsave $(GRASP)/bin; chmod +x $(GRASP)/bin/rsave

EXE : $(BIN)/rmixextract $(BIN)/rcsfblock  \
      $(BIN)/rwfnmchfmcdf $(BIN)/readrwf \
      $(BIN)/rlevels $(BIN)/rtabtrans1 $(BIN)/rotate \
      $(BIN)/rwfnrotate \
      $(BIN)/ris_lsj $(BIN)/risosupport \
      $(BIN)/rcsfexcitation $(BIN)/HF $(BIN)/rwfnplot $(BIN)/rcsfsplit \
      $(BIN)/rtablevels $(BIN)/rtabtransE1 $(BIN)/rmixreduce \
		$(BIN)/rtabtrans2\
		$(BIN)/rlevelseV\
   

OBJ : rmixextract.o rcsfblock.o  rwfnmchfmcdf.o readrwf.o\
      rlevels.o rtabtrans1.o rotate.o format_mix.o rlevelsj.o \
      rwfnrotate.o ris_lsj.o risosupport.o rcsfexcitation.o HF.o \
      rwfnplot.o rcsfsplit.o rtablevels.o rtabtransE1.o \
		rmixreduce.o rtabtrans2.o rlevelseV.o
      


$(BIN)/format_mix  : format_mix.o
	$(FC) -o $(BIN)/format_mix  $(FC_LD) -L$(LIBDIR) \
                 format_mix.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rmixextract: rmixextract.o
	$(FC) -o $(BIN)/rmixextract $(FC_LD) -L$(LIBDIR)  rmixextract.o \
                 $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rcsfblock: rcsfblock.o
	$(FC) -o $(BIN)/rcsfblock $(FC_LD)  -L$(LIBDIR) rcsfblock.o \
                 $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rwfnmchfmcdf: rwfnmchfmcdf.o
	$(FC) -o $(BIN)/rwfnmchfmcdf $(FC_LD) -L$(LIBDIR) rwfnmchfmcdf.o \
                 $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/readrwf: readrwf.o
	$(FC) -o $(BIN)/readrwf $(FC_LD) -L$(LIBDIR) readrwf.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rlevels: rlevels.o
	$(FC) -o $(BIN)/rlevels $(FC_LD) -L$(LIBDIR) rlevels.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rtabtrans1: rtabtrans1.o
	$(FC) -o $(BIN)/rtabtrans1 $(FC_LD) -L$(LIBDIR) rtabtrans1.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rlevelsj: rlevelsj.o
	$(FC) -o $(BIN)/rlevelsj $(FC_LD) -L$(LIBDIR) rlevelsj.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rotate: rotate.o
	$(FC) -o $(BIN)/rotate $(FC_LD) -L$(LIBDIR) rotate.o $(BD) $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rwfnrotate: rwfnrotate.o
	$(FC) -o $(BIN)/rwfnrotate  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                rwfnrotate.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/ris_lsj: ris_lsj.o
	$(FC) $(FC_FLAGS) -o $(BIN)/ris_lsj  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                ris_lsj.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/risosupport: risosupport.o
	$(FC) -o $(BIN)/risosupport  $(FC_LD) -L$(LIBDIR) -I $(GRASP_INCLUDES) \
                risosupport.o $(BD) $(LIB) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(BIN)/rcsfexcitation: rcsfexcitation.o
	$(FC) -o $(BIN)/rcsfexcitation rcsfexcitation.o 

$(BIN)/HF: HF.o
	$(FC) -o $(BIN)/HF HF.o 

$(BIN)/rwfnplot: rwfnplot.o
	$(FC) -o $(BIN)/rwfnplot rwfnplot.o 

$(BIN)/rcsfsplit: rcsfsplit.o
	$(FC) -o $(BIN)/rcsfsplit rcsfsplit.o 

$(BIN)/rtablevels: rtablevels.o
	$(FC) -o $(BIN)/rtablevels rtablevels.o 

$(BIN)/rtabtransE1: rtabtransE1.o
	$(FC) -o $(BIN)/rtabtransE1 rtabtransE1.o 

$(BIN)/rmixreduce: rmixreduce.o
	$(FC) -o $(BIN)/rmixreduce rmixreduce.o 

$(BIN)/rtabtrans2: rtabtrans2.o
	$(FC) -o $(BIN)/rtabtrans2 rtabtrans2.o 

$(BIN)/rlevelseV: rlevelseV.o
	$(FC) -o $(BIN)/rlevelseV $(FC_LD) -L$(LIBDIR) rlevelseV.o $(LIB) \
                 -L ${LAPACK_DIR} ${LAPACK_LIB}


.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

.f90.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

