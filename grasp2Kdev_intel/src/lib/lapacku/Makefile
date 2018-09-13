LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/liblapacku.a
LIBOBJ = $(OBJ)

FC_FLAGS=-O0

OBJ = \
xerbla.o slamch.o dlamc1.o  dlamc2.o  dlamc3.o  dlamc4.o  dlamc5.o  dlamch.o  ilaenv.o  lsame.o

$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
#	lorder $(LIBOBJ) | tsort | xargs ar cr $(LIBA) #
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)
	rm -f *.o

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core 

