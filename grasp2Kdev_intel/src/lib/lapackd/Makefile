LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/liblapackd.a
LIBOBJ = $(OBJ)

OBJ = \
dlae2.o   dlaebz.o  dlaev2.o  dlagtf.o  dlagts.o  dlansp.o  dlanst.o  dlapy2.o\
dlarf.o   dlarfg.o  dlarnv.o  dlartg.o  dlaruv.o  dlascl.o  dlaset.o  dlasr.o\
dlasrt.o  dlassq.o  dopgtr.o  dopmtr.o  dorg2l.o  dorg2r.o  dspevx.o  dsptrd.o\
dstebz.o  dstein.o  dsteqr.o  dsterf.o

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

