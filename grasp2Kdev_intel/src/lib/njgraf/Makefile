LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libnjgraf.a
LIBOBJ = $(OBJ) 

OBJ = \
 bubble.o  cut1l.o diagrm.o lolpop.o otherj.o  printj.o  settab.o  triang.o \
 change.o  cut2l.o dracah.o neibor.o phase.o   rdiag.o   sprate.o  var.o \
 chklp1.o  cutnl.o gensum.o njgraf.o phase2.o  search.o  square.o  way.o \
 chvar.o   delta.o intab.o  ordtri.o polygn.o  setdm.o   trdel.o   zero.o

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

