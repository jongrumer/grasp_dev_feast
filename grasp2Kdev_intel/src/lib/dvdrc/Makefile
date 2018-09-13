LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libdvdson.a
LIBOBJ = $(OBJ)

OBJ = dvdson.o iniest.o gdvd.o
SRC = dvdson.f iniest.f gdvd.f

$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
#	lorder $(LIBOBJ) | tsort | xargs ar cr $(LIBA) #
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)
#	rm -f *.o

TAR = `date '+%y_%m_%d:%H'`_grasp_dvdson.tar   # format: year_month_date_hour
archive :
	tar -cvf $(TAR) *.f Makefile

.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core 

