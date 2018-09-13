LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libdvdson_mpi.a
LIBOBJ = $(OBJ)

OBJ = dvdson_mpi.o iniest.o gdvd.o

$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
#	lorder $(LIBOBJ) | tsort | xargs ar cr $(LIBA) #
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)
	rm -f *.o

.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -o $@

clean:
	-rm -f *.o core 

