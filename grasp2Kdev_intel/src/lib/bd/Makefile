LIBDIR = $(GRASP)/lib
LIBO = $(LIBDIR)/graspBD.o
LIBOBJ = $(OBJ)

$(LIBO) : BD.f
	$(FC) -c $(FC_FLAGS) BD.f
	mv -f BD.o $(LIBO)
	@echo The object $(LBF) is now up-to-date
	#cp lapack_d2s.cld $(LIBDIR)
	#cp lapack_d2s.map $(LIBDIR)
clean:
	-rm -f *.o core #${GRASP}/lib/graspBD.o

.f.o:
	$(FC) -c $(FC_FLAGS) $<  -o $@
