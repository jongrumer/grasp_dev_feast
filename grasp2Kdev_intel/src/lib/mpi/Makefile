LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libmpiu.a
# MPI_INC = ${MPI_INCLUDES}


LIBOBJ =  \
   cpath_$(FC_MALLOC).o cslhmpi.o iniestmpi.o lodcslmpi.o \
   lodrwfmpi.o mpiu.o setisompi.o setrwfmpi.o \
   spicmvmpi.o sys_chdir_$(FC_MALLOC).o sys_getwd_$(FC_MALLOC).o  \
   sys_mkdir_$(FC_MALLOC).o

$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
#	lorder $(LIBOBJ) | tsort | xargs ar cr $(LIBA)
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)
	rm -f *.o

.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

