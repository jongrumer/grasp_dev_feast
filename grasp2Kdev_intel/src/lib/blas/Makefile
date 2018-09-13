LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/libblas.a
LIBOBJ = $(OBJ) 

OBJ = \
   dasum.o    daxpy.o    dcopy.o    ddot.o     dgbmv.o    dgemm.o    dgemv.o \
   dger.o     dmach.o    dnrm2.o    drot.o     drotg.o    dsbmv.o    dscal.o \
   dspmv.o    dspr.o     dspr2.o    dswap.o    dsymm.o    dsymv.o    dsyr.o \
   dsyr2.o    dsyr2k.o   dsyrk.o    dtbmv.o    dtbsv.o    dtpmv.o    dtpsv.o \
   dtrmm.o    dtrmv.o    dtrsm.o    dtrsv.o    gather.o   idamax.o   iinit.o \
   scatter.o  iinit.o    dgthr.o    droti.o    dsctr.o

$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
#	lorder $(LIBOBJ) | tsort | xargs ar cr $(LIBA) #
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)
	rm -f *.o

TAR = `date '+%y_%m_%d:%H'`_grasp_blas.tar   # format: year_month_date_hour
archive :
	 tar -cvf $(TAR) *.f Makefile
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

clean:
	-rm -f *.o core 

