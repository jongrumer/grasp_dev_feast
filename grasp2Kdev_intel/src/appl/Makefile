SUBDIR = erwf gencsl geniso  \
         jjgen mrgcsl \
         cndens \
         xcsl \
         rscf jj2lsj \
         hfszeeman \
         biotra bioscl jjreduce_breit \
         mcp rci rhfs rhfs_lsj ris \
         mcp3_mpi rscf2_mpi rci_mpi biotra3_mpi bioscl3_mpi zerofirstorder rhfs3_mpi rdensity

TARGETS	 = install

$(TARGETS):
	@for i in $(SUBDIR) ;\
	do \
		echo "Current directory: $$i" ; \
	                cd $$i ; \
			make; \
echo "....Leaving: $$i" ; echo; echo;\
		cd .. ; \
	done

clean :
	@for i in $(SUBDIR) ;\
	do \
	echo "Current directory: $$i" ; \
           cd $$i ; \
           make clean; \
           cd .. ; \
        done

