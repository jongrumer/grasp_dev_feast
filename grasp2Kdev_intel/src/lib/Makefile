

SUBDIR = bd dvdrc njgraf lib92 lapackd lapacku blas mpi librang # librang_bio
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
clean:
	@for i in $(SUBDIR) ;\
        do \
                echo "Current directory: $$i" ; \
                        cd $$i ; \
                        $(MAKE) $@ ; \
                cd .. ; \
        done


