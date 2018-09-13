LIBDIR = $(GRASP)/lib
LIBA = $(LIBDIR)/librang.a
LIBOBJ = $(OBJ) 
OBJ = \
   a1jj.o awp1.o awp1jjg.o c0t5s.o c1e0sm.o c1e1sm.o cle0sm.o coulomb.o \
   diaga1.o diaga2.o \
   diaga3.o diaga4.o diaga5.o eile.o el1.o \
   el2.o el3.o el31.o el32.o el33.o el4.o el41.o el5.o el51.o \
   el52.o el53.o gg1112.o gg1122.o gg1222.o gg1233.o \
   gg1234.o gg12.o glcons.o Gracah1.o itjj2.o \
   itjj3.o itjj.o itrexg.o ittk.o ixjtik.o izas1.o jfaze.o \
   jthn.o mes.o Nine0.o Nine.o nmtejj.o \
   onescalar.o onescalar1.o onescalar2.o oneparticlejj.o \
   oneparticlejj1.o oneparticlejj2.o \
   perko1.o perko2.o \
   reco.o recoonescalar.o reco2.o rec3.o \
   reco3.o reco4.o recop00.o recop1.o recop2.o \
   ribjj.o Rmeajj11.o Rmeajj9.o Rmeajj.o \
   Rmew1jj.o Rmew3jj.o Rmew5jj.o Rmew7bjj.o Rmew7jj.o rumtjj.o \
   Rwjj.o sixj1.o sixj2.o sixj35.o sixj3.o sixj4.o sixj5.o \
   sixj.o suwjj.o termjj.o rkco_gg.o w1jjg.o wap1.o wap1jjg.o \
   wj1.o ww1.o

$(LIBA) : $(LIBOBJ)
	@echo " Building " $(LIBA)
#	lorder $(LIBOBJ) | tsort | xargs ar cr $(LIBA) #
	ar -curvs $(LIBA) $(LIBOBJ)
	ranlib $(LIBA)
#	rm -f *.o

.f.o:
	$(FC) -c $(FC_FLAGS) -I $(GRASP_INCLUDES) $< -o $@

clean:
	-rm -f *.o core 

