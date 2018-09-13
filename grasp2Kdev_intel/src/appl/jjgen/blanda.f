*     last edited Januar 2, 1997
      subroutine Blanda(org,varmax,lock,minJ,maxJ,skal,nmax,low,
     :                 posn,posl,lim,dubbel,first)
      integer fil_1,fil_2
      parameter (fil_1=7, fil_2=8)
      integer org(1:15,0:10),antel(1:15,0:10),start(1:15,0:10),skal,cf
      integer ansats(1:15,0:10,0:1),varupp(1:15,0:10),varned(1:15,0:10)
      integer an10,an20,an21,an30,an31,an32,an40,an41,an42,an43,k
      integer an50,an51,an52,an53,an54,an60,an61,an62,an63,an64,an65
      integer an70,an71,an72,an73,an74,an75,an76,stopp(1:15,0:10)
      integer an80,an81,an82,an83,an84,an85,an86,an87,low(1:15,0:10)
      integer an90,an91,an92,an93,an94,an95,an96,an97,an98
      integer anA0,anA1,anA2,anA3,anA4,anA5,anA6,anA7,anA8,anA9,nmax
      integer plus21,plus31,plus32,plus41,plus42,plus43,plus51,plus52
      integer plus53,plus54,plus61,plus62,plus63,plus64,plus65,plus71
      integer plus72,plus73,plus74,plus75,plus76,plus81,plus82,plus83
      integer plus84,plus85,plus86,plus87,plus91,plus92,plus93,plus94
      integer plus95,plus96,plus97,plus98,plusA1,plusA2,plusA3,plusA4
      integer plusA5,plusA6,plusA7,plusA8,plusA9
      integer varmax,par0,par,resS,resL,i,j,antal,posn(110),posl(110)
      integer lim(15),steg(1:15,0:10),dum,ras1,ras3,elar
      integer rasett,rastre,ras(15,0:10),minJ,maxJ
      logical lock(1:15,0:10),first,finns
      logical dubbel(1:15,0:10),napp
      integer plusba,plusca,plusda,plusea,plusfa
      integer plusb1,plusb2,plusb3,plusb4,plusb5,plusb6,plusb7,plusb8,
     :        plusb9
      integer plusc1,plusc2,plusc3,plusc4,plusc5,plusc6,plusc7,plusc8,
     :        plusc9
      integer plusd1,plusd2,plusd3,plusd4,plusd5,plusd6,plusd7,plusd8,
     :        plusd9
      integer pluse1,pluse2,pluse3,pluse4,pluse5,pluse6,pluse7,pluse8,
     :         pluse9
      integer plusf1,plusf2,plusf3,plusf4,plusf5,plusf6,plusf7,plusf8,
     :         plusf9
      integer anba,anca,anda,anea,anfa
      integer anb0,anb1,anb2,anb3,anb4,anb5,anb6,anb7,anb8,anb9
      integer anc0,anc1,anc2,anc3,anc4,anc5,anc6,anc7,anc8,anc9
      integer and0,and1,and2,and3,and4,and5,and6,and7,and8,and9
      integer ane0,ane1,ane2,ane3,ane4,ane5,ane6,ane7,ane8,ane9
      integer anf0,anf1,anf2,anf3,anf4,anf5,anf6,anf7,anf8,anf9

      cf    = 0
      antal = 0
      par0  = 0
      finns = .FALSE.
      do 1 i=1,nmax
         do 1 j=0,min(10,i-1)
            if (dubbel(i,j)) then
               steg(i,j) = -2
            else
               steg(i,j) = -1
            endif
            antal = antal + org(i,j)
   1        par0  = mod(par0+j*org(i,j),2)
      if (nmax.LT.15) then
         do 2 i=nmax+1,15
            do 2 j=0,min(10,i-1)
   2           steg(i,j) = -1
      endif
*     1s
      call Slug(1,0,varmax,varupp,varned,ansats,org,lock(1,0),
     :                      dubbel,low,start(1,0),stopp(1,0))
      do 100 an10 = start(1,0),stopp(1,0),steg(1,0)
         antel(1,0) = an10
         if (antel(1,0).LE.antal .AND. antel(1,0).GE.lim(1)) then
            ansats(1,0,0) = an10
*     2s
      call Slug(2,0,varmax,varupp,varned,ansats,org,lock(2,0),
     :                      dubbel,low,start(2,0),stopp(2,0))
      do 200 an20 = start(2,0),stopp(2,0),steg(2,0)
         antel(2,0) = an20 + antel(1,0)
         if (antel(2,0).LE.antal) then
            ansats(2,0,0) = an20
*     2p
      call Slug(2,1,varmax,varupp,varned,ansats,org,lock(2,1),
     :                      dubbel,low,start(2,1),stopp(2,1))
      do 210 an21 = start(2,1),stopp(2,1),steg(2,1)
         antel(2,1) = an21 + antel(2,0)
         if (antel(2,1).LE.antal .AND. antel(2,1).GE.lim(2)) then 
            do 211 plus21 = min(an21,4),max(an21-2,0),-1
               ansats(2,1,1) = plus21
               ansats(2,1,0) = an21-plus21
*     3s
      call Slug(3,0,varmax,varupp,varned,ansats,org,lock(3,0),
     :                      dubbel,low,start(3,0),stopp(3,0))
      do 300 an30 = start(3,0),stopp(3,0),steg(3,0)
         antel(3,0) = an30 + antel(2,1)
         if (antel(3,0).LE.antal) then 
            ansats(3,0,0) = an30
*     3p
      call Slug(3,1,varmax,varupp,varned,ansats,org,lock(3,1),
     :                      dubbel,low,start(3,1),stopp(3,1))
      do 310 an31 = start(3,1),stopp(3,1),steg(3,1)
         antel(3,1) = an31 + antel(3,0)
         if (antel(3,1).LE.antal) then 
            do 311 plus31 = min(an31,4),max(an31-2,0),-1
               ansats(3,1,1) = plus31
               ansats(3,1,0) = an31-plus31
*     3d
      call Slug(3,2,varmax,varupp,varned,ansats,org,lock(3,2),
     :                      dubbel,low,start(3,2),stopp(3,2))
      do 320 an32 = start(3,2),stopp(3,2),steg(3,2)
         antel(3,2) = an32 + antel(3,1)
         if (antel(3,2).LE.antal .AND. antel(3,2).GE.lim(3)) then 
            do 321 plus32 = min(an32,6),max(an32-4,0),-1
               ansats(3,2,1) = plus32
               ansats(3,2,0) = an32-plus32
*     4s
      call Slug(4,0,varmax,varupp,varned,ansats,org,lock(4,0),
     :                      dubbel,low,start(4,0),stopp(4,0))
      do 400 an40 = start(4,0),stopp(4,0),steg(4,0)
         antel(4,0) = an40 + antel(3,2)
         if (antel(4,0).LE.antal) then 
            ansats(4,0,0) = an40
*     4p
      call Slug(4,1,varmax,varupp,varned,ansats,org,lock(4,1),
     :                      dubbel,low,start(4,1),stopp(4,1))
      do 410 an41 = start(4,1),stopp(4,1),steg(4,1)
         antel(4,1) = an41 + antel(4,0)
         if (antel(4,1).LE.antal) then 
            do 411 plus41 = min(an41,4),max(an41-2,0),-1
               ansats(4,1,1) = plus41
               ansats(4,1,0) = an41-plus41
*     4d
      call Slug(4,2,varmax,varupp,varned,ansats,org,lock(4,2),
     :                      dubbel,low,start(4,2),stopp(4,2))
      do 420 an42 = start(4,2),stopp(4,2),steg(4,2)
         antel(4,2) = an42 + antel(4,1)
         if (antel(4,2).LE.antal) then 
            do 421 plus42 = min(an42,6),max(an42-4,0),-1
               ansats(4,2,1) = plus42
               ansats(4,2,0) = an42-plus42
*     4f
      call Slug(4,3,varmax,varupp,varned,ansats,org,lock(4,3),
     :                      dubbel,low,start(4,3),stopp(4,3))
      do 430 an43 = start(4,3),stopp(4,3),steg(4,3)
         antel(4,3) = an43 + antel(4,2)
         if (antel(4,3).LE.antal .AND. antel(4,3).GE.lim(4)) then 
            do 431 plus43 = min(an43,8),max(an43-6,0),-1
               ansats(4,3,1) = plus43
               ansats(4,3,0) = an43-plus43
*     5s
      call Slug(5,0,varmax,varupp,varned,ansats,org,lock(5,0),
     :                      dubbel,low,start(5,0),stopp(5,0))
      do 500 an50 = start(5,0),stopp(5,0),steg(5,0)
         antel(5,0) = an50 + antel(4,3)
         if (antel(5,0).LE.antal) then 
            ansats(5,0,0) = an50
*     5p
      call Slug(5,1,varmax,varupp,varned,ansats,org,lock(5,1),
     :                      dubbel,low,start(5,1),stopp(5,1))
      do 510 an51 = start(5,1),stopp(5,1),steg(5,1)
         antel(5,1) = an51 + antel(5,0)
         if (antel(5,1).LE.antal) then 
            do 511 plus51 = min(an51,4),max(an51-2,0),-1
               ansats(5,1,1) = plus51
               ansats(5,1,0) = an51-plus51
*     5d
      call Slug(5,2,varmax,varupp,varned,ansats,org,lock(5,2),
     :                      dubbel,low,start(5,2),stopp(5,2))
      do 520 an52 = start(5,2),stopp(5,2),steg(5,2)
         antel(5,2) = an52 + antel(5,1)
         if (antel(5,2).LE.antal) then 
            do 521 plus52 = min(an52,6),max(an52-4,0),-1
               ansats(5,2,1) = plus52
               ansats(5,2,0) = an52-plus52

*     5f
      call Slug(5,3,varmax,varupp,varned,ansats,org,lock(5,3),
     :                      dubbel,low,start(5,3),stopp(5,3))
      do 530 an53 = start(5,3),stopp(5,3),steg(5,3)
         antel(5,3) = an53 + antel(5,2)
         if (antel(5,3).LE.antal) then 
            do 531 plus53 = min(an53,8),max(an53-6,0),-1
               ansats(5,3,1) = plus53
               ansats(5,3,0) = an53-plus53
*     5g
      call Slug(5,4,varmax,varupp,varned,ansats,org,lock(5,4),
     :                      dubbel,low,start(5,4),stopp(5,4))
      do 540 an54 = start(5,4),stopp(5,4),steg(5,4)
         antel(5,4) = an54 + antel(5,3)
         if (antel(5,4).LE.antal .AND. antel(5,4).GE.lim(5)) then 
            do 541 plus54 = min(an54,10),max(an54-8,0),-1
               ansats(5,4,1) = plus54
               ansats(5,4,0) = an54-plus54
*     6s
      call Slug(6,0,varmax,varupp,varned,ansats,org,lock(6,0),
     :                      dubbel,low,start(6,0),stopp(6,0))
      do 600 an60 = start(6,0),stopp(6,0),steg(6,0)
         antel(6,0) = an60 + antel(5,4)
         if (antel(6,0).LE.antal .AND. ansats(5,4,1).LE.2) then 
            ansats(6,0,0) = an60
*     6p
      call Slug(6,1,varmax,varupp,varned,ansats,org,lock(6,1),
     :                      dubbel,low,start(6,1),stopp(6,1))
      do 610 an61 = start(6,1),stopp(6,1),steg(6,1)
         antel(6,1) = an61 + antel(6,0)
         if (antel(6,1).LE.antal) then 
            do 611 plus61 = min(an61,4),max(an61-2,0),-1
               ansats(6,1,1) = plus61
               ansats(6,1,0) = an61-plus61
*     6d
      call Slug(6,2,varmax,varupp,varned,ansats,org,lock(6,2),
     :                      dubbel,low,start(6,2),stopp(6,2))
      do 620 an62 = start(6,2),stopp(6,2),steg(6,2)
         antel(6,2) = an62 + antel(6,1)
         if (antel(6,2).LE.antal) then 
            do 621 plus62 = min(an62,6),max(an62-4,0),-1
               ansats(6,2,1) = plus62
               ansats(6,2,0) = an62-plus62
*     6f
      call Slug(6,3,varmax,varupp,varned,ansats,org,lock(6,3),
     :                      dubbel,low,start(6,3),stopp(6,3))
      do 630 an63 = start(6,3),stopp(6,3),steg(6,3)
         antel(6,3) = an63 + antel(6,2)
         if (antel(6,3).LE.antal) then 
            do 631 plus63 = min(an63,8),max(an63-6,0),-1
               ansats(6,3,1) = plus63
               ansats(6,3,0) = an63-plus63
*     6g
      call Slug(6,4,varmax,varupp,varned,ansats,org,lock(6,4),
     :                      dubbel,low,start(6,4),stopp(6,4))
      do 640 an64 = start(6,4),stopp(6,4),steg(6,4)
         antel(6,4) = an64 + antel(6,3)
         if (antel(6,4).LE.antal) then 
            do 641 plus64 = min(an64,10),max(an64-8,0),-1
               ansats(6,4,1) = plus64
               ansats(6,4,0) = an64-plus64
*     6h
      call Slug(6,5,varmax,varupp,varned,ansats,org,lock(6,5),
     :                      dubbel,low,start(6,5),stopp(6,5))
      do 650 an65 = start(6,5),stopp(6,5),steg(6,5)
         antel(6,5) = an65 + antel(6,4)
         if (antel(6,5).LE.antal .AND. ansats(6,4,1).LE.2 .AND.
     :                                antel(6,5).GE.lim(6)) then 
            do 651 plus65 = min(an65,12),max(an65-10,0),-1
               ansats(6,5,1) = plus65
               ansats(6,5,0) = an65-plus65
*     7s
      call Slug(7,0,varmax,varupp,varned,ansats,org,lock(7,0),
     :                      dubbel,low,start(7,0),stopp(7,0))
      do 700 an70 = start(7,0),stopp(7,0),steg(7,0)
         antel(7,0) = an70 + antel(6,5)
         if (antel(7,0).LE.antal .AND. ansats(6,5,1).LE.2 .AND.
     :                                   ansats(6,5,0).LE.2) then 
            ansats(7,0,0) = an70
*     7p
      call Slug(7,1,varmax,varupp,varned,ansats,org,lock(7,1),
     :                      dubbel,low,start(7,1),stopp(7,1))
      do 710 an71 = start(7,1),stopp(7,1),steg(7,1)
         antel(7,1) = an71 + antel(7,0)
         if (antel(7,1).LE.antal) then 
            do 711 plus71 = min(an71,4),max(an71-2,0),-1
               ansats(7,1,1) = plus71
               ansats(7,1,0) = an71-plus71
*     7d
      call Slug(7,2,varmax,varupp,varned,ansats,org,lock(7,2),
     :                      dubbel,low,start(7,2),stopp(7,2))
      do 720 an72 = start(7,2),stopp(7,2),steg(7,2)
         antel(7,2) = an72 + antel(7,1)
         if (antel(7,2).LE.antal) then 
            do 721 plus72 = min(an72,6),max(an72-4,0),-1
               ansats(7,2,1) = plus72
               ansats(7,2,0) = an72-plus72
*     7f
      call Slug(7,3,varmax,varupp,varned,ansats,org,lock(7,3),
     :                      dubbel,low,start(7,3),stopp(7,3))
      do 730 an73 = start(7,3),stopp(7,3),steg(7,3)
         antel(7,3) = an73 + antel(7,2)
         if (antel(7,3).LE.antal) then 
            do 731 plus73 = min(an73,8),max(an73-6,0),-1
               ansats(7,3,1) = plus73
               ansats(7,3,0) = an73-plus73
*     7g
      call Slug(7,4,varmax,varupp,varned,ansats,org,lock(7,4),
     :                      dubbel,low,start(7,4),stopp(7,4))
      do 740 an74 = start(7,4),stopp(7,4),steg(7,4)
         antel(7,4) = an74 + antel(7,3)
         if (antel(7,4).LE.antal) then 
            do 741 plus74 = min(an74,10),max(an74-8,0),-1
               ansats(7,4,1) = plus74
               ansats(7,4,0) = an74-plus74
*     7h
      call Slug(7,5,varmax,varupp,varned,ansats,org,lock(7,5),
     :                      dubbel,low,start(7,5),stopp(7,5))
      do 750 an75 = start(7,5),stopp(7,5),steg(7,5)
         antel(7,5) = an75 + antel(7,4)
         if (antel(7,5).LE.antal .AND. ansats(7,4,1).LE.2) then 
            do 751 plus75 = min(an75,12),max(an75-10,0),-1
               ansats(7,5,1) = plus75
               ansats(7,5,0) = an75-plus75
*     7i
      call Slug(7,6,varmax,varupp,varned,ansats,org,lock(7,6),
     :                      dubbel,low,start(7,6),stopp(7,6))
      do 760 an76 = start(7,6),stopp(7,6),steg(7,6)
         antel(7,6) = an76 + antel(7,5)
         if (antel(7,6).LE.antal .AND. ansats(7,5,1).LE.2 .AND.
     :            ansats(7,5,0).LE.2 .AND. antel(7,6).GE.lim(7)) then 
            do 761 plus76 = min(an76,14),max(an76-12,0),-1
               ansats(7,6,1) = plus76
               ansats(7,6,0) = an76-plus76
*     8s
      call Slug(8,0,varmax,varupp,varned,ansats,org,lock(8,0),
     :                      dubbel,low,start(8,0),stopp(8,0))
      do 800 an80 = start(8,0),stopp(8,0),steg(8,0)
         antel(8,0) = an80 + antel(7,6)
         if (antel(8,0).LE.antal .AND. ansats(7,6,1).LE.2 .AND.
     :                                   ansats(7,6,0).LE.2) then 
            ansats(8,0,0) = an80
*     8p
      call Slug(8,1,varmax,varupp,varned,ansats,org,lock(8,1),
     :                      dubbel,low,start(8,1),stopp(8,1))
      do 810 an81 = start(8,1),stopp(8,1),steg(8,1)
         antel(8,1) = an81 + antel(8,0)
         if (antel(8,1).LE.antal) then 
            do 811 plus81 = min(an81,4),max(an81-2,0),-1
               ansats(8,1,1) = plus81
               ansats(8,1,0) = an81-plus81
*     8d
      call Slug(8,2,varmax,varupp,varned,ansats,org,lock(8,2),
     :                      dubbel,low,start(8,2),stopp(8,2))
      do 820 an82 = start(8,2),stopp(8,2),steg(8,2)
         antel(8,2) = an82 + antel(8,1)
         if (antel(8,2).LE.antal) then 
            do 821 plus82 = min(an82,6),max(an82-4,0),-1
               ansats(8,2,1) = plus82
               ansats(8,2,0) = an82-plus82
*     8f
      call Slug(8,3,varmax,varupp,varned,ansats,org,lock(8,3),
     :                      dubbel,low,start(8,3),stopp(8,3))
      do 830 an83 = start(8,3),stopp(8,3),steg(8,3)
         antel(8,3) = an83 + antel(8,2)
         if (antel(8,3).LE.antal) then 
            do 831 plus83 = min(an83,8),max(an83-6,0),-1
               ansats(8,3,1) = plus83
               ansats(8,3,0) = an83-plus83
*     8g
      call Slug(8,4,varmax,varupp,varned,ansats,org,lock(8,4),
     :                      dubbel,low,start(8,4),stopp(8,4))

      do 840 an84 = start(8,4),stopp(8,4),steg(8,4)
         antel(8,4) = an84 + antel(8,3)
         if (antel(8,4).LE.antal) then 
            do 841 plus84 = min(an84,10),max(an84-8,0),-1
               ansats(8,4,1) = plus84
               ansats(8,4,0) = an84-plus84
*     8h
      call Slug(8,5,varmax,varupp,varned,ansats,org,lock(8,5),
     :                      dubbel,low,start(8,5),stopp(8,5))
      do 850 an85 = start(8,5),stopp(8,5),steg(8,5)
         antel(8,5) = an85 + antel(8,4)
         if (antel(8,5).LE.antal .AND. ansats(8,4,1).LE.2) then 
            do 851 plus85 = min(an85,12),max(an85-10,0),-1
               ansats(8,5,1) = plus85
               ansats(8,5,0) = an85-plus85
*     8i
      call Slug(8,6,varmax,varupp,varned,ansats,org,lock(8,6),
     :                      dubbel,low,start(8,6),stopp(8,6))
      do 860 an86 = start(8,6),stopp(8,6),steg(8,6)
         antel(8,6) = an86 + antel(8,5)
         if (antel(8,6).LE.antal .AND. ansats(8,5,1).LE.2 .AND.
     :                                    ansats(8,5,0).LE.2) then 
            do 861 plus86 = min(an86,14),max(an86-12,0),-1
               ansats(8,6,1) = plus86
               ansats(8,6,0) = an86-plus86
*     8k
      call Slug(8,7,varmax,varupp,varned,ansats,org,lock(8,7),
     :                      dubbel,low,start(8,7),stopp(8,7))
      do 870 an87 = start(8,7),stopp(8,7),steg(8,7)
         antel(8,7) = an87 + antel(8,6)
         if (antel(8,7).LE.antal .AND. ansats(8,6,1).LE.2 .AND.
     :            ansats(8,6,0).LE.2 .AND. antel(8,7).GE.lim(8)) then 
            do 871 plus87 = min(an87,16),max(an87-14,0),-1
               ansats(8,7,1) = plus87
               ansats(8,7,0) = an87-plus87
*     9s
      call Slug(9,0,varmax,varupp,varned,ansats,org,lock(9,0),
     :                      dubbel,low,start(9,0),stopp(9,0))
      do 900 an90 = start(9,0),stopp(9,0),steg(9,0)
         antel(9,0) = an90 + antel(8,7)
         if (antel(9,0).LE.antal .AND. ansats(8,7,1).LE.2 .AND.
     :                                      ansats(8,7,0).LE.2) then 
            ansats(9,0,0) = an90
*     9p
      call Slug(9,1,varmax,varupp,varned,ansats,org,lock(9,1),
     :                      dubbel,low,start(9,1),stopp(9,1))
      do 910 an91 = start(9,1),stopp(9,1),steg(9,1)
         antel(9,1) = an91 + antel(9,0)
         if (antel(9,1).LE.antal) then 
            do 911 plus91 = min(an91,4),max(an91-2,0),-1
               ansats(9,1,1) = plus91
               ansats(9,1,0) = an91-plus91
*     9d
      call Slug(9,2,varmax,varupp,varned,ansats,org,lock(9,2),
     :                      dubbel,low,start(9,2),stopp(9,2))
      do 920 an92 = start(9,2),stopp(9,2),steg(9,2)
         antel(9,2) = an92 + antel(9,1)
         if (antel(9,2).LE.antal) then 
            do 921 plus92 = min(an92,6),max(an92-4,0),-1
               ansats(9,2,1) = plus92
               ansats(9,2,0) = an92-plus92
*     9f
      call Slug(9,3,varmax,varupp,varned,ansats,org,lock(9,3),
     :                      dubbel,low,start(9,3),stopp(9,3))
      do 930 an93 = start(9,3),stopp(9,3),steg(9,3)
         antel(9,3) = an93 + antel(9,2)
         if (antel(9,3).LE.antal) then 
            do 931 plus93 = min(an93,8),max(an93-6,0),-1
               ansats(9,3,1) = plus93
               ansats(9,3,0) = an93-plus93
*     9g
      call Slug(9,4,varmax,varupp,varned,ansats,org,lock(9,4),
     :                      dubbel,low,start(9,4),stopp(9,4))
      do 940 an94 = start(9,4),stopp(9,4),steg(9,4)
         antel(9,4) = an94 + antel(9,3)
         if (antel(9,4).LE.antal) then 
            do 941 plus94 = min(an94,10),max(an94-8,0),-1
               ansats(9,4,1) = plus94
               ansats(9,4,0) = an94-plus94
*     9h
      call Slug(9,5,varmax,varupp,varned,ansats,org,lock(9,5),
     :                      dubbel,low,start(9,5),stopp(9,5))
      do 950 an95 = start(9,5),stopp(9,5),steg(9,5)
         antel(9,5) = an95 + antel(9,4)
         if (antel(9,5).LE.antal .AND. ansats(9,4,1).LE.2) then 
            do 951 plus95 = min(an95,12),max(an95-10,0),-1
               ansats(9,5,1) = plus95
               ansats(9,5,0) = an95-plus95
*     9i
      call Slug(9,6,varmax,varupp,varned,ansats,org,lock(9,6),
     :                      dubbel,low,start(9,6),stopp(9,6))
      do 960 an96 = start(9,6),stopp(9,6),steg(9,6)
         antel(9,6) = an96 + antel(9,5)
         if (antel(9,6).LE.antal .AND. ansats(9,5,1).LE.2 .AND.
     :                                    ansats(9,5,0).LE.2) then 
            do 961 plus96 = min(an96,14),max(an96-12,0),-1
               ansats(9,6,1) = plus96
               ansats(9,6,0) = an96-plus96
*     9k
      call Slug(9,7,varmax,varupp,varned,ansats,org,lock(9,7),
     :                      dubbel,low,start(9,7),stopp(9,7))
      do 970 an97 = start(9,7),stopp(9,7),steg(9,7)
         antel(9,7) = an97 + antel(9,6)
         if (antel(9,7).LE.antal .AND. ansats(9,6,1).LE.2 .AND.
     :                                    ansats(9,6,0).LE.2) then 
            do 971 plus97 = min(an97,16),max(an97-14,0),-1
               ansats(9,7,1) = plus97
               ansats(9,7,0) = an97-plus97
*     9l
      call Slug(9,8,varmax,varupp,varned,ansats,org,lock(9,8),
     :                      dubbel,low,start(9,8),stopp(9,8))
      do 980 an98 = start(9,8),stopp(9,8),steg(9,8)
         antel(9,8) = an98 + antel(9,7)
         if (antel(9,8).LE.antal .AND. ansats(9,7,1).LE.2 .AND.
     :            ansats(9,7,0).LE.2 .AND. antel(9,8).GE.lim(9)) then 
            do 981 plus98 = min(an98,18),max(an98-16,0),-1
               ansats(9,8,1) = plus98
               ansats(9,8,0) = an98-plus98
*     10s
      call Slug(10,0,varmax,varupp,varned,ansats,org,lock(10,0),
     :                      dubbel,low,start(10,0),stopp(10,0))
      do 1000 anA0 = start(10,0),stopp(10,0),steg(10,0)
         antel(10,0) = anA0 + antel(9,8)
         if (antel(10,0).LE.antal .AND. ansats(9,8,1).LE.2 .AND.
     :                                   ansats(9,8,0).LE.2) then 
            ansats(10,0,0) = anA0
*     10p
      call Slug(10,1,varmax,varupp,varned,ansats,org,lock(10,1),
     :                      dubbel,low,start(10,1),stopp(10,1))
      do 1010 anA1 = start(10,1),stopp(10,1),steg(10,1)
         antel(10,1) = anA1 + antel(10,0)
         if (antel(10,1).LE.antal) then
            do 1011 plusA1 = min(anA1,4),max(anA1-2,0),-1
               ansats(10,1,1) = plusA1
               ansats(10,1,0) = anA1-plusA1
*     10d
      call Slug(10,2,varmax,varupp,varned,ansats,org,lock(10,2),
     :                      dubbel,low,start(10,2),stopp(10,2))
      do 1020 anA2 = start(10,2),stopp(10,2),steg(10,2)
         antel(10,2) = anA2 + antel(10,1)
         if (antel(10,2).LE.antal) then 
            do 1021 plusA2 = min(anA2,6),max(anA2-4,0),-1
               ansats(10,2,1) = plusA2
               ansats(10,2,0) = anA2-plusA2
*     10f
      call Slug(10,3,varmax,varupp,varned,ansats,org,lock(10,3),
     :                      dubbel,low,start(10,3),stopp(10,3))
      do 1030 anA3 = start(10,3),stopp(10,3),steg(10,3)
         antel(10,3) = anA3 + antel(10,2)
         if (antel(10,3).LE.antal) then 
            do 1031 plusA3 = min(anA3,8),max(anA3-6,0),-1
               ansats(10,3,1) = plusA3
               ansats(10,3,0) = anA3-plusA3
*     10g
      call Slug(10,4,varmax,varupp,varned,ansats,org,lock(10,4),
     :                      dubbel,low,start(10,4),stopp(10,4))
      do 1040 anA4 = start(10,4),stopp(10,4),steg(10,4)
         antel(10,4) = anA4 + antel(10,3)
         if (antel(10,4).LE.antal) then 
            do 1041 plusA4 = min(anA4,10),max(anA4-8,0),-1
               ansats(10,4,1) = plusA4
               ansats(10,4,0) = anA4-plusA4
*     10h
      call Slug(10,5,varmax,varupp,varned,ansats,org,lock(10,5),
     :                      dubbel,low,start(10,5),stopp(10,5))
      do 1050 anA5 = start(10,5),stopp(10,5),steg(10,5)
         antel(10,5) = anA5 + antel(10,4)
         if (antel(10,5).LE.antal .AND. ansats(10,4,1).LE.2) then 
            do 1051 plusA5 = min(anA5,12),max(anA5-10,0),-1
               ansats(10,5,1) = plusA5
               ansats(10,5,0) = anA5-plusA5
*     10i
      call Slug(10,6,varmax,varupp,varned,ansats,org,lock(10,6),
     :                      dubbel,low,start(10,6),stopp(10,6))
      do 1060 anA6 = start(10,6),stopp(10,6),steg(10,6)
         antel(10,6) = anA6 + antel(10,5)
         if (antel(10,6).LE.antal .AND. ansats(10,5,1).LE.2 .AND.
     :                                        ansats(10,5,0).LE.2) then 
            do 1061 plusA6 = min(anA6,14),max(anA6-12,0),-1
               ansats(10,6,1) = plusA6
               ansats(10,6,0) = anA6-plusA6
*     10k
      call Slug(10,7,varmax,varupp,varned,ansats,org,lock(10,7),
     :                      dubbel,low,start(10,7),stopp(10,7))
      do 1070 anA7 = start(10,7),stopp(10,7),steg(10,7)
         antel(10,7) = anA7 + antel(10,6)
         if (antel(10,7).LE.antal .AND. ansats(10,6,1).LE.2 .AND.
     :                                         ansats(10,6,0).LE.2) then 
            do 1071 plusA7 = min(anA7,16),max(anA7-14,0),-1
               ansats(10,7,1) = plusA7
               ansats(10,7,0) = anA7-plusA7
*     10l
      call Slug(10,8,varmax,varupp,varned,ansats,org,lock(10,8),
     :                      dubbel,low,start(10,8),stopp(10,8))
      do 1080 anA8 = start(10,8),stopp(10,8),steg(10,8)
         antel(10,8) = anA8 + antel(10,7)
         if (antel(10,8).LE.antal .AND. ansats(10,7,1).LE.2 .AND.
     :                                         ansats(10,7,0).LE.2) then 
            do 1081 plusA8 = min(anA8,18),max(anA8-16,0),-1
               ansats(10,8,1) = plusA8
               ansats(10,8,0) = anA8-plusA8
*     10m
      call Slug(10,9,varmax,varupp,varned,ansats,org,lock(10,9),
     :                      dubbel,low,start(10,9),stopp(10,9))
      do 1090 anA9 = start(10,9),stopp(10,9),steg(10,9)
         antel(10,9) = anA9 + antel(10,8)
         if (antel(10,9).LE.antal .AND. ansats(10,8,1).LE.2 .AND.
     :            ansats(10,8,0).LE.2 .AND. antel(10,9).GE.lim(10)) then 
            do 1091 plusA9 = min(anA9,20),max(anA9-18,0),-1
               ansats(10,9,1) = plusA9
               ansats(10,9,0) = anA9-plusA9  
*     11s
      call Slug(11,0,varmax,varupp,varned,ansats,org,lock(11,0),
     :                      dubbel,low,start(11,0),stopp(11,0))
      do 1100 anB0 = start(11,0),stopp(11,0),steg(11,0)
         antel(11,0) = anB0 + antel(10,9)
         if (antel(11,0).LE.antal .AND. ansats(10,9,1).LE.2 .AND.
     :                                   ansats(10,9,0).LE.2) then 
            ansats(11,0,0) = anB0
*     11p
      call Slug(11,1,varmax,varupp,varned,ansats,org,lock(11,1),
     :                      dubbel,low,start(11,1),stopp(11,1))
      do 1110 anB1 = start(11,1),stopp(11,1),steg(11,1)
         antel(11,1) = anB1 + antel(11,0)
         if (antel(11,1).LE.antal) then 
            do 1111 plusB1 = min(anB1,4),max(anB1-2,0),-1
               ansats(11,1,1) = plusB1
               ansats(11,1,0) = anB1-plusB1
*     11d
      call Slug(11,2,varmax,varupp,varned,ansats,org,lock(11,2),
     :                      dubbel,low,start(11,2),stopp(11,2))
      do 1120 anB2 = start(11,2),stopp(11,2),steg(11,2)
         antel(11,2) = anB2 + antel(11,1)
         if (antel(11,2).LE.antal) then 
            do 1121 plusB2 = min(anB2,6),max(anB2-4,0),-1
               ansats(11,2,1) = plusB2
               ansats(11,2,0) = anB2-plusB2
*     11f
      call Slug(11,3,varmax,varupp,varned,ansats,org,lock(11,3),
     :                      dubbel,low,start(11,3),stopp(11,3))
      do 1130 anB3 = start(11,3),stopp(11,3),steg(11,3)
         antel(11,3) = anB3 + antel(11,2)
         if (antel(11,3).LE.antal) then 
            do 1131 plusB3 = min(anB3,8),max(anB3-6,0),-1
               ansats(11,3,1) = plusB3
               ansats(11,3,0) = anB3-plusB3
*     11g
      call Slug(11,4,varmax,varupp,varned,ansats,org,lock(11,4),
     :                      dubbel,low,start(11,4),stopp(11,4))
      do 1140 anB4 = start(11,4),stopp(11,4),steg(11,4)
         antel(11,4) = anB4 + antel(11,3)
         if (antel(11,4).LE.antal) then 
            do 1141 plusB4 = min(anB4,10),max(anB4-8,0),-1
               ansats(11,4,1) = plusB4
               ansats(11,4,0) = anB4-plusB4
*     11h
      call Slug(11,5,varmax,varupp,varned,ansats,org,lock(11,5),
     :                      dubbel,low,start(11,5),stopp(11,5))
      do 1150 anB5 = start(11,5),stopp(11,5),steg(11,5)
         antel(11,5) = anB5 + antel(11,4)
         if (antel(11,5).LE.antal .AND. ansats(11,4,1).LE.2) then 
            do 1151 plusB5 = min(anB5,12),max(anB5-10,0),-1
               ansats(11,5,1) = plusB5
               ansats(11,5,0) = anB5-plusB5
*     11i
      call Slug(11,6,varmax,varupp,varned,ansats,org,lock(11,6),
     :                      dubbel,low,start(11,6),stopp(11,6))
      do 1160 anB6 = start(11,6),stopp(11,6),steg(11,6)
         antel(11,6) = anB6 + antel(11,5)
         if (antel(11,6).LE.antal .AND. ansats(11,5,1).LE.2 .AND.
     :                                    ansats(11,5,0).LE.2) then         
            do 1161 plusB6 = min(anB6,14),max(anB6-12,0),-1
               ansats(11,6,1) = plusB6
               ansats(11,6,0) = anB6-plusB6
*     11k
      call Slug(11,7,varmax,varupp,varned,ansats,org,lock(11,7),
     :                      dubbel,low,start(11,7),stopp(11,7))
      do 1170 anB7 = start(11,7),stopp(11,7),steg(11,7)
         antel(11,7) = anB7 + antel(11,6)
         if (antel(11,7).LE.antal .AND. ansats(11,6,1).LE.2 .AND.
     :                                    ansats(11,6,0).LE.2) then 
            do 1171 plusB7 = min(anB7,16),max(anB7-14,0),-1
               ansats(11,7,1) = plusB7
               ansats(11,7,0) = anB7-plusB7
*     11l
      call Slug(11,8,varmax,varupp,varned,ansats,org,lock(11,8),
     :                      dubbel,low,start(11,8),stopp(11,8))
      do 1180 anB8 = start(11,8),stopp(11,8),steg(11,8)
         antel(11,8) = anB8 + antel(11,7)
         if (antel(11,8).LE.antal .AND. ansats(11,7,1).LE.2 .AND.
     :                                    ansats(11,7,0).LE.2) then 
            do 1181 plusB8 = min(anB8,18),max(anB8-16,0),-1
               ansats(11,8,1) = plusB8
               ansats(11,8,0) = anB8-plusB8
*     11m
      call Slug(11,9,varmax,varupp,varned,ansats,org,lock(11,9),
     :                      dubbel,low,start(11,9),stopp(11,9))
      do 1190 anB9 = start(11,9),stopp(11,9),steg(11,9)
         antel(11,9) = anB9 + antel(11,8)
         if (antel(11,9).LE.antal .AND. ansats(11,8,1).LE.2 .AND.
     :                                       ansats(11,8,0).LE.2) then 
            do 1191 plusB9 = min(anB9,20),max(anB9-18,0),-1
               ansats(11,9,1) = plusB9
               ansats(11,9,0) = anB9-plusB9
*     11n
      call Slug(11,10,varmax,varupp,varned,ansats,org,lock(11,10),
     :                      dubbel,low,start(11,10),stopp(11,10))
      do 11100 anBA = start(11,10),stopp(11,10),steg(11,10)
         antel(11,10) = anBA + antel(11,9)
         if (antel(11,10).LE.antal .AND. ansats(11,9,1).LE.2 .AND.
     :         ansats(11,9,0).LE.2 .AND. antel(11,10).GE.lim(11)) then 
            do 11101 plusBA = min(anBA,22),max(anBA-20,0),-1
               ansats(11,10,1) = plusBA
               ansats(11,10,0) = anBA-plusBA
*     12s
      call Slug(12,0,varmax,varupp,varned,ansats,org,lock(12,0),
     :                      dubbel,low,start(12,0),stopp(12,0))
      do 1200 anC0 = start(12,0),stopp(12,0),steg(12,0)
         antel(12,0) = anC0 + antel(11,10)
         if (antel(12,0).LE.antal .AND. ansats(11,10,1).LE.2 .AND.
     :                                   ansats(11,10,0).LE.2) then 
            ansats(12,0,0) = anC0
*     12p
      call Slug(12,1,varmax,varupp,varned,ansats,org,lock(12,1),
     :                      dubbel,low,start(12,1),stopp(12,1))
      do 1210 anC1 = start(12,1),stopp(12,1),steg(12,1)
         antel(12,1) = anC1 + antel(12,0)
         if (antel(12,1).LE.antal) then 
            do 1211 plusC1 = min(anC1,4),max(anC1-2,0),-1
               ansats(12,1,1) = plusC1
               ansats(12,1,0) = anC1-plusC1
*     12d
      call Slug(12,2,varmax,varupp,varned,ansats,org,lock(12,2),
     :                      dubbel,low,start(12,2),stopp(12,2))
      do 1220 anC2 = start(12,2),stopp(12,2),steg(12,2)
         antel(12,2) = anC2 + antel(12,1)
         if (antel(12,2).LE.antal) then 
            do 1221 plusC2 = min(anC2,6),max(anC2-4,0),-1
               ansats(12,2,1) = plusC2
               ansats(12,2,0) = anC2-plusC2
*     12f
      call Slug(12,3,varmax,varupp,varned,ansats,org,lock(12,3),
     :                      dubbel,low,start(12,3),stopp(12,3))
      do 1230 anC3 = start(12,3),stopp(12,3),steg(12,3)
         antel(12,3) = anC3 + antel(12,2)
         if (antel(12,3).LE.antal) then 
            do 1231 plusC3 = min(anC3,8),max(anC3-6,0),-1
               ansats(12,3,1) = plusC3
               ansats(12,3,0) = anC3-plusC3
*     12g
      call Slug(12,4,varmax,varupp,varned,ansats,org,lock(12,4),
     :                      dubbel,low,start(12,4),stopp(12,4))
      do 1240 anC4 = start(12,4),stopp(12,4),steg(12,4)
         antel(12,4) = anC4 + antel(12,3)
         if (antel(12,4).LE.antal) then 
            do 1241 plusC4 = min(anC4,10),max(anC4-8,0),-1
               ansats(12,4,1) = plusC4
               ansats(12,4,0) = anC4-plusC4
*     12h
      call Slug(12,5,varmax,varupp,varned,ansats,org,lock(12,5),
     :                      dubbel,low,start(12,5),stopp(12,5))
      do 1250 anC5 = start(12,5),stopp(12,5),steg(12,5)
         antel(12,5) = anC5 + antel(12,4)
         if (antel(12,5).LE.antal .AND. ansats(12,4,1).LE.2) then 
            do 1251 plusC5 = min(anC5,12),max(anC5-10,0),-1
               ansats(12,5,1) = plusC5
               ansats(12,5,0) = anC5-plusC5
*     12i
      call Slug(12,6,varmax,varupp,varned,ansats,org,lock(12,6),
     :                      dubbel,low,start(12,6),stopp(12,6))
      do 1260 anC6 = start(12,6),stopp(12,6),steg(12,6)
         antel(12,6) = anC6 + antel(12,5)
         if (antel(12,6).LE.antal .AND. ansats(12,5,1).LE.2 .AND.
     :                                       ansats(12,5,0).LE.2) then 
            do 1261 plusC6 = min(anC6,14),max(anC6-12,0),-1
               ansats(12,6,1) = plusC6
               ansats(12,6,0) = anC6-plusC6
*     12k
      call Slug(12,7,varmax,varupp,varned,ansats,org,lock(12,7),
     :                      dubbel,low,start(12,7),stopp(12,7))
      do 1270 anC7 = start(12,7),stopp(12,7),steg(12,7)
         antel(12,7) = anC7 + antel(12,6)
         if (antel(12,7).LE.antal .AND. ansats(12,6,1).LE.2 .AND.
     :                                         ansats(12,6,0).LE.2) then 
            do 1271 plusC7 = min(anC7,16),max(anC7-14,0),-1
               ansats(12,7,1) = plusC7
               ansats(12,7,0) = anC7-plusC7
*     12l
      call Slug(12,8,varmax,varupp,varned,ansats,org,lock(12,8),
     :                      dubbel,low,start(12,8),stopp(12,8))
      do 1280 anC8 = start(12,8),stopp(12,8),steg(12,8)
         antel(12,8) = anC8 + antel(12,7)
         if (antel(12,8).LE.antal .AND. ansats(12,7,1).LE.2 .AND.
     :                                         ansats(12,7,0).LE.2) then 
            do 1281 plusC8 = min(anC8,18),max(anC8-16,0),-1
               ansats(12,8,1) = plusC8
               ansats(12,8,0) = anC8-plusC8
*     12m
      call Slug(12,9,varmax,varupp,varned,ansats,org,lock(12,9),
     :                      dubbel,low,start(12,9),stopp(12,9))
      do 1290 anC9 = start(12,9),stopp(12,9),steg(12,9)
         antel(12,9) = anC9 + antel(12,8)
         if (antel(12,9).LE.antal .AND. ansats(12,8,1).LE.2 .AND.
     :                                         ansats(12,8,0).LE.2) then 
            do 1291 plusC9 = min(anC9,20),max(anC9-18,0),-1
               ansats(12,9,1) = plusC9
               ansats(12,9,0) = anC9-plusC9
*     12n
      call Slug(12,10,varmax,varupp,varned,ansats,org,lock(12,10),
     :                      dubbel,low,start(12,10),stopp(12,10))
      do 12100 anCA = start(12,10),stopp(12,10),steg(12,10)
         antel(12,10) = anCA + antel(12,9)
         if (antel(12,10).LE.antal .AND. ansats(12,9,1).LE.2 .AND.
     :           ansats(12,9,0).LE.2 .AND. antel(12,10).GE.lim(12)) then 
            do 12101 plusCA = min(anCA,22),max(anCA-20,0),-1
               ansats(12,10,1) = plusCA
               ansats(12,10,0) = anCA-plusCA
*     13s
      call Slug(13,0,varmax,varupp,varned,ansats,org,lock(13,0),
     :                      dubbel,low,start(13,0),stopp(13,0))
      do 1300 anD0 = start(13,0),stopp(13,0),steg(13,0)
         antel(13,0) = anD0 + antel(12,10)
         if (antel(13,0).LE.antal .AND. ansats(12,10,1).LE.2 .AND.
     :                                   ansats(12,10,0).LE.2) then 
            ansats(13,0,0) = anD0
*     13p
      call Slug(13,1,varmax,varupp,varned,ansats,org,lock(13,1),
     :                      dubbel,low,start(13,1),stopp(13,1))
      do 1310 anD1 = start(13,1),stopp(13,1),steg(13,1)
         antel(13,1) = anD1 + antel(13,0)
         if (antel(13,1).LE.antal) then 
            do 1311 plusD1 = min(anD1,4),max(anD1-2,0),-1
               ansats(13,1,1) = plusD1
               ansats(13,1,0) = anD1-plusD1
*     13d
      call Slug(13,2,varmax,varupp,varned,ansats,org,lock(13,2),
     :                      dubbel,low,start(13,2),stopp(13,2))
      do 1320 anD2 = start(13,2),stopp(13,2),steg(13,2)
         antel(13,2) = anD2 + antel(13,1)
         if (antel(13,2).LE.antal) then 
            do 1321 plusD2 = min(anD2,6),max(anD2-4,0),-1
               ansats(13,2,1) = plusD2
               ansats(13,2,0) = anD2-plusD2
*     13f
      call Slug(13,3,varmax,varupp,varned,ansats,org,lock(13,3),
     :                      dubbel,low,start(13,3),stopp(13,3))
      do 1330 anD3 = start(13,3),stopp(13,3),steg(13,3)
         antel(13,3) = anD3 + antel(13,2)
         if (antel(13,3).LE.antal) then 
            do 1331 plusD3 = min(anD3,8),max(anD3-6,0),-1
               ansats(13,3,1) = plusD3
               ansats(13,3,0) = anD3-plusD3
*     13g
      call Slug(13,4,varmax,varupp,varned,ansats,org,lock(13,4),
     :                      dubbel,low,start(13,4),stopp(13,4))
      do 1340 anD4 = start(13,4),stopp(13,4),steg(13,4)
         antel(13,4) = anD4 + antel(13,3)
         if (antel(13,4).LE.antal) then 
            do 1341 plusD4 = min(anD4,10),max(anD4-8,0),-1
               ansats(13,4,1) = plusD4
               ansats(13,4,0) = anD4-plusD4
*     13h
      call Slug(13,5,varmax,varupp,varned,ansats,org,lock(13,5),
     :                      dubbel,low,start(13,5),stopp(13,5))
      do 1350 anD5 = start(13,5),stopp(13,5),steg(13,5)
         antel(13,5) = anD5 + antel(13,4)
         if (antel(13,5).LE.antal .AND. ansats(13,4,1).LE.2) then 
            do 1351 plusD5 = min(anD5,12),max(anD5-10,0),-1
               ansats(13,5,1) = plusD5
               ansats(13,5,0) = anD5-plusD5
*     13i
      call Slug(13,6,varmax,varupp,varned,ansats,org,lock(13,6),
     :                      dubbel,low,start(13,6),stopp(13,6))
      do 1360 anD6 = start(13,6),stopp(13,6),steg(13,6)
         antel(13,6) = anD6 + antel(13,5)
         if (antel(13,6).LE.antal .AND. ansats(13,5,1).LE.2 .AND.
     :                                    ansats(13,5,0).LE.2) then 
            do 1361 plusD6 = min(anD6,14),max(anD6-12,0),-1
               ansats(13,6,1) = plusD6
               ansats(13,6,0) = anD6-plusD6
*     13k
      call Slug(13,7,varmax,varupp,varned,ansats,org,lock(13,7),
     :                      dubbel,low,start(13,7),stopp(13,7))
      do 1370 anD7 = start(13,7),stopp(13,7),steg(13,7)
         antel(13,7) = anD7 + antel(13,6)
         if (antel(13,7).LE.antal .AND. ansats(13,6,1).LE.2 .AND.
     :                                    ansats(13,6,0).LE.2) then 
            do 1371 plusD7 = min(anD7,16),max(anD7-14,0),-1
               ansats(13,7,1) = plusD7
               ansats(13,7,0) = anD7-plusD7
*     13l
      call Slug(13,8,varmax,varupp,varned,ansats,org,lock(13,8),
     :                      dubbel,low,start(13,8),stopp(13,8))
      do 1380 anD8 = start(13,8),stopp(13,8),steg(13,8)
         antel(13,8) = anD8 + antel(13,7)
         if (antel(13,8).LE.antal .AND. ansats(13,7,1).LE.2 .AND.
     :                                    ansats(13,7,0).LE.2) then 
            do 1381 plusD8 = min(anD8,18),max(anD8-16,0),-1
               ansats(13,8,1) = plusD8
               ansats(13,8,0) = anD8-plusD8
*     13m
      call Slug(13,9,varmax,varupp,varned,ansats,org,lock(13,9),
     :                      dubbel,low,start(13,9),stopp(13,9))
      do 1390 anD9 = start(13,9),stopp(13,9),steg(13,9)
         antel(13,9) = anD9 + antel(13,8)
         if (antel(13,9).LE.antal .AND. ansats(13,8,1).LE.2 .AND.
     :                                    ansats(13,8,0).LE.2) then 
            do 1391 plusD9 = min(anD9,20),max(anD9-18,0),-1
               ansats(13,9,1) = plusD9
               ansats(13,9,0) = anD9-plusD9
*     13n
      call Slug(13,10,varmax,varupp,varned,ansats,org,lock(13,10),
     :                      dubbel,low,start(13,10),stopp(13,10))
      do 13100 anDA = start(13,10),stopp(13,10),steg(13,10)
         antel(13,10) = anDA + antel(13,9)
         if (antel(13,10).LE.antal .AND. ansats(13,9,1).LE.2 .AND.
     :       ansats(13,9,0).LE.2 .AND. antel(13,10).GE.lim(13)) then 
            do 13101 plusDA = min(anDA,22),max(anDA-20,0),-1
               ansats(13,10,1) = plusDA
               ansats(13,10,0) = anDA-plusDA
*     14s
      call Slug(14,0,varmax,varupp,varned,ansats,org,lock(14,0),
     :                      dubbel,low,start(14,0),stopp(14,0))
      do 1400 anE0 = start(14,0),stopp(14,0),steg(14,0)
         antel(14,0) = anE0 + antel(13,10)
         if (antel(14,0).LE.antal .AND. ansats(13,10,1).LE.2 .AND.
     :                                   ansats(13,10,0).LE.2) then 
            ansats(14,0,0) = anE0
*     14p
      call Slug(14,1,varmax,varupp,varned,ansats,org,lock(14,1),
     :                      dubbel,low,start(14,1),stopp(14,1))
      do 1410 anE1 = start(14,1),stopp(14,1),steg(14,1)
         antel(14,1) = anE1 + antel(14,0)
         if (antel(14,1).LE.antal) then 
            do 1411 plusE1 = min(anE1,4),max(anE1-2,0),-1
               ansats(14,1,1) = plusE1
               ansats(14,1,0) = anE1-plusE1
*     14d
      call Slug(14,2,varmax,varupp,varned,ansats,org,lock(14,2),
     :                      dubbel,low,start(14,2),stopp(14,2))
      do 1420 anE2 = start(14,2),stopp(14,2),steg(14,2)
         antel(14,2) = anE2 + antel(14,1)
         if (antel(14,2).LE.antal) then 
            do 1421 plusE2 = min(anE2,6),max(anE2-4,0),-1
               ansats(14,2,1) = plusE2
               ansats(14,2,0) = anE2-plusE2
*     14f
      call Slug(14,3,varmax,varupp,varned,ansats,org,lock(14,3),
     :                      dubbel,low,start(14,3),stopp(14,3))
      do 1430 anE3 = start(14,3),stopp(14,3),steg(14,3)
         antel(14,3) = anE3 + antel(14,2)
         if (antel(14,3).LE.antal) then 
            do 1431 plusE3 = min(anE3,8),max(anE3-6,0),-1
               ansats(14,3,1) = plusE3
               ansats(14,3,0) = anE3-plusE3
*     14g
      call Slug(14,4,varmax,varupp,varned,ansats,org,lock(14,4),
     :                      dubbel,low,start(14,4),stopp(14,4))
      do 1440 anE4 = start(14,4),stopp(14,4),steg(14,4)
         antel(14,4) = anE4 + antel(14,3)
         if (antel(14,4).LE.antal) then 
            do 1441 plusE4 = min(anE4,10),max(anE4-8,0),-1
               ansats(14,4,1) = plusE4
               ansats(14,4,0) = anE4-plusE4
*     14h
      call Slug(14,5,varmax,varupp,varned,ansats,org,lock(14,5),
     :                      dubbel,low,start(14,5),stopp(14,5))
      do 1450 anE5 = start(14,5),stopp(14,5),steg(14,5)
         antel(14,5) = anE5 + antel(14,4)
         if (antel(14,5).LE.antal .AND. ansats(14,4,1).LE.2) then 
            do 1451 plusE5 = min(anE5,12),max(anE5-10,0),-1
               ansats(14,5,1) = plusE5
               ansats(14,5,0) = anE5-plusE5
*     14i
      call Slug(14,6,varmax,varupp,varned,ansats,org,lock(14,6),
     :                      dubbel,low,start(14,6),stopp(14,6))
      do 1460 anE6 = start(14,6),stopp(14,6),steg(14,6)
         antel(14,6) = anE6 + antel(14,5)
         if (antel(14,6).LE.antal .AND. ansats(14,5,1).LE.2 .AND.
     :       ansats(14,5,0).LE.2) then 
            do 1461 plusE6 = min(anE6,14),max(anE6-12,0),-1
               ansats(14,6,1) = plusE6
               ansats(14,6,0) = anE6-plusE6
*     14k
      call Slug(14,7,varmax,varupp,varned,ansats,org,lock(14,7),
     :                      dubbel,low,start(14,7),stopp(14,7))
      do 1470 anE7 = start(14,7),stopp(14,7),steg(14,7)
         antel(14,7) = anE7 + antel(14,6)
         if (antel(14,7).LE.antal .AND. ansats(14,6,1).LE.2 .AND.
     :                                    ansats(14,6,0).LE.2) then 
            do 1471 plusE7 = min(anE7,16),max(anE7-14,0),-1
               ansats(14,7,1) = plusE7
               ansats(14,7,0) = anE7-plusE7
*     14l
      call Slug(14,8,varmax,varupp,varned,ansats,org,lock(14,8),
     :                      dubbel,low,start(14,8),stopp(14,8))
      do 1480 anE8 = start(14,8),stopp(14,8),steg(14,8)
         antel(14,8) = anE8 + antel(14,7)
         if (antel(14,8).LE.antal .AND. ansats(14,7,1).LE.2 .AND.
     :                                    ansats(14,7,0).LE.2) then 
            do 1481 plusE8 = min(anE8,18),max(anE8-16,0),-1
               ansats(14,8,1) = plusE8
               ansats(14,8,0) = anE8-plusE8
*     14m
      call Slug(14,9,varmax,varupp,varned,ansats,org,lock(14,9),
     :                      dubbel,low,start(14,9),stopp(14,9))
      do 1490 anE9 = start(14,9),stopp(14,9),steg(14,9)
         antel(14,9) = anE9 + antel(14,8)
         if (antel(14,9).LE.antal .AND. ansats(14,8,1).LE.2 .AND.
     :                                    ansats(14,8,0).LE.2) then 
            do 1491 plusE9 = min(anE9,20),max(anE9-18,0),-1
               ansats(14,9,1) = plusE9
               ansats(14,9,0) = anE9-plusE9
*     14n
      call Slug(14,10,varmax,varupp,varned,ansats,org,lock(14,10),
     :                      dubbel,low,start(14,10),stopp(14,10))
      do 14100 anEA = start(14,10),stopp(14,10),steg(14,10)
         antel(14,10) = anEA + antel(14,9)
         if (antel(14,10).LE.antal .AND. ansats(14,9,1).LE.2 .AND.
     :          ansats(14,9,0).LE.2 .AND. antel(14,10).GE.lim(14)) then 
            do 14101 plusEA = min(anEA,22),max(anEA-20,0),-1
               ansats(14,10,1) = plusEA
               ansats(14,10,0) = anEA-plusEA
*     15s
      call Slug(15,0,varmax,varupp,varned,ansats,org,lock(15,0),
     :                      dubbel,low,start(15,0),stopp(15,0))
      do 1500 anF0 = start(15,0),stopp(15,0),steg(15,0)
         antel(15,0) = anF0 + antel(14,10)
         if (antel(15,0).LE.antal .AND. ansats(14,10,1).LE.2 .AND.
     :                                   ansats(14,10,0).LE.2) then 
            ansats(15,0,0) = anF0
*     15p
      call Slug(15,1,varmax,varupp,varned,ansats,org,lock(15,1),
     :                      dubbel,low,start(15,1),stopp(15,1))
      do 1510 anF1 = start(15,1),stopp(15,1),steg(15,1)
         antel(15,1) = anF1 + antel(15,0)
         if (antel(15,1).LE.antal) then 
            do 1511 plusF1 = min(anF1,4),max(anF1-2,0),-1
               ansats(15,1,1) = plusF1
               ansats(15,1,0) = anF1-plusF1
*     15d
      call Slug(15,2,varmax,varupp,varned,ansats,org,lock(15,2),
     :                      dubbel,low,start(15,2),stopp(15,2))
      do 1520 anF2 = start(15,2),stopp(15,2),steg(15,2)
         antel(15,2) = anF2 + antel(15,1)
         if (antel(15,2).LE.antal) then 
            do 1521 plusF2 = min(anF2,6),max(anF2-4,0),-1
               ansats(15,2,1) = plusF2
               ansats(15,2,0) = anF2-plusF2
*     15f
      call Slug(15,3,varmax,varupp,varned,ansats,org,lock(15,3),
     :                      dubbel,low,start(15,3),stopp(15,3))
      do 1530 anF3 = start(15,3),stopp(15,3),steg(15,3)
         antel(15,3) = anF3 + antel(15,2)
         if (antel(15,3).LE.antal) then 
            do 1531 plusF3 = min(anF3,8),max(anF3-6,0),-1
               ansats(15,3,1) = plusF3
               ansats(15,3,0) = anF3-plusF3
*     15g
      call Slug(15,4,varmax,varupp,varned,ansats,org,lock(15,4),
     :                      dubbel,low,start(15,4),stopp(15,4))
      do 1540 anF4 = start(15,4),stopp(15,4),steg(15,4)
         antel(15,4) = anF4 + antel(15,3)
         if (antel(15,4).LE.antal) then 
            do 1541 plusF4 = min(anF4,10),max(anF4-8,0),-1
               ansats(15,4,1) = plusF4
               ansats(15,4,0) = anF4-plusF4
*     15h
      call Slug(15,5,varmax,varupp,varned,ansats,org,lock(15,5),
     :                      dubbel,low,start(15,5),stopp(15,5))
      do 1550 anF5 = start(15,5),stopp(15,5),steg(15,5)
         antel(15,5) = anF5 + antel(15,4)
         if (antel(15,5).LE.antal .AND. ansats(15,4,1).LE.2) then 
            do 1551 plusF5 = min(anF5,12),max(anF5-10,0),-1
               ansats(15,5,1) = plusF5
               ansats(15,5,0) = anF5-plusF5
*     15i
      call Slug(15,6,varmax,varupp,varned,ansats,org,lock(15,6),
     :                      dubbel,low,start(15,6),stopp(15,6))
      do 1560 anF6 = start(15,6),stopp(15,6),steg(15,6)
         antel(15,6) = anF6 + antel(15,5)
         if (antel(15,6).LE.antal .AND. ansats(15,5,1).LE.2 .AND.
     :                                    ansats(15,5,0).LE.2) then 
            do 1561 plusF6 = min(anF6,14),max(anF6-12,0),-1
               ansats(15,6,1) = plusF6
               ansats(15,6,0) = anF6-plusF6
*     15k
      call Slug(15,7,varmax,varupp,varned,ansats,org,lock(15,7),
     :                      dubbel,low,start(15,7),stopp(15,7))
      do 1570 anF7 = start(15,7),stopp(15,7),steg(15,7)
         antel(15,7) = anF7 + antel(15,6)
         if (antel(15,7).LE.antal .AND. ansats(15,6,1).LE.2 .AND.
     :                                    ansats(15,6,0).LE.2) then 
            do 1571 plusF7 = min(anF7,16),max(anF7-14,0),-1
               ansats(15,7,1) = plusF7
               ansats(15,7,0) = anF7-plusF7
*     15l
      call Slug(15,8,varmax,varupp,varned,ansats,org,lock(15,8),
     :                      dubbel,low,start(15,8),stopp(15,8))
      do 1580 anF8 = start(15,8),stopp(15,8),steg(15,8)
         antel(15,8) = anF8 + antel(15,7)
         if (antel(15,8).LE.antal .AND. ansats(15,7,1).LE.2 .AND.
     :                                    ansats(15,7,0).LE.2) then
            do 1581 plusF8 = min(anF8,18),max(anF8-16,0),-1
               ansats(15,8,1) = plusF8
               ansats(15,8,0) = anF8-plusF8
*     15m
      call Slug(15,9,varmax,varupp,varned,ansats,org,lock(15,9),
     :                      dubbel,low,start(15,9),stopp(15,9))
      do 1590 anF9 = start(15,9),stopp(15,9),steg(15,9)
         antel(15,9) = anF9 + antel(15,8)
         if (antel(15,9).LE.antal .AND. ansats(15,8,1).LE.2 .AND.
     :                                    ansats(15,8,0).LE.2) then
            do 1591 plusF9 = min(anF9,20),max(anF9-18,0),-1
               ansats(15,9,1) = plusF9
               ansats(15,9,0) = anF9-plusF9
*     15n
      call Slug(15,10,varmax,varupp,varned,ansats,org,lock(15,10),
     :                      dubbel,low,start(15,10),stopp(15,10))
      do 15100 anFA = start(15,10),stopp(15,10),steg(15,10)
         antel(15,10) = anFA + antel(15,9)
         if (antel(15,10).EQ.antal .AND. ansats(15,9,1).LE.2 .AND.
     :                                    ansats(15,9,0).LE.2) then
            do 15101 plusFA = min(anFA,22),max(anFA-20,0),-1
               ansats(15,10,1) = plusFA
               ansats(15,10,0) = anFA-plusFA
               if (ansats(15,10,1).LE.2 .AND. ansats(15,10,0).LE.2) then
                  par = 0
                  elar = 0
                  do 6 i=1,15
                     do 6 j=0,min(10,i-1)
                        do 6 k=0,min(j,1)
                           elar = elar + ansats(i,j,k)
    6                      par = mod(par+j*ansats(i,j,k),2)
                  if (par.EQ.par0) then
                     if (elar.EQ.antal) then 
                        call 
     :               Gen(ansats,posn,posl,skal,cf,first,minJ,maxJ,par0)
                     else
                        write(*,*) 'FEL'
                     endif
                  endif
               endif
15101       continue
         endif
15100 continue
 1591 continue
      endif
 1590 continue
 1581 continue
      endif
 1580 continue
 1571 continue
      endif
 1570 continue
 1561 continue
      endif
 1560 continue
 1551 continue
      endif
 1550 continue
 1541 continue
      endif
 1540 continue
 1531 continue
      endif
 1530 continue
 1521 continue
      endif
 1520 continue
 1511 continue
      endif
 1510 continue
      endif
 1500 continue
14101 continue
      endif
14100 continue
 1491 continue
      endif
 1490 continue
 1481 continue
      endif
 1480 continue
 1471 continue
      endif
 1470 continue
 1461 continue
      endif
 1460 continue
 1451 continue
      endif
 1450 continue
 1441 continue
      endif
 1440 continue
 1431 continue
      endif
 1430 continue
 1421 continue
      endif
 1420 continue
 1411 continue
      endif
 1410 continue
      endif
 1400 continue
13101 continue
      endif
13100 continue
 1391 continue
      endif
 1390 continue
 1381 continue
      endif
 1380 continue
 1371 continue
      endif
 1370 continue
 1361 continue
      endif
 1360 continue
 1351 continue
      endif
 1350 continue
 1341 continue
      endif
 1340 continue
 1331 continue
      endif
 1330 continue
 1321 continue
      endif
 1320 continue
 1311 continue
      endif
 1310 continue
      endif
 1300 continue
12101 continue
      endif
12100 continue
 1291 continue
      endif
 1290 continue
 1281 continue
      endif
 1280 continue
 1271 continue
      endif
 1270 continue
 1261 continue
      endif
 1260 continue
 1251 continue
      endif
 1250 continue
 1241 continue
      endif
 1240 continue
 1231 continue
      endif
 1230 continue
 1221 continue
      endif
 1220 continue
 1211 continue
      endif
 1210 continue
      endif
 1200 continue
11101 continue
      endif
11100 continue
 1191 continue
      endif
 1190 continue
 1181 continue
      endif
 1180 continue
 1171 continue
      endif
 1170 continue
 1161 continue
      endif
 1160 continue
 1151 continue
      endif
 1150 continue
 1141 continue
      endif
 1140 continue
 1131 continue
      endif
 1130 continue
 1121 continue
      endif
 1120 continue
 1111 continue
      endif
 1110 continue
      endif
 1100 continue
 1091 continue
      endif
 1090 continue
 1081 continue
      endif
 1080 continue
 1071 continue
      endif
 1070 continue
 1061 continue
      endif
 1060 continue
 1051 continue
      endif
 1050 continue
 1041 continue
      endif
 1040 continue
 1031 continue
      endif
 1030 continue
 1021 continue
      endif
 1020 continue
 1011 continue
      endif
 1010 continue
      endif
 1000 continue
  981 continue
      endif
  980 continue
  971 continue
      endif
  970 continue
  961 continue
      endif
  960 continue
  951 continue
      endif
  950 continue
  941 continue
      endif
  940 continue
  931 continue
      endif
  930 continue
  921 continue
      endif
  920 continue
  911 continue
      endif
  910 continue
      endif
  900 continue
  871 continue
      endif
  870 continue
  861 continue
      endif
  860 continue
  851 continue
      endif
  850 continue
  841 continue
      endif
  840 continue
  831 continue
      endif
  830 continue
  821 continue
      endif
  820 continue
  811 continue
      endif
  810 continue
      endif
  800 continue
  761 continue
      endif
  760 continue
  751 continue
      endif
  750 continue
  741 continue
      endif
  740 continue
  731 continue
      endif
  730 continue
  721 continue
      endif
  720 continue
  711 continue
      endif
  710 continue
      endif
  700 continue
  651 continue
      endif
  650 continue
  641 continue
      endif
  640 continue
  631 continue
      endif
  630 continue
  621 continue
      endif
  620 continue
  611 continue
      endif
  610 continue
      endif
  600 continue
  541 continue
      endif
  540 continue
  531 continue
      endif
  530 continue
  521 continue
      endif
  520 continue
  511 continue
      endif
  510 continue
      endif
  500 continue
  431 continue
      endif
  430 continue
  421 continue
      endif
  420 continue
  411 continue
      endif
  410 continue
      endif
  400 continue
  321 continue
      endif
  320 continue
  311 continue
      endif
  310 continue
      endif
  300 continue
  211 continue
      endif
  210 continue
      endif
  200 continue
      endif
  100 continue
      if (first) then
         rewind(fil_1)
      else
         rewind(fil_2)
      endif
      if (cf.EQ.0) then
         write(*,1005) 'No configuration state has been generated.'
      elseif (cf.EQ.1) then
         write(*,1005) 'One configuration state has been generated.'
      elseif (cf.LT.10) then
         write(*,1001) cf,' configuration states have been generated.'
      elseif (cf.LT.100) then
         write(*,1002) cf,' configuration states have been generated.'
      elseif (cf.LT.1000) then
         write(*,1003) cf,' configuration states have been generated.'
      elseif (cf.LT.10000) then
         write(*,1004) cf,' configuration states have been generated.'
      elseif (cf.LT.100000) then
         write(*,1006) cf,' configuration states have been generated.'
      else
         write(*,*) cf,' configuration states have been generated.'
      endif
C 1000 format(A)
 1001 format(' ',I1,A)
 1002 format(' ',I2,A)
 1003 format(' ',I3,A)
 1004 format(' ',I4,A)
 1005 format(' ',A)
 1006 format(' ',I5,A)
 5000 format(11I2)
      return
      end
