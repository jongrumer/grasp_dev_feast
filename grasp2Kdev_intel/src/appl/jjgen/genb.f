*     last edited July 31, 1996
      subroutine Gen(ansats,posn,posl,skal,cf,first,minJ,maxJ,par) 
      integer fil_1,fil_2
      parameter (fil_1=7, fil_2=8)
      integer ansats(1:15,0:10,0:1),koppl(0:10,0:11,0:1)
      integer JKVANT(0:10,0:1,0:5,20),antmax(0:10,0:1),par
      integer minJ,maxJ,pos,i,n,l,k,J(20),JK(20),orbit(20),antel(20)
      integer i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16
      integer i17,i18,i19,i20,plus(20),S(20),resJ
      integer JK1,JK2,JK3,JK4,JK5,JK6,JK7,JK8,JK9,JK10,JK11,JK12,JK13
      integer JK14,JK15,JK16,JK17,JK18,fil,antko(20),skal,cf
      integer posn(110),posl(110),SENIOR(0:10,0:1,0:5,20),n1,n10
      character rad1*200,rad2*200,rad3*200
      character*2 L1(0:10,0:1)
      logical first
      data (L1(i,0),i=0,10)
     :          /'s ','p-','d-','f-','g-','h-','i-','k-','l-','m-','n-'/
      data (L1(i,1),i=0,10)
     :          /'s ','p ','d ','f ','g ','h ','i ','k ','l ','m ','n '/
*     The value of antmax(l-number,x) is the maximum number of electrons
*     in the orbital, x represents +/- coupling of s- and l- number
      data (antmax(i,0),i=0,10)/2,2,4,6,8,10,12,14,16,18,20/
      data (antmax(i,1),i=0,10)/0,4,6,8,10,12,14,16,18,20,22/
*     The value of koppl(l-number,number of electrons,x) is the number of
*     possible couplings for a certain orbital. If the orbital is
*     populated with more than half of the maximal number of electrons
*     the index "number of electrons" should be substituted with
*     "antmax(l-number) - number of electrons".
      data (koppl(0,i,0),i=0,1)       / 1, 1/
*     l=0
      data (koppl(1,i,0),i=0,1)       / 1, 1/
      data (koppl(1,i,1),i=0,2)       / 1, 1, 2/
*     l=1
      data (koppl(2,i,0),i=0,2)       / 1, 1, 2/
      data (koppl(2,i,1),i=0,3)       / 1, 1, 3, 3/
*     l=2
      data (koppl(3,i,0),i=0,3)       / 1, 1, 3, 3/
      data (koppl(3,i,1),i=0,4)       / 1, 1, 4, 6, 8/
*     l=3
      data (koppl(4,i,0),i=0,4)       / 1, 1, 4, 6, 8/
      data (koppl(4,i,1),i=0,5)       / 1, 1, 5,10,16,20/
*     l=4
      data (koppl(5,i,0),i=0,5)       / 1, 1, 5,10,16,20/
      data (koppl(5,i,1),i=0,2)       / 1, 1, 6/
*     l=5 
      data (koppl(6,i,0),i=0,2)       / 1, 1, 6/
      data (koppl(6,i,1),i=0,2)       / 1, 1, 7/
*     l=6 
      data (koppl(7,i,0),i=0,2)       / 1, 1, 7/
      data (koppl(7,i,1),i=0,2)       / 1, 1, 8/
*     l=7 
      data (koppl(8,i,0),i=0,2)       / 1, 1, 8/
      data (koppl(8,i,1),i=0,2)       / 1, 1, 9/
*     l=8 
      data (koppl(9,i,0),i=0,2)       / 1, 1, 9/
      data (koppl(9,i,1),i=0,2)       / 1, 1,10/
*     l=9
      data (koppl(10,i,0),i=0,2)      / 1, 1,10/
      data (koppl(10,i,1),i=0,2)      / 1, 1,11/ 
*     l=10

*  JKVANT(l-number, +/-, number of electrons, coupling number) is 2*J-number

      data  JKVANT(0,0,0,1)           / 0/     
*     data  SENIOR(0,0,0,1)           / 0/
      data  SENIOR(0,0,0,1)           /-1/
*     l=0 #=0
      data  JKVANT(0,0,1,1)           / 1/
*     data  SENIOR(0,0,1,1)           / 1/
      data  SENIOR(0,0,1,1)           /-1/
*     l=0 #=1
      data  JKVANT(1,0,0,1)           / 0/
*     data  SENIOR(1,0,0,1)           / 0/
      data  SENIOR(1,0,0,1)           /-1/
*     l=1 #=0 -
      data  JKVANT(1,0,1,1)           / 1/
*     data  SENIOR(1,0,1,1)           / 1/
      data  SENIOR(1,0,1,1)           /-1/
*     l=1 #=1 -
      data  JKVANT(1,1,0,1)           / 0/
*     data  SENIOR(1,1,0,1)           / 0/
      data  SENIOR(1,1,0,1)           /-1/
*     l=1 #=0 +
      data  JKVANT(1,1,1,1)           / 3/
*     data  SENIOR(1,1,1,1)           / 1/
      data  SENIOR(1,1,1,1)           /-1/
*     l=1 #=1 +
      data (JKVANT(1,1,2,i),i=1,2)    / 0, 4/
*     data (SENIOR(1,1,2,i),i=1,2)    / 0, 2/
      data (SENIOR(1,1,2,i),i=1,2)    /-1,-1/
*     l=1 #=2 +
      data  JKVANT(2,0,0,1)           / 0/
*     data  SENIOR(2,0,0,1)           / 0/
      data  SENIOR(2,0,0,1)           /-1/
*     l=2 #=0 -
      data  JKVANT(2,0,1,1)           / 3/
*     data  SENIOR(2,0,1,1)           / 1/
      data  SENIOR(2,0,1,1)           /-1/
*     l=2 #=1 -
      data (JKVANT(2,0,2,i),i=1,2)    / 0, 4/
*     data (SENIOR(2,0,2,i),i=1,2)    / 0, 2/
      data (SENIOR(2,0,2,i),i=1,2)    /-1,-1/
*     l=2 #=2 -
      data  JKVANT(2,1,0,1)           / 0/
*     data  SENIOR(2,1,0,1)           / 0/
      data  SENIOR(2,1,0,1)           /-1/
*     l=2 #=0 +
      data  JKVANT(2,1,1,1)           / 5/
*     data  SENIOR(2,1,1,1)           / 1/
      data  SENIOR(2,1,1,1)           /-1/
*     l=2 #=1 +
      data (JKVANT(2,1,2,i),i=1,3)    / 0, 4, 8/
*     data (SENIOR(2,1,2,i),i=1,3)    / 0, 2, 2/
      data (SENIOR(2,1,2,i),i=1,3)    /-1,-1,-1/
*     l=2 #=2 +
      data (JKVANT(2,1,3,i),i=1,3)    / 5, 3, 9/
*     data (SENIOR(2,1,3,i),i=1,3)    / 1, 3, 3/
      data (SENIOR(2,1,3,i),i=1,3)    /-1,-1,-1/
*     l=2 #=3 +
      data  JKVANT(3,0,0,1)           / 0/
*     data  SENIOR(3,0,0,1)           / 0/
      data  SENIOR(3,0,0,1)           /-1/
*     l=3 #=0 -
      data  JKVANT(3,0,1,1)           / 5/
*     data  SENIOR(3,0,1,1)           / 1/
      data  SENIOR(3,0,1,1)           /-1/
*     l=3 #=1 -
      data (JKVANT(3,0,2,i),i=1,3)    / 0, 4, 8/
*     data (SENIOR(3,0,2,i),i=1,3)    / 0, 2, 2/
      data (SENIOR(3,0,2,i),i=1,3)    /-1,-1,-1/
*     l=3 #=2 -
      data (JKVANT(3,0,3,i),i=1,3)    / 5, 3, 9/
*     data (SENIOR(3,0,3,i),i=1,3)    / 1, 3, 3/
      data (SENIOR(3,0,3,i),i=1,3)    /-1,-1,-1/
*     l=3 #=3 -
      data  JKVANT(3,1,0,1)           / 0/
*     data  SENIOR(3,1,0,1)           / 0/
      data  SENIOR(3,1,0,1)           /-1/
*     l=3 #=0 +
      data  JkVANT(3,1,1,1)           / 7/
*     data  SENIOR(3,1,1,1)           / 1/
      data  SENIOR(3,1,1,1)           /-1/
*     l=3 #=1 +
      data (JKVANT(3,1,2,i),i=1,4)    / 0, 4, 8,12/
*     data (SENIOR(3,1,2,i),i=1,4)    / 0, 2, 2, 2/
      data (SENIOR(3,1,2,i),i=1,4)    /-1,-1,-1,-1/
*     l=3 #=2 +
      data (JKVANT(3,1,3,i),i=1,6)    / 7, 3, 5, 9,11,15/
*     data (SENIOR(3,1,3,i),i=1,6)    / 1, 3, 3, 3, 3, 3/
      data (SENIOR(3,1,3,i),i=1,6)    /-1,-1,-1,-1,-1,-1/
*     l=3 #=3 +
      data (JKVANT(3,1,4,i),i=1,8)    / 0, 4, 8,12, 4, 8,10,16/
*     data (SENIOR(3,1,4,i),i=1,8)    / 0, 2, 2, 2, 4, 4, 4, 4/
      data (SENIOR(3,1,4,i),i=1,8)    /-1, 2, 2,-1, 4, 4,-1,-1/
*     l=3 #=4 +
      data  JKVANT(4,0,0,1)           / 0/
*     data  SENIOR(4,0,0,1)           / 0/
      data  SENIOR(4,0,0,1)           /-1/
*     l=4 #=0 -
      data  JkVANT(4,0,1,1)           / 7/
*     data  SENIOR(4,0,1,1)           / 1/
      data  SENIOR(4,0,1,1)           /-1/
*     l=4 #=1 -
      data (JKVANT(4,0,2,i),i=1,4)    / 0, 4, 8,12/
*     data (SENIOR(4,0,2,i),i=1,4)    / 0, 2, 2, 2/
      data (SENIOR(4,0,2,i),i=1,4)    /-1,-1,-1,-1/
*     l=4 #=2 -
      data (JKVANT(4,0,3,i),i=1,6)    / 7, 3, 5, 9,11,15/
*     data (SENIOR(4,0,3,i),i=1,6)    / 1, 3, 3, 3, 3, 3/
      data (SENIOR(4,0,3,i),i=1,6)    /-1,-1,-1,-1,-1,-1/
*     l=4 #=3 -
      data (JKVANT(4,0,4,i),i=1,8)    / 0, 4, 8,12, 4, 8,10,16/
*     data (SENIOR(4,0,4,i),i=1,8)    / 0, 2, 2, 2, 4, 4, 4, 4/
      data (SENIOR(4,0,4,i),i=1,8)    /-1, 2, 2,-1, 4, 4,-1,-1/
*     l=4 #=4 -
      data  JKVANT(4,1,0,1)           / 0/
*     data  SENIOR(4,1,0,1)           / 0/
      data  SENIOR(4,1,0,1)           /-1/
*     l=4 #=0 +
      data  JKVANT(4,1,1,1)           / 9/
*     data  SENIOR(4,1,1,1)           / 1/
      data  SENIOR(4,1,1,1)           /-1/
*     l=4 #=1 +
      data (JKVANT(4,1,2,i),i=1,5)    / 0, 4, 8,12,16/
*     data (SENIOR(4,1,2,i),i=1,5)    / 0, 2, 2, 2, 2/
      data (SENIOR(4,1,2,i),i=1,5)    /-1,-1,-1,-1,-1/
*     l=4 #=2 +
      data (JKVANT(4,1,3,i),i=1,10)   / 9, 3, 5, 7, 9,11,13,15,17,21/
*     data (SENIOR(4,1,3,i),i=1,10)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3/
      data (SENIOR(4,1,3,i),i=1,10)   / 1,-1,-1,-1, 3,-1,-1,-1,-1,-1/
*     l=4 #=3 + 
      data (JKVANT(4,1,4,i),i=1,16)   / 0, 4, 8,12,16, 0, 4, 6, 8,10,12, 
     :                                 14,16,18,20,24/
*     data (SENIOR(4,1,4,i),i=1,16)   / 0, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 
*    :                                  4, 4, 4, 4, 4/
      data (SENIOR(4,1,4,i),i=1,16)   / 0, 2, 2, 2, 2, 4, 4,-1, 4,-1, 4, 
     :                                 -1, 4,-1,-1,-1/
*     l=4 #=4 + 
      data (JKVANT(4,1,5,i),i=1,20)   / 9, 3, 5, 7, 9,11,13,15,17,21, 1,
     :                                  5, 7, 9, 11,13,15,17,19,25/
*     data (SENIOR(4,1,5,i),i=1,20)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5,
*    :                                  5, 5, 5, 5, 5, 5, 5, 5, 5/
      data (SENIOR(4,1,5,i),i=1,20)   / 1,-1, 3, 3, 3, 3, 3, 3, 3,-1,-1,
     :                                  5, 5, 5, 5, 5, 5, 5,-1,-1/
*     l=4 #=5 +
      data  JKVANT(5,0,0,1)           / 0/
*     data  SENIOR(5,0,0,1)           / 0/
      data  SENIOR(5,0,0,1)           /-1/
*     l=5 #=0 -
      data  JKVANT(5,0,1,1)           / 9/
*     data  SENIOR(5,0,1,1)           / 1/
      data  SENIOR(5,0,1,1)           /-1/
*     l=5 #=1 -
      data (JKVANT(5,0,2,i),i=1,5)    / 0, 4, 8,12,16/
*     data (SENIOR(5,0,2,i),i=1,5)    / 0, 2, 2, 2, 2/
      data (SENIOR(5,0,2,i),i=1,5)    /-1,-1,-1,-1,-1/
*     l=5 #=2 -
      data (JKVANT(5,0,3,i),i=1,10)   / 9, 3, 5, 7, 9,11,13,15,17,21/
*     data (SENIOR(5,0,3,i),i=1,10)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3/
      data (SENIOR(5,0,3,i),i=1,10)   / 1,-1,-1,-1, 3,-1,-1,-1,-1,-1/
*     l=5 #=3 - 
      data (JKVANT(5,0,4,i),i=1,16)   / 0, 4, 8,12,16, 0, 4, 6, 8,10,12, 
     :                                 14,16,18,20,24/
*     data (SENIOR(5,0,4,i),i=1,16)   / 0, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 
*    :                                  4, 4, 4, 4, 4/
      data (SENIOR(5,0,4,i),i=1,16)   / 0, 2, 2, 2, 2, 4, 4,-1, 4,-1, 4,
     :                                 -1, 4,-1,-1,-1/ 
*     l=5 #=4 - 
      data (JKVANT(5,0,5,i),i=1,20)   / 9, 3, 5, 7, 9,11,13,15,17,21, 1,
     :                                  5, 7, 9, 11,13,15,17,19,25/
*     data (SENIOR(5,0,5,i),i=1,20)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5,
*    :                                  5, 5, 5, 5, 5, 5, 5, 5, 5/
      data (SENIOR(5,0,5,i),i=1,20)   / 1,-1, 3, 3, 3, 3, 3, 3, 3,-1,-1,
     :                                  5, 5, 5, 5, 5, 5, 5,-1,-1/
*     l=5 #=5 - 
      data  JKVANT(5,1,0,1)           / 0/
*     data  SENIOR(5,1,0,1)           / 0/
      data  SENIOR(5,1,0,1)           /-1/
*     l=5 #=0 +
      data  JKVANT(5,1,1,1)           /11/
*     data  SENIOR(5,1,1,1)           / 1/
      data  SENIOR(5,1,1,1)           /-1/
*     l=5 #=1 +
      data (JKVANT(5,1,2,i),i=1,6)    / 0, 4, 8,12,16,20/
*     data (SENIOR(5,1,2,i),i=1,6)    / 0, 2, 2, 2, 2, 2/
      data (SENIOR(5,1,2,i),i=1,6)    /-1,-1,-1,-1,-1,-1/
*     l=5 #=2 +
      data  JKVANT(6,0,0,1)           / 0/
c     data  SENIOR(6,0,0,1)           / 0/
      data  SENIOR(6,0,0,1)           /-1/
*     l=6 #=0 -
      data  JKVANT(6,0,1,1)           /11/
c     data  SENIOR(6,0,1,1)           / 1/
      data  SENIOR(6,0,1,1)           /-1/
*     l=6 #=1 -
      data (JKVANT(6,0,2,i),i=1,6)    / 0, 4, 8,12,16,20/
c     data (SENIOR(6,0,2,i),i=1,6)    / 0, 2, 2, 2, 2, 2/
      data (SENIOR(6,0,2,i),i=1,6)    /-1,-1,-1,-1,-1,-1/
*     l=6 #=2 -
      data  JKVANT(6,1,0,1)           / 0/
c     data  SENIOR(6,1,0,1)           / 0/
      data  SENIOR(6,1,0,1)           /-1/
*     l=6 #=0 +
      data  JKVANT(6,1,1,1)           /13/
c     data  SENIOR(6,1,1,1)           / 1/
      data  SENIOR(6,1,1,1)           /-1/
*     l=6 #=1 +
      data (JKVANT(6,1,2,i),i=1,7)    / 0, 4, 8,12,16,20,24/
c     data (SENIOR(6,1,2,i),i=1,7)    / 0, 2, 2, 2, 2, 2, 2/
      data (SENIOR(6,1,2,i),i=1,7)    /-1,-1,-1,-1,-1,-1,-1/
*     l=6 #=2 +
      data  JKVANT(7,0,0,1)           / 0/
c     data  SENIOR(7,0,0,1)           / 0/
      data  SENIOR(7,0,0,1)           /-1/
*     l=7 #=0 -
      data  JKVANT(7,0,1,1)           /13/
c     data  SENIOR(7,0,1,1)           / 1/
      data  SENIOR(7,0,1,1)           /-1/
*     l=7 #=1 -
      data (JKVANT(7,0,2,i),i=1,7)    / 0, 4, 8,12,16,20,24/
c     data (SENIOR(7,0,2,i),i=1,7)    / 0, 2, 2, 2, 2, 2, 2/
      data (SENIOR(7,0,2,i),i=1,7)    /-1,-1,-1,-1,-1,-1,-1/
*     l=7 #=2 -
      data  JKVANT(7,1,0,1)           / 0/
c     data  SENIOR(7,1,0,1)           / 0/
      data  SENIOR(7,1,0,1)           /-1/
*     l=7 #=0 +
      data  JKVANT(7,1,1,1)           /15/
c     data  SENIOR(7,1,1,1)           / 1/
      data  SENIOR(7,1,1,1)           /-1/
*     l=7 #=1 +
      data (JKVANT(7,1,2,i),i=1,8)    / 0, 4, 8,12,16,20,24,28/
c     data (SENIOR(7,1,2,i),i=1,8)    / 0, 2, 2, 2, 2, 2, 2, 2/
      data (SENIOR(7,1,2,i),i=1,8)    /-1,-1,-1,-1,-1,-1,-1,-1/
*     l=7 #=2 +
      data  JKVANT(8,0,0,1)           / 0/
c     data  SENIOR(8,0,0,1)           / 0/
      data  SENIOR(8,0,0,1)           /-1/
*     l=8 #=0 -
      data  JKVANT(8,0,1,1)           /15/
c     data  SENIOR(8,0,1,1)           / 1/
      data  SENIOR(8,0,1,1)           /-1/
*     l=8 #=1 -
      data (JKVANT(8,0,2,i),i=1,8)    / 0, 4, 8,12,16,20,24,28/
c     data (SENIOR(8,0,2,i),i=1,8)    / 0, 2, 2, 2, 2, 2, 2, 2/
      data (SENIOR(8,0,2,i),i=1,8)    /-1,-1,-1,-1,-1,-1,-1,-1/
*     l=8 #=2 -
      data  JKVANT(8,1,0,1)           / 0/
c     data  SENIOR(8,1,0,1)           / 0/
      data  SENIOR(8,1,0,1)           /-1/
*     l=8 #=0 +
      data  JKVANT(8,1,1,1)           /17/
c     data  SENIOR(8,1,1,1)           / 1/
      data  SENIOR(8,1,1,1)           /-1/
*     l=8 #=1 +
      data (JKVANT(8,1,2,i),i=1,9)    / 0, 4, 8,12,16,20,24,28,32/
c     data (SENIOR(8,1,2,i),i=1,9)    / 0, 2, 2, 2, 2, 2, 2, 2, 2/
      data (SENIOR(8,1,2,i),i=1,9)    /-1,-1,-1,-1,-1,-1,-1,-1,-1/
*     l=8 #=2 +
      data  JKVANT(9,0,0,1)           / 0/
c     data  SENIOR(9,0,0,1)           / 0/
      data  SENIOR(9,0,0,1)           /-1/
*     l=9 #=0 -
      data  JKVANT(9,0,1,1)           /17/
c     data  SENIOR(9,0,1,1)           / 1/
      data  SENIOR(9,0,1,1)           /-1/
*     l=9 #=1 -
      data (JKVANT(9,0,2,i),i=1,9)    / 0, 4, 8,12,16,20,24,28,32/
c     data (SENIOR(9,0,2,i),i=1,9)    / 0, 2, 2, 2, 2, 2, 2, 2, 2/
      data (SENIOR(9,0,2,i),i=1,9)    /-1,-1,-1,-1,-1,-1,-1,-1,-1/
*     l=9 #=2 -
      data  JKVANT(9,1,0,1)           / 0/
c     data  SENIOR(9,1,0,1)           / 0/
      data  SENIOR(9,1,0,1)           /-1/
*     l=9 #=0 +
      data  JKVANT(9,1,1,1)           /19/
c     data  SENIOR(9,1,1,1)           / 1/
      data  SENIOR(9,1,1,1)           /-1/
*     l=9 #=1 +
      data (JKVANT(9,1,2,i),i=1,10)   / 0, 4, 8,12,16,20,24,28,32,36/
c     data (SENIOR(9,1,2,i),i=1,10)   / 0, 2, 2, 2, 2, 2, 2, 2, 2, 2/
      data (SENIOR(9,1,2,i),i=1,10)   /-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
*     l=9 #=2 +
      data  JKVANT(10,0,0,1)          / 0/
c     data  SENIOR(10,0,0,1)          / 0/
      data  SENIOR(10,0,0,1)          /-1/
*     l=10 #=0 -
      data  JKVANT(10,0,1,1)          /19/
c     data  SENIOR(10,0,1,1)          / 1/
      data  SENIOR(10,0,1,1)          /-1/
*     l=10 #=1 -
      data (JKVANT(10,0,2,i),i=1,10)   / 0, 4, 8,12,16,20,24,28,32,36/
c     data (SENIOR(10,0,2,i),i=1,10)   / 0, 2, 2, 2, 2, 2, 2, 2, 2, 2/
      data (SENIOR(10,0,2,i),i=1,10)   /-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
*     l=10 #=2 -
      data  JKVANT(10,1,0,1)          / 0/
c     data  SENIOR(10,1,0,1)          / 0/
      data  SENIOR(10,1,0,1)          /-1/
*     l=10 #=0 +
      data  JKVANT(10,1,1,1)          /21/
c     data  SENIOR(10,1,1,1)          / 1/
      data  SENIOR(10,1,1,1)          /-1/
*     l=10 #=1 +
      data (JKVANT(10,1,2,i),i=1,11)   / 0, 4, 8,12,16,20,24,28,32,36,
     :                                  40/
c     data (SENIOR(10,1,2,i),i=1,11)   / 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
c    :                                   2/
      data (SENIOR(10,1,2,i),i=1,11)   /-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 
     :                                  -1/
*     l=10 #=2 +
      if (first) then
         fil = fil_1
      else
         fil = fil_2
      endif
      do 12 i=1,20
   12    antko(i) = 1
      pos = 0
      do 20 i=1,110
         n = posn(i)
         l = posl(i)
!Jacek mailed the fix 98-10-29
         do 20 k=0,min(L,1)
         !do 20 k=0,min(n-1,1)
            if (ansats(n,l,k).NE.0) then
               rad1(pos*9+1:pos*9+9)  = '         '
               if (n.LT.10) then
                  rad1(pos*9+3:pos*9+3) = char(48+n)
               else
                  n1 = mod(n,10)
                  n10 = n/10
                  rad1(pos*9+2:pos*9+2) = char(48+n10)
                  rad1(pos*9+3:pos*9+3) = char(48+n1)
               endif
               rad1(pos*9+4:pos*9+5) = L1(l,k)
               rad1(pos*9+6:pos*9+9) = '(  )'
               if (ansats(n,l,k).GE.10) then
                  rad1(pos*9+7:pos*9+8) = char(ansats(n,l,k)/10+48)
               else
                  rad1(pos*9+7:pos*9+7) = ' '
               endif
               rad1(pos*9+8:pos*9+8) = char(mod(ansats(n,l,k),10)+48)
               pos = pos + 1
               if (pos.GT.skal) then
		  write(*,*) 'More than 20 subshells'
		  return
	       endif
               orbit(pos) = l
               antel(pos) = min(ansats(n,l,k),
     :                                        antmax(l,k)-ansats(n,l,k))
               antko(pos) = koppl(l,antel(pos),k)
               plus(pos) = k
            endif
   20 continue

      if (pos.EQ.0) return
      do 900 i1=1,antko(1)
       do 900 i2=1,antko(2)
        do 900 i3=1,antko(3)
         do 900 i4=1,antko(4)
          do 900 i5=1,antko(5)
           do 900 i6=1,antko(6)
            do 900 i7=1,antko(7)
             do 900 i8=1,antko(8)
              do 900 i9=1,antko(9)
               do 900 i10=1,antko(10)
                do 900 i11=1,antko(11)
                 do 900 i12=1,antko(12)
                  do 900 i13=1,antko(13)
                   do 900 i14=1,antko(14)
                    do 900 i15=1,antko(15)
                     do 900 i16=1,antko(16)
                      do 900 i17=1,antko(17)
                       do 900 i18=1,antko(18)
                        do 900 i19=1,antko(19)
                         do 900 i20=1,antko(20)

                          J(1) = JKVANT(orbit(1),plus(1),antel(1),i1)
                          S(1) = SENIOR(orbit(1),plus(1),antel(1),i1)
                          if (pos.EQ.1) then
                           if (J(1).GE.minJ .AND. J(1).LE.maxJ) then
                            call Kopp1(pos,rad2,J,S,antko)
                            call Kopp2(pos,rad3,J,J,par,antko)
                            write(fil,999) rad1(1:9)
                            write(fil,999) rad2(1:9)
                            write(fil,999) rad3(1:11)
                            cf = cf + 1
                           endif
                          else                                          pos>1 
                          
                          do 25 resJ=minJ,maxJ,2                      
                           JK(pos-1) = resJ
                           J(2) = JKVANT(orbit(2),plus(2),antel(2),i2)
                           S(2) = SENIOR(orbit(2),plus(2),antel(2),i2)
                           if (pos.EQ.2) then
                            if (resJ.GE.abs(J(1)-J(2)) .AND. 
     :                                   resJ.LE.J(1)+J(2)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:18)
                             write(fil,999) rad2(1:18)
                             write(fil,999) rad3(1:20)
                             cf = cf + 1
                            endif
                           else                                         pos>2

                          J(3) = JKVANT(orbit(3),plus(3),antel(3),i3)
                          S(3) = SENIOR(orbit(3),plus(3),antel(3),i3)
                          do 30 JK1=abs(J(1)-J(2)),J(1)+J(2),2
                           JK(1) = JK1
                           if (pos.EQ.3) then
                            if (resJ.GE.abs(JK1-J(3)) .AND. 
     :                                        resJ.LE.JK1+J(3)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:27)
                             write(fil,999) rad2(1:27)
                             write(fil,999) rad3(1:29)
                             cf = cf + 1
                            endif
                           else                                         pos>3

                          J(4) = JKVANT(orbit(4),plus(4),antel(4),i4)
                          S(4) = SENIOR(orbit(4),plus(4),antel(4),i4)
                          do 40 JK2=abs(JK1-J(3)),JK1+J(3),2
                           JK(2) = JK2
                           if (pos.EQ.4) then
                            if (resJ.GE.abs(JK2-J(4)) .AND. 
     :                                     resJ.LE.JK2+J(4)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:36)
                             write(fil,999) rad2(1:36)
                             write(fil,999) rad3(1:38)
                             cf = cf + 1
                            endif
                           else                                         pos>4

                          J(5) = JKVANT(orbit(5),plus(5),antel(5),i5)
                          S(5) = SENIOR(orbit(5),plus(5),antel(5),i5)
                          do 50 JK3=abs(JK2-J(4)),JK2+J(4),2
                           JK(3) = JK3
                           if(pos.EQ.5) then
                            if (resJ.GE.abs(JK3-J(5)) .AND. 
     :                                     resJ.LE.JK3+J(5)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:45)
                             write(fil,999) rad2(1:45)
                             write(fil,999) rad3(1:47)
                             cf = cf + 1
                            endif
                           else                                         pos>5

                          J(6) = JKVANT(orbit(6),plus(6),antel(6),i6)
                          S(6) = SENIOR(orbit(6),plus(6),antel(6),i6)
                          do 60 JK4=abs(JK3-J(5)),JK3+J(5),2
                           JK(4) = JK4
                           if(pos.EQ.6) then
                            if (resJ.GE.abs(JK4-J(6)) .AND.
     :                                      resJ.LE.JK4+J(6)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:54)
                             write(fil,999) rad2(1:54)
                             write(fil,999) rad3(1:56)
                             cf = cf + 1
                            endif
                           else                                         pos>6

                          J(7) = JKVANT(orbit(7),plus(7),antel(7),i7)
                          S(7) = SENIOR(orbit(7),plus(7),antel(7),i7)
                          do 70 JK5=abs(JK4-J(6)),JK4+J(6),2
                           JK(5) = JK5
                           if(pos.EQ.7) then
                            if (resJ.GE.abs(JK5-J(7)) .AND.
     :                                      resJ.LE.JK5+J(7)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:63)
                             write(fil,999) rad2(1:63)
                             write(fil,999) rad3(1:65)
                             cf = cf + 1
                            endif
                           else                                         pos>7
 
                          J(8) = JKVANT(orbit(8),plus(8),antel(8),i8)
                          S(8) = SENIOR(orbit(8),plus(8),antel(8),i8)
                          do 80 JK6=abs(JK5-J(7)),JK5+J(7),2
                           JK(6) = JK6
                           if(pos.EQ.8) then
                            if (resJ.GE.abs(JK6-J(8)) .AND.
     :                                    resJ.LE.JK6+J(8)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:72)
                             write(fil,999) rad2(1:72)
                             write(fil,999) rad3(1:74)
                             cf = cf + 1
                            endif
                           else                                         pos>8
 
                          J(9) = JKVANT(orbit(9),plus(9),antel(9),i9)
                          S(9) = SENIOR(orbit(9),plus(9),antel(9),i9)
                          do 90 JK7=abs(JK6-J(8)),JK6+J(8),2
                           JK(7) = JK7
                           if(pos.EQ.9) then
                            if (resJ.GE.abs(JK7-J(9)) .AND.
     :                                   resJ.LE.JK7+J(9)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:81)
                             write(fil,999) rad2(1:81)
                             write(fil,999) rad3(1:83)
                             cf = cf + 1
                            endif
                           else                                         pos>9

                          J(10) = 
     :                          JKVANT(orbit(10),plus(10),antel(10),i10)
                          S(10) = 
     :                          SENIOR(orbit(10),plus(10),antel(10),i10)
                          do 100 JK8=abs(JK7-J(9)),JK7+J(9),2
                           JK(8) = JK8
                           if(pos.EQ.10) then
                            if (resJ.GE.abs(JK8-J(10)) .AND.
     :                                   resJ.LE.JK8+J(10)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:90)
                             write(fil,999) rad2(1:90)
                             write(fil,999) rad3(1:92)
                             cf = cf + 1
                            endif
                           else                                         pos>10

                          J(11) = 
     :                          JKVANT(orbit(11),plus(11),antel(11),i11)
                          S(11) = 
     :                          SENIOR(orbit(11),plus(11),antel(11),i11)
                          do 110 JK9=abs(JK8-J(10)),JK8+J(10),2
                           JK(9) = JK9
                           if(pos.EQ.11) then
                            if (resJ.GE.abs(JK9-J(11)) .AND.
     :                                    resJ.LE.JK9+J(11)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:99)
                             write(fil,999) rad2(1:99)
                             write(fil,999) rad3(1:101)
                             cf = cf + 1
                            endif
                           else                                         pos>11
 
                          J(12) = 
     :                          JKVANT(orbit(12),plus(12),antel(12),i12)
                          S(12) = 
     :                          SENIOR(orbit(12),plus(12),antel(12),i12)
                          do 120 JK10=abs(JK9-J(11)),JK9+J(11),2
                           JK(10) = JK10
                           if(pos.EQ.12) then
                            if (resJ.GE.abs(JK10-J(12)) .AND.
     :                                   resJ.LE.JK10+J(12)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:108)
                             write(fil,999) rad2(1:108)
                             write(fil,999) rad3(1:110)
                             cf = cf + 1
                            endif
                           else                                         pos>12

                          J(13) = 
     :                          JKVANT(orbit(13),plus(13),antel(13),i13)
                          S(13) = 
     :                          SENIOR(orbit(13),plus(13),antel(13),i13)
                          do 130 JK11=abs(JK10-J(12)),JK10+J(12),2
                           JK(11) = JK11
                           if(pos.EQ.13) then
                            if (resJ.GE.abs(JK11-J(13)) .AND.
     :                                   resJ.LE.JK11+J(13)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:117)
                             write(fil,999) rad2(1:117)
                             write(fil,999) rad3(1:119)
                             cf = cf + 1
                            endif
                           else                                         pos>13

                          J(14) = 
     :                          JKVANT(orbit(14),plus(14),antel(14),i14)
                          S(14) = 
     :                          SENIOR(orbit(14),plus(14),antel(14),i14)
                          do 140 JK12=abs(JK11-J(13)),JK11+J(13),2
                           JK(12) = JK12
                           if(pos.EQ.14) then
                            if (resJ.GE.abs(JK12-J(14)) .AND.
     :                                   resJ.LE.JK12+J(14)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:126)
                             write(fil,999) rad2(1:126)
                             write(fil,999) rad3(1:128)
                             cf = cf + 1
                            endif
                           else                                         pos>14

                          J(15) = 
     :                          JKVANT(orbit(15),plus(15),antel(15),i15)
                          S(15) = 
     :                          SENIOR(orbit(15),plus(15),antel(15),i15)
                          do 150 JK13=abs(JK12-J(14)),JK12+J(14),2
                           JK(13) = JK13
                           if(pos.EQ.15) then
                            if (resJ.GE.abs(JK13-J(15)) .AND.
     :                                   resJ.LE.JK13+J(15)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:135)
                             write(fil,999) rad2(1:135)
                             write(fil,999) rad3(1:137)
                             cf = cf + 1
                            endif
                           else                                         pos>15
 
                          J(16) = 
     :                          JKVANT(orbit(16),plus(16),antel(16),i16)
                          S(16) = 
     :                          SENIOR(orbit(16),plus(16),antel(16),i16)
                          do 160 JK14=abs(JK13-J(15)),JK13+J(15),2
                           JK(14) = JK14
                           if(pos.EQ.16) then
                            if (resJ.GE.abs(JK14-J(16)) .AND.
     :                                   resJ.LE.JK14+J(16)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:144)
                             write(fil,999) rad2(1:144)
                             write(fil,999) rad3(1:146)
                             cf = cf + 1
                            endif
                           else                                         pos>16

                          J(17) = 
     :                          JKVANT(orbit(17),plus(17),antel(17),i17)
                          S(17) = 
     :                          SENIOR(orbit(17),plus(17),antel(17),i17)
                          do 170 JK15=abs(JK14-J(16)),JK14+J(16),2
                           JK(15) = JK15
                           if(pos.EQ.17) then
                            if (resJ.GE.abs(JK15-J(17)) .AND.
     :                                     resJ.LE.JK15+J(17)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:153)
                             write(fil,999) rad2(1:153)
                             write(fil,999) rad3(1:155)
                             cf = cf + 1
                            endif
                           else                                             pos>17
 
                          J(18) = 
     :                          JKVANT(orbit(18),plus(18),antel(18),i18)
                          S(18) = 
     :                          SENIOR(orbit(18),plus(18),antel(18),i18)
                          do 180 JK16=abs(JK15-J(17)),JK15+J(17),2
                           JK(16) = JK16
                           if(pos.EQ.18) then
                            if (resJ.GE.abs(JK16-J(18)) .AND.
     :                                   resJ.LE.JK16+J(18)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:162)
                             write(fil,999) rad2(1:162)
                             write(fil,999) rad3(1:164)
                             cf = cf + 1
                            endif
                           else                                         pos>18
 
                          J(19) = 
     :                          JKVANT(orbit(19),plus(19),antel(19),i19)
                          S(19) = 
     :                          SENIOR(orbit(19),plus(19),antel(19),i19)
                          do 190 JK17=abs(JK16-J(18)),JK16+J(18),2
                           JK(17) = JK17
                           if(pos.EQ.19) then
                            if (resJ.GE.abs(JK17-J(19)) .AND.
     :                                   resJ.LE.JK17+J(19)) then
                             call Kopp1(pos,rad2,J,S,antko)
                             call Kopp2(pos,rad3,JK,J,par,antko)
                             write(fil,999) rad1(1:171)
                             write(fil,999) rad2(1:171)
                             write(fil,999) rad3(1:173)
                             cf = cf + 1
                            endif
                           else                                         pos=last
 
                          J(20) = 
     :                          JKVANT(orbit(20),plus(20),antel(20),i20)
                          S(20) = 
     :                          SENIOR(orbit(20),plus(20),antel(20),i20)
                          do 200 JK18=abs(JK17-J(19)),JK17+J(19),2
                           if (resJ.GE.abs(JK18-J(20)) .AND.
     :                                          resJ.LE.JK18+J(20)) then
                            JK(18) = JK18
                            call Kopp1(pos,rad2,J,S,antko)
                            call Kopp2(pos,rad3,JK,J,par,antko)
                            write(fil,999) rad1(1:180)
                            write(fil,999) rad2(1:180)
                            write(fil,999) rad3(1:182)
                            cf = cf + 1
                           endif
  200                     continue
                          endif                                         pos=19?
  190                     continue
                          endif                                         pos=18?
  180                     continue
                          endif                                         pos=17?
  170                     continue
                          endif                                         pos=16?
  160                     continue
                          endif                                         pos=15?
  150                     continue
                          endif                                         pos=14?
  140                     continue
                          endif                                         pos=13?
  130                     continue
                          endif                                         pos=12?
  120                     continue
                          endif                                         pos=11?
  110                     continue
                          endif                                         pos=10?
  100                     continue
                          endif                                         pos=9?
   90                     continue
                          endif                                         pos=8?
   80                     continue
                          endif                                         pos=7?
   70                     continue
                          endif                                         pos=6?
   60                     continue
                          endif                                         pos=5?
   50                     continue
                          endif                                         pos=4?
   40                     continue
                          endif                                         pos=3?
   30                     continue
                          endif                                         pos=2?
   25                     continue
                          endif                                         pos=1?
  900         continue
  999 format(2A)
      return
      end
