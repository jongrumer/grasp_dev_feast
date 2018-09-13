
*     last edited September 23, 1995
      logical function Lika(pop0,pop1)
      integer i,j,pop0(15,0:10,0:1),pop1(15,0:10,0:1),k
      logical dum
      dum = .TRUE.
      do 10 i=1,15
         do 10 j=0,min(10,i-1)
            do 10 k=0,1
               dum = dum .AND. pop0(i,j,k).EQ.pop1(i,j,k)
   10          if (.NOT.dum) goto 20
   20 Lika = dum
      return
      end
