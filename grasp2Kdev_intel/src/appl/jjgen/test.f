* last edited October 31, 1996
      subroutine Test(p1,p2,pop1,pop2,nmax)
      logical p1,p2
      integer pop1(15,0:10,0:1),pop2(15,0:10,0:1),n,l,nmax,k

      p1 = .TRUE.
      p2 = .TRUE.
      do 10 n=1,nmax
         do 10 l=0,min(10,n-1)
            if ((pop1(n,l,1)+pop1(n,l,0)) .LT.
     :             (pop2(n,l,1)+pop2(n,l,0))) then
               p1 = .false.
               return
            elseif ((pop1(n,l,1)+pop1(n,l,0)) .GT. 
     :              (pop2(n,l,1)+pop2(n,l,0))) then
               p2 = .false.
               return
            elseif (pop1(n,l,1) .LT. pop2(n,l,1)) then
               p1 = .false.
               return
            elseif (pop1(n,l,1) .GT. pop2(n,l,1)) then
               p2 = .false.
               return
            elseif (pop1(n,l,0) .LT. pop2(n,l,0)) then
               p1 = .false.
               return
            elseif (pop1(n,l,0) .GT. pop2(n,l,0)) then
               p2 = .false.
               return
            endif
   10 continue
      return
      end
