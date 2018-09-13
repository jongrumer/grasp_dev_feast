*     last edited Februar 20, 1996
      subroutine Lockad(closed,med,slut,expand)
      integer fil_1,fil_2
      parameter (fil_1=7, fil_2=8)
      logical closed(15,0:10),med(15,0:10),slut,expand
      character rad*1000,orb(0:10)
      integer i,j,n,l
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      if (expand) then
         read(fil_2,*,end=40)
         read(fil_2,100,end=40) rad
      else
         read(fil_1,*,end=40)
         read(fil_1,100,end=40) rad
      endif
      do 10 n=1,15
         do 10 l=1,min(10,n-1)
   10       closed(n,l) = .FALSE.
      do 30 i=0,205
         j = i*5
         n = ichar(rad(j+3:j+3)) - ichar('0')
         if (n.GE.1 .AND. n.LE.15) then
            do 20 l=0,min(10,n-1)
               if (rad(j+4:j+4).EQ.orb(l)) then
                  closed(n,l) = .TRUE.
                  goto 30
               endif
   20       continue
         else
            goto 31
         endif
   30 continue
   31 continue
      if (expand) then
         read(fil_2,*,end=40)
         read(fil_2,100,end=40) rad
      else
         read(fil_1,*,end=40)
         read(fil_1,100,end=40) rad
      endif
      do 110 n=1,15
         do 110 l=1,min(10,n-1)
  110       med(n,l) = .FALSE.
      do 130 i=0,205
         j = i*5
         n = ichar(rad(j+3:j+3)) - ichar('0')
         if (n.GE.1 .AND. n.LE.15) then
            do 120 l=0,min(10,n-1)
               if (rad(j+4:j+4).EQ.orb(l)) then
                  med(n,l) = .TRUE.
                  goto 130
               endif
  120       continue
         else
            return
         endif
  130 continue

      return
   40 slut = .TRUE.
      return
  100 format(A)
      end
