*     last edited July 30, 1996
      subroutine Fivelines(org,locked,closed,first,posn,posl)
      integer org(1:15,0:10),posn(110),posl(110)
      logical locked(1:15,0:10),closed(1:15,0:10),first
      integer k,i,j,start,stopp
      character rad*1000,orb(0:10)
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/

      write(*,*) 'In fivelines'
      if (.NOT.first) then
C        open(unit=8,file='fil2.dat',status='unknown')
         open(unit=8,status='scratch')
         write(8,999)
         write(8,999)
         write(8,999)
      else
         open(unit=7,file='fil1.dat',status='unknown')
c        open(unit=7,status='scratch')
         write(7,999) 'Core subshells:'
         do 10 i=1,1000
   10       rad(i:i) = ' '
         start =-2
         stopp = 0
         do 20 k=1,110
            i = posn(k)
            j = posl(k)
            if (closed(i,j)) then
               start                = start + 5
               stopp                = stopp + 5
               if (i.LT.10) then
                  rad(start:start)     = CHAR(ICHAR('0')+i)
               else
                  rad(start-1:start-1) = '1'
                  rad(start:start)     = CHAR(ICHAR('0')+i-10)
               endif
               rad(start+1:start+1) = orb(j)
               org(i,j)             = 0
               if (j.GE.1) then
                  rad(start+2:start+2) = '-'
                  start                = start + 5
                  stopp                = stopp + 5
                  if (i.LT.10) then
                     rad(start:start)     = CHAR(ICHAR('0')+i)
                  else
                     rad(start-1:start-1) = '1'
                     rad(start:start)     = CHAR(ICHAR('0')+i-10)
                  endif
                  rad(start+1:start+1) = orb(j)
               endif
            endif
   20    continue
         if (stopp .EQ. 0) then
            write(7,999)
         else
            write(7,999) rad(1:stopp)
         endif
         write(7,999) 'Peel subshells:'
      endif
      do 30 i=1,1000
   30    rad(i:i) = ' '
      start =-2
      stopp = 0
      do 40 k=1,110
         i=posn(k)
         j=posl(k)
         if ((.NOT.(org(i,j).EQ.0 .AND. locked(i,j))) .AND. 
     :                              (.NOT.closed(i,j))) then
            start                = start + 5
            stopp                = stopp + 5
            if (i.LT.10) then
               rad(start:start)     = CHAR(ICHAR('0')+i)
            else
               rad(start-1:start-1) = '1'
               rad(start:start)     = CHAR(ICHAR('0')+i-10)
            endif
            rad(start+1:start+1) = orb(j)
            if (j.GE.1) then
               rad(start+2:start+2) = '-'
               start                = start + 5
               stopp                = stopp + 5
               if (i.LT.10) then
                  rad(start:start)     = CHAR(ICHAR('0')+i)
               else
                  rad(start-1:start-1) = '1'
                  rad(start:start)     = CHAR(ICHAR('0')+i-10)
               endif
               rad(start+1:start+1) = orb(j)
            endif
C           write(*,*) i,rad(1:100)    
         endif
   40 continue
      if (first) then
         if (stopp .EQ. 0) then
            write(7,999)
         else
            write(7,999) rad(1:stopp)
         endif
         write(7,999) 'CSF(s):'
      else
         if (stopp .EQ. 0) then
            write(8,999)
         else
            write(8,999) rad(1:stopp)
         endif
         write(8,999) 'CSF(s):'
      endif
  999 format(A)
      return
      end
