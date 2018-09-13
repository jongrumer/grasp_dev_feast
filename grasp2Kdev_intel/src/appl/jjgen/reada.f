*     last edited July 31, 1996
      subroutine Reada(rad1,pop,skal,slut)
      character rad1*200,orb(0:10)
      integer   skal,i,j,k,n,l,antal,pop(15,0:10,0:1),stopp
      logical   slut
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      slut = .false.
      do 5 n=1,15
         do 5 l=0,min(10,n-1)
            do 5 k=0,1
    5          pop(n,l,k) = 0
      stopp = skal-1
      do 10 i=0,stopp
         j = 9*i
         if (rad1(j+3:j+3).EQ.' ') return
         skal = i + 1
         slut = .TRUE.
         n = ichar(rad1(j+3:j+3)) - ichar('0')
         if (rad1(j+2:j+2).EQ.'1') n = n+10
         if (n.LE.15 .AND. n.GE.1) then
            do 15 l=0,min(10,n-1)
               if (rad1(j+4:j+4).EQ.orb(l)) then
                  slut = .FALSE.
                  if (rad1(j+7:j+7).EQ.' ' .OR. rad1(j+7:j+7).EQ.'0')
     :                                                             then
                     antal = 0
                  else
                     antal = 10*(ichar(rad1(j+7:j+7)) - ichar('0'))
                  endif
                  antal = antal + ichar(rad1(j+8:j+8)) - ichar('0')
                  if (antal.GT.4*l+2) then
                     slut = .TRUE.
                     return
                  endif
                  if (rad1(j+5:j+5).EQ.'-' .or. l.EQ.0) then
                     pop(n,l,0) = antal
                  else
                     pop(n,l,1) = antal
                  endif
                  goto 10
               endif
   15       continue
         else
            slut = .TRUE.
            return
         endif
   10 continue
      return
      end
