*     last edited Februar 20, 1996
      subroutine Fivefirst(slut1,slut2,posn,posl)
      character rad0*1000,rad1*1000,rad2*1000
      integer k,i,j,n,l,pos,stopp,posn(110),posl(110)
      logical slut1,slut2,med(1:15,0:10,0:1)
      character orb(0:10)
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/

      read(7,999,end=100)
      write(9,999) 'Core subshells:'
      read(7,999,end=101) rad0
      stopp = 0
      do 1 i=0,210
         if (ichar(rad0(stopp+3:stopp+3)).GE.ichar('0') .AND.
     :       ichar(rad0(stopp+3:stopp+3)).LE.ichar('9') ) then
            stopp = stopp + 5
         else
            goto 2
         endif
    1 continue
    2 if (stopp.NE.0) then
         write(9,999) rad0(1:stopp)
      else
         write(9,999)
      endif
      write(9,999) 'Peel subshells:'
      read(7,999,end=102)
      read(7,999,end=102) rad1
      read(7,999,end=102)
      if (.NOT. slut2) then
         read(8,999,end=90)
         read(8,999,end=90)
         read(8,999,end=90)
         read(8,999,end=90) rad2
         read(8,999,end=90)
      endif
      if (.NOT.slut2) then
         do 10 i=1,15
            do 10 j=0,min(10,i-1)
               med(i,j,0) = .FALSE.
   10          med(i,j,1) = .FALSE.
         do 30 i=1,205
            pos = 5*i
            n = ichar(rad1(pos-2:pos-2)) - ichar('0')
            if (rad1(pos-3:pos-3).EQ.'1') n = n+10
            l = -1
            if (n.GE.1 .AND. n.LE.15) then
               do 20 j=0,min(10,n-1)
   20             if (rad1(pos-1:pos-1).EQ.orb(j)) l=j
            endif
            if (l.EQ.-1) goto 40
            if (rad1(pos:pos).EQ.'-' .OR. l.EQ.0) then
               med(n,l,0) = .TRUE.
            else
               med(n,l,1) = .TRUE.
            endif
   30    continue
   40    do 60 i=1,205
            pos = 5*i
            n = ichar(rad2(pos-2:pos-2)) - ichar('0')
            if (rad2(pos-3:pos-3).EQ.'1') n = n+10
            l = -1
            if (n.GE.1 .AND. n.LE.15) then
               do 50 j=0,min(10,n-1)
   50             if (rad2(pos-1:pos-1).EQ.orb(j)) l=j
            endif
            if (l.EQ.-1) goto 70
            if (rad2(pos:pos).EQ.'-' .OR. l.EQ.0) then
               med(n,l,0) = .TRUE.
            else
               med(n,l,1) = .TRUE.
            endif
   60    continue
   70    pos = 3
         do 80 k=1,110
            i = posn(k)
            j = posl(k)
            if (med(i,j,0)) then
               rad0(pos-2:pos+2) = '     '
               if (i.LT.10) then
                  rad0(pos:pos) = char(i + ichar('0'))
               else
                  rad0(pos:pos) = char(i + ichar('0') - 10)
                  rad0(pos-1:pos-1) = '1'
               endif
               rad0(pos+1:pos+1) = orb(j)
               if (j.NE.0) rad0(pos+2:pos+2) = '-'
               pos = pos + 5
            endif
            if (med(i,j,1)) then
               rad0(pos-2:pos+2) = '     '
               if (i.LT.10) then
                  rad0(pos:pos) = char(i + ichar('0'))
               else
                  rad0(pos:pos) = char(i + ichar('0') - 10)
                  rad0(pos-1:pos-1) = '1'
               endif
               rad0(pos+1:pos+1) = orb(j)
               pos = pos + 5
            endif
   80    continue
         write(9,999) rad0(1:pos-3)
         write(9,999) 'CSF(s):'
         return
      endif
   90 slut2 = .TRUE.
      stopp = 0
      do 91 i=0,210
         if (ichar(rad1(stopp+3:stopp+3)).GE.ichar('0') .AND.
     :       ichar(rad1(stopp+3:stopp+3)).LE.ichar('9') ) then
            stopp = stopp + 5
         else
            goto 92
         endif
   91 continue
   92 if (stopp.NE.0) then
CPJ         write(9,999) rad1(1:stopp)
         do i = 1,1000
           if (rad1(i:i).eq.':') rad1(i-1:i) = '10'
           if (rad1(i:i).eq.';') rad1(i-1:i) = '11'
           if (rad1(i:i).eq.'<') rad1(i-1:i) = '12'
           if (rad1(i:i).eq.'=') rad1(i-1:i) = '13'
           if (rad1(i:i).eq.'>') rad1(i-1:i) = '14'
           if (rad1(i:i).eq.'?') rad1(i-1:i) = '15'
         end do
         write(9,999) trim(rad1)
CPJ
      else
         write(9,999)
      endif
      write(9,999) 'CSF(s):'
      return
  100 slut1 = .TRUE.
      write(9,999) 'Core subshells:'
      read(8,999,end=200) rad0
  101 if (.NOT.slut1) then
         slut1 = .TRUE.
         read(8,999,end=200) rad0
      endif
      stopp = 0
      do 111 i=0,210
         if (ichar(rad0(stopp+3:stopp+3)).GE.ichar('0') .AND.
     :       ichar(rad0(stopp+3:stopp+3)).LE.ichar('9') ) then
            stopp = stopp + 5
         else
            goto 112
         endif
  111 continue
  112 if (stopp.NE.0) then
         write(9,999) rad0(1:stopp)
      else
         write(9,999)
      endif
      write(9,999) 'Peel subshells:'
  102 if (.NOT.slut1) then
         slut1 = .TRUE.
         read(8,999,end=200)
         read(8,999,end=200) rad2
      endif
      stopp = 0
      do 121 i=0,210
         if (ichar(rad2(stopp+3:stopp+3)).GE.ichar('0') .AND.
     :       ichar(rad2(stopp+3:stopp+3)).LE.ichar('9') ) then
            stopp = stopp + 5
         else
            goto 122
         endif
  121 continue
  122 if (stopp.NE.0) then
         write(9,999) rad2(1:stopp)
      else
         write(9,999)
      endif
      write(9,999) 'CSF(s):'
      return
  200 slut2 = .TRUE.
      return
  999 format(A)
      end
