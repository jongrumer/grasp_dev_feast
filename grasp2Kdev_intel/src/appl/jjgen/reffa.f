*     last edited July 31, 1996
      subroutine Reffa(posn,posl)
      integer logfil,reffil
      parameter (logfil=31, reffil=18)
      integer posn(110),posl(110),stat(15,0:10),i,n,l,num
      character orb(0:10)*1,M*1,X*1,Y*1
      logical OK
      data (orb(i),i=0,10) /'s','p','d','f','g','h','i','k','l','m','n'/

      do 10 n=1,15
         do 10 l=0,min(n-1,10)
   10       stat(n,l) = 0
      write(*,200)
     :       'Default, reverse, symmetry or user specified ordering?',
     :       ' (*/r/s/u)'
      read(*,1000) X
      if (X.EQ.'u' .OR. X.EQ.'U') then
         write(logfil,*) 'User specified ordering.'
         inquire(file='clist.ref',exist=OK)
         if (OK) open(unit=reffil,status='old',file='clist.ref')
         l   = -1
         num = 1
         if (.NOT.OK) then
            write(*,200) 'No reference file present! ',
     :              'The couplings will appear in standard order.'
         else
            write(*,200) 'Reference file present!'
   20       read(reffil,1000,end=40) M,X,Y
               n = ichar(M) - ichar('0')
               if (X.GE.'0' .AND. X.LE.'9') then
                  n = n*10 + ichar(X) - ichar('0')
                  X = Y
               endif
               do 30 i=0,10
                  if (orb(i).EQ.X) l=i
   30          continue
               if (l.EQ.-1 .OR. n.LT.0 .OR. n.GT.15 .OR.
     :                                      n.LE.l .OR. l.GT.10) goto 40
               if (stat(n,l).NE.0) then
                  write(*,200)
     :                 'The same orbital appeared more than once!'
                  l = -1
                  goto 20
               endif
               posn(num) = n
               posl(num) = l
               stat(n,l) = num
               num       = num + 1
               l         = -1
               goto 20
   40       continue
            if (num.EQ.1) then
               write(*,200) 'The program failed reading the order of ',
     :                 'the customized coupling scheme.'
            else
               write(*,200) 'The couplings will ',
     :                      'be made in the following customized order:'
               if (num.EQ.2) then
                  write(*,100) posn(1),orb(posl(1))
               else
                  write(*,100) posn(1),orb(posl(1)),
     :                      (',',posn(i),orb(posl(i)),i=2,num-1)
               endif
            endif
         endif
         do 50 n=1,15
            do 50 l=0,min(n-1,10)
               if (stat(n,l).EQ.0) then
                  posn(num) = n
                  posl(num) = l
                  num       = num + 1
               endif
   50    continue
         close(reffil)
         write(*,200)
      elseif (X.EQ.'s' .OR. X.EQ.'S') then
         write(logfil,*) 'Symmetry ordering.'
         num = 1
         do 60 l=0,10
            do 60 n=l+1,15
               posn(num) = n
               posl(num) = l
   60          num       = num + 1
      elseif (X.EQ.'r' .OR. X.EQ.'R') then
         write(logfil,*) 'Reverse ordering.'
         num = 1
         do 80 n=15,1,-1
            do 80 l=min(n-1,10),0,-1
               posn(num) = n
               posl(num) = l
   80          num       = num + 1
      else
         write(logfil,*) 'Standard ordering.'
         num = 1
         do 70 n=1,15
            do 70 l=0,min(n-1,10)
               posn(num) = n
               posl(num) = l
   70          num       = num + 1
      endif
      return
  100 format(' ',110(I2,A,A))
  200 format(' ',2A)
 1000 format(3A)
      end
