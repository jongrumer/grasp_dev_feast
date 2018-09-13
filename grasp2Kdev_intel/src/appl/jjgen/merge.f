*     last edited July 31, 1996
      subroutine Merge(single,posn,posl,ii)
      integer fil_1,fil_2,utfil,nyfil,posn(110),posl(110)
      parameter (fil_1=7, fil_2=8, utfil=9, nyfil=13)
      character rad11*200,rad12*200,rad21*200,rad22*200
      character rad31*200,rad32*200
      logical p1,p2,slut1,slut2,Lika,single
      integer pop1(15,0:10,0:1),pop2(15,0:10,0:1),skal1,skal2
      integer i,j,k,cf,stopp1,stopp2,popo(15,0:10,0:1)

      ii1=mod(ii,2)
      if (ii1.eq.0) then
        open(unit=utfil,file='rcsf.out',status='unknown')
      else
        open(unit=utfil,file='fil1.dat',status='unknown')
      endif
      open(unit=nyfil,file='clist.new',status='unknown')
      slut1 = .FALSE.
      slut2 = single
      cf    = 0
      call Fivefirst(slut1,slut2,posn,posl)
      skal1 = 20
      skal2 = 20
      call Lasa1(fil_1,rad11,pop1,skal1,slut1)
      call Lasa1(fil_2,rad12,pop2,skal2,slut2)
   10 if (.NOT.slut1 .AND. .NOT.slut2) then
         call Test(p1,p2,pop1,pop2,15)  
         if (p1) then
            do 20 i=1,15
               do 20 j=0,min(10,i-1)
                  do 20 k=0,1
   20                popo(i,j,k) = pop1(i,j,k)
            stopp1 = max(1,9*skal1)
            stopp2 = 9*skal1 + 2
   30       call Lasa2(fil_1,rad21,rad31,stopp1,slut1)
            if (.NOT.slut1) then
               write(utfil,999) rad11(1:stopp1)
C               write(utfil,999) 'Hej1'
               write(utfil,999) rad21(1:stopp1)
               write(utfil,999) rad31(1:stopp2)
               cf = cf + 1
            endif
            skal1 = 20
            call Lasa1(fil_1,rad11,pop1,skal1,slut1)
            if (.NOT.slut1) then
               if (Lika(popo,pop1)) goto 30
            endif
            if (p2) then
   40          call Lasa2(fil_2,rad22,rad32,stopp1,slut2)
               skal2 = 20
               call Lasa1(fil_2,rad12,pop2,skal2,slut2)
               if (.NOT.slut2) then
                  if (Lika(popo,pop2)) goto 40
               endif
            endif
            goto  10
         elseif (p2) then
            do 50 i=1,15
               do 50 j=0,min(10,i-1)
                  do 50 k=0,1
   50                popo(i,j,k) = pop2(i,j,k)
            stopp1 = max(1,9*skal2)
            stopp2 = 9*skal2 + 2
   60       call Lasa2(fil_2,rad22,rad32,stopp1,slut2)
            if (.NOT.slut2) then
               write(utfil,999) rad12(1:stopp1)
C               write(utfil,999) 'Hej2'
               write(utfil,999) rad22(1:stopp1)
               write(utfil,999) rad32(1:stopp2)
               write(nyfil,999) rad12(1:stopp1)
C               write(utfil,999) 'Hej3'
               write(nyfil,999) rad22(1:stopp1)
               write(nyfil,999) rad32(1:stopp2)
               cf = cf + 1
            endif
            skal2 = 20
            call Lasa1(fil_2,rad12,pop2,skal2,slut2)
            if (.NOT.slut2) then
               if (Lika(popo,pop2)) goto 60
            endif
            goto 10
         else
            write(*,*) 'fatal error'
            stop
         endif
      elseif (.NOT.slut1 .AND. slut2) then
   70    stopp1 = max(1,9*skal1)
         stopp2 = 9*skal1 + 2
         call Lasa2(fil_1,rad21,rad31,stopp1,slut1)
C         write(*,*) 'file1',fil_1
         if (.NOT.slut1) then
            write(utfil,999) rad11(1:stopp1)
C               write(utfil,999) 'Hej4'
            write(utfil,999) rad21(1:stopp1)
            write(utfil,999) rad31(1:stopp2)
            cf = cf + 1
         endif
         skal1 = 20
         call Lasa1(fil_1,rad11,pop1,skal1,slut1)
         if (.NOT.slut1) goto 70
      elseif (slut1 .AND. .NOT.slut2) then
   80    stopp1 = max(1,9*skal2)
         stopp2 = 9*skal2 + 2
         call Lasa2(fil_2,rad22,rad32,stopp1,slut2)
         if (.NOT.slut2) then
            write(utfil,999) rad12(1:stopp1)
C               write(utfil,999) 'Hej5'
            write(utfil,999) rad22(1:stopp1)
            write(utfil,999) rad32(1:stopp2)
            write(nyfil,999) rad12(1:stopp1)
C               write(utfil,999) 'Hej6'
            write(nyfil,999) rad22(1:stopp1)
            write(nyfil,999) rad32(1:stopp2)
            cf = cf + 1
         endif
         skal2 = 20
         call Lasa1(fil_2,rad12,pop2,skal2,slut2)
         if (.NOT.slut2) goto 80
      endif
      close(fil_1)
      close(fil_2)
      close(utfil)
      close(nyfil)
      if (cf.EQ.0) then
         write(*,105) 'No configuration state in the final list.'
      elseif (cf.EQ.1) then
         write(*,105) 'One configuration state in the final list.'
      elseif (cf.LT.10) then
         write(*,101) cf,' configuration states in the final list.'
      elseif (cf.LT.100) then
         write(*,102) cf,' configuration states in the final list.'
      elseif (cf.LT.1000) then
         write(*,103) cf,' configuration states in the final list.'
      elseif (cf.LT.10000) then
         write(*,104) cf,' configuration states in the final list.'
      elseif (cf.LT.100000) then
         write(*,106) cf,' configuration states in the final list.'
      else
         write(*,*) cf,' configuration states in the final list.'
      endif
      return
  101 format(' ',I1,A)
  102 format(' ',I2,A)
  103 format(' ',I3,A)
  104 format(' ',I4,A)
  105 format(' ',A)
  106 format(' ',I5,A)
  999 format(A)
      end
