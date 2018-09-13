*     last edited July 31, 1996
      subroutine Mergeb(antal)
      logical p1,p2,slut1,slut2
      integer pop1(15,0:10,0:1),pop2(15,0:10,0:1),popo(15,0:10,0:1)
      integer i,j,k,antal
      slut1 = .FALSE.
      slut2 = .FALSE.
      antal = 0
      open(unit=22,status='scratch')
      do 1 i=1,15
         read(20,5000,end=2) (pop1(i,j,0),j=0,min(10,i-1))
    1    read(20,5000,end=2) (pop1(i,j,1),j=0,min(10,i-1))
      goto 3
    2 slut1 = .TRUE.
    3 do 4 i=1,15
         read(21,5000,end=5) (pop2(i,j,0),j=0,min(10,i-1))
    4    read(21,5000,end=5) (pop2(i,j,1),j=0,min(10,i-1))
      goto 6
    5 slut2 = .TRUE.
    6 continue
   10 if (.NOT.slut1 .AND. .NOT.slut2) then
         call Test(p1,p2,pop1,pop2,15)
         if (p1) then
            do 20 i=1,15
               do 20 j=0,min(10,i-1)
                  do 20 k=0,1
   20                popo(i,j,k) = pop1(i,j,k)
            do 120 i=1,15
               write(22,5000) (pop1(i,j,0),j=0,min(10,i-1))
  120          write(22,5000) (pop1(i,j,1),j=0,min(10,i-1))
            do 121 i=1,15
               read(20,5000,end=21) (pop1(i,j,0),j=0,min(10,i-1))
  121          read(20,5000,end=21) (pop1(i,j,1),j=0,min(10,i-1))  
            goto 22
   21       slut1 = .TRUE.
   22       continue
            if (p2) then
            do 122 i=1,15
               read(21,5000,end=23) (pop2(i,j,0),j=0,min(10,i-1))
  122          read(21,5000,end=23) (pop2(i,j,1),j=0,min(10,i-1))
               goto 10
   23          slut2 = .TRUE.
            endif
         elseif (p2) then
            do 50 i=1,15
               do 50 j=0,min(10,i-1)         
                  do 50 k=0,1
   50                popo(i,j,k) = pop2(i,j,k)
            if (.NOT.slut2) then
            do 51 i=1,15
               write(22,5000) (pop2(i,j,0),j=0,min(10,i-1))
   51          write(22,5000) (pop2(i,j,1),j=0,min(10,i-1))
            endif
            do 52 i=1,15
               read(21,5000,end=53) (pop2(i,j,0),j=0,min(10,i-1))
   52          read(21,5000,end=53) (pop2(i,j,1),j=0,min(10,i-1))
            goto 10
   53       slut2 = .TRUE.
         else
            write(*,*) 'fatal error'
            stop
         endif
         goto 10
      elseif (.NOT.slut1 .AND. slut2) then
   70    continue
         do 170 i=1,15
           write(22,5000) (pop1(i,j,0),j=0,min(10,i-1))
  170      write(22,5000) (pop1(i,j,1),j=0,min(10,i-1))
         do 171 i=1,15
            read(20,5000,end=71) (pop1(i,j,0),j=0,min(10,i-1))
  171       read(20,5000,end=71) (pop1(i,j,1),j=0,min(10,i-1))
         goto 70
   71    slut1 = .TRUE.
      elseif (slut1 .AND. .NOT.slut2) then
   80    continue
         do 180 i=1,15
            write(22,5000) (pop2(i,j,0),j=0,min(10,i-1))
  180       write(22,5000) (pop2(i,j,1),j=0,min(10,i-1))
         do 181 i=1,15
            read(21,5000,end=81) (pop2(i,j,0),j=0,min(10,i-1))
  181       read(21,5000,end=81) (pop2(i,j,1),j=0,min(10,i-1))
         goto 80
   81    slut2 = .TRUE.
      endif
      rewind(22)
      close(20)
      close(21)
      open(unit=20,status='scratch')
  580 continue
      do 581 i=1,15
         read(22,5000,end=999) (pop2(i,j,0),j=0,min(10,i-1))
  581    read(22,5000,end=999) (pop2(i,j,1),j=0,min(10,i-1))
      do 582 i=1,15
         write(20,5000) (pop2(i,j,0),j=0,min(10,i-1))
  582    write(20,5000) (pop2(i,j,1),j=0,min(10,i-1))
      antal = antal + 1
      goto 580
  999 close(22)
      rewind(20)
      return
 5000 format(11I2)
      end
