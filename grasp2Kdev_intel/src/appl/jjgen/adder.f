*     last edited July 31, 1996
      subroutine Adder(closed,med,slut,anel,par,expand)
      integer fil_1,fil_2
      parameter (fil_1=7, fil_2=8)
      integer pop(1:15,0:10,0:1),skal,anel,par,i,j,kl,nr
      logical closed(1:15,0:10),med(1:15,0:10),slut,finns,expand
      character rad1*500,rad2*500,rad3*500
      skal = 20
      inquire(file='clist.inp',exist=finns)
      if (finns) then
         if (.NOT. expand) then
            open(unit=fil_1,file='clist.inp',status='old')
         else 
            open(unit=fil_2,file='clist.inp',status='old')
         endif
         slut = .FALSE.
         call Lockad(closed,med,slut,expand)
         if (.NOT.slut) then
            if (expand) then
               read(fil_2,*,end=99) 
               call Lasa1(fil_2,rad1,pop,skal,slut)
            else
               read(fil_1,*,end=99) 
               call Lasa1(fil_1,rad1,pop,skal,slut)
            endif
         endif
         if (.NOT.slut) then
            anel = 0
            par  = 0
            do 10 i=1,15
               do 10 j=0,min(10,i-1)
                  if (closed(i,j)) then
                     anel = anel + 2 + 4*j
                  else
                     anel = anel + pop(i,j,0) + pop(i,j,1)
                     par  = mod(par+j*(pop(i,j,0) + pop(i,j,1)),2)
                  endif
   10       continue
            if (expand) then
               read(fil_2,100,end=99) rad2
               read(fil_2,100,end=99) rad3
            else
               read(fil_1,100,end=99) rad2
               read(fil_1,100,end=99) rad3
            endif
            kl = skal*9
            if (rad3(kl:kl).NE.'/') then
               if (rad3(kl:kl).NE.' ') then
                  nr = 10*(ichar(rad3(kl:kl)) - ichar('0'))
               else
                  nr = 0
               endif
               kl = kl+1
               if (rad3(kl:kl).NE.' ') 
     :              nr = nr + (ichar(rad3(kl:kl)) - ichar('0'))
            else
               kl = skal*9-2
               if (rad3(kl:kl).NE.' ') then
                  nr = 10*(ichar(rad3(kl:kl)) - ichar('0'))
               else
                  nr = 0
               endif
               kl = kl+1
               nr = nr + ichar(rad3(kl:kl)) - ichar('0')
            endif 
         endif
         if (expand) then
            rewind(fil_2)
         else
            rewind(fil_1)
         endif
      else
         slut = .TRUE.
      endif
      return
   99 slut = .TRUE.
      if (expand) then
         close(fil_2)
      else
         close(fil_1)
      endif
      return
  100 format(A)
      end
