*     last edited July 31, 1996
      subroutine Lasa1(fil,rad,pop,skal,slut)
      character rad*200
      integer   i,fil,pop(15,0:10,0:1),skal
      logical   slut
      if (.NOT.slut) then
         do 5 i=1,200
   5        rad(i:i) = ' '
         read(fil,999,end=10) rad
         call Reada(rad,pop,skal,slut)
         return
      endif
   10 slut = .TRUE.
  999 format(A)
      return
      end
