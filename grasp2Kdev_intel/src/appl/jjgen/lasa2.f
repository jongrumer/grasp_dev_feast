*     last edited July 30, 1996
      subroutine Lasa2(fil,rad2,rad3,stopp,slut)
      character rad2*200,rad3*200
      integer   fil,stopp
      logical   slut
      if (.NOT.slut) then
         read(fil,999,end=10) rad2
C        read(fil,999,end=10) rad2(1:stopp)
         read(fil,999,end=10) rad3
C        read(fil,999,end=10) rad3(1:stopp+4)
         return
      endif
   10 slut = .TRUE.
  999 format(A)
      return
      end
