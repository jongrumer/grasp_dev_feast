************************************************************************
*                                                                      *
      SUBROUTINE setNCF
*                                                                      *
*   Determines NCF from the  .csl  file. 
*   the  .csl  file must be OPENed before the call to setNCF
*   and  .csl   has to be REWIND-ed after the call to setNCF
*                                                                      *
************************************************************************
*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     implicit none
      IMPLICIT REAL*8          (A-H,O-Z)

      integer :: NCF 
      integer       :: i,n
      character*1000 :: string
      character*600  :: csf1,csf2,csf3

      COMMON /ORB2/NCF,NW,PNTRIQ

! Read header 

      do i = 1,5
        read(21,'(a)') string
      end do

! Read CSFs 

      n = 0
      do 
        read(21,'(a)',end = 99) csf1
       if (csf1(1:2) .ne. ' *') then
        read(21,'(a)') csf2
        read(21,'(a)') csf3
         n = n + 1
       else
         cycle
       endif
      end do

  99  continue

      NCF = n 

      print *, '  setNCF: NCF = ', NCF

      return
      end 
