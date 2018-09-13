************************************************************************
*                                                                      *
        SUBROUTINE genintbreit2wrap (myid, nprocs, j2max)

*   Written by Per Jonsson                Last revision: October 2014  *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      POINTER (PINDT1,INDTP1(*))
      POINTER (PVALT1,VALTP1(*))
      POINTER (PINDT2,INDTP2(*))
      POINTER (PVALT2,VALTP2(*))

      POINTER (PINDT3,INDTP3DUMMY(*))
      POINTER (PVALT3,VALTP3DUMMY(*))
      POINTER (PINDT4,INDTP4DUMMY(*))
      POINTER (PVALT4,VALTP4DUMMY(*))
      POINTER (PINDT5,INDTP5DUMMY(*))
      POINTER (PVALT5,VALTP5DUMMY(*))
      POINTER (PINDT6,INDTP6DUMMY(*))
      POINTER (PVALT6,VALTP6DUMMY(*))

      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)

!-----------------------------------------------------------------------

      CALL genintbreit2 ((myid), (nprocs), N, j2max)

! Gather integrals (and their indeces) from- and send to- all nodes

      CALL gisummpi (INDTP2, N)
      CALL gdsummpi (VALTP2, N)

      RETURN
      END
