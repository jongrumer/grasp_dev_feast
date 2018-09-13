************************************************************************
*                                                                      *
        SUBROUTINE genintrkwrap (myid, nprocs, j2max)

*   Written by     Xinghong He            Last revision: 12 Jun 1998   *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      POINTER (PCTEVLRK,VALTEIRK(*))
      POINTER (PCTEILRK, INDTEIRK(*))
      COMMON/CTEILSRK/PCTEILRK,PCTEVLRK

!-----------------------------------------------------------------------

      CALL genintrk ((myid), (nprocs), N, j2max)

! Gather integrals (and their indeces) from- and send to- all nodes

      CALL gisummpi (INDTEIRK, N)
      CALL gdsummpi (VALTEIRK, N)

      RETURN
      END
