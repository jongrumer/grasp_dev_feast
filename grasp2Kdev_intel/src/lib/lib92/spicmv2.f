************************************************************************
      SUBROUTINE spicmv2 (N,M,B,C)
      IMPLICIT REAL*8          (A-H,O-Z)

*  Modified from the mpi version spicmvmpi.f by simply removing things
*  required by mpi communication.

*  Xinghong He  98-07-29

************************************************************************

      integer*8 nelmnt

      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))

      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT

      DIMENSION B(N,M),C(N,M)
!-----------------------------------------------------------------------
*
*   Initialise the result matrix; note that this is specific to the
*   data structure of DVDSON --- no overdimensioning
*
      CALL DINIT (N*M, 0.D0, C, 1)

      ibeg = 1
      DO ICOL = 1, N
            !IBEG = IENDC(ICOL-1)+1
            !IEND = IENDC(ICOL)
         IEND = IENDC(ICOL)
         NELC = IEND - IBEG + 1
         DO IV = 1, M
            DIAG =  C(ICOL,IV) + EMT(iend)*B(icol,IV)
            CALL DMERGE (NELC-1,B(1,IV),C(1,IV),
     :                   IROW(IBEG),EMT(IBEG),B(ICOL,IV),DL)
            C(ICOL,IV) = DIAG + DL
         ENDDO
         ibeg = iend + 1
      ENDDO

      RETURN
      END
