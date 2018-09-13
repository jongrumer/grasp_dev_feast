************************************************************************
      SUBROUTINE SPICMVmpi (N,M,B,C)
      IMPLICIT REAL*8          (A-H,O-Z)

*   This routine now works for both rscfmpivu and rcimpivu. By removing
*   the include statement and the call to gdsummpi and setting myid=0,
*   nprocs=1, it also works for the corresponding serial program.
*   Matrix is stored in the mode of upper-triangle-by columns, or
*   you can say lower-triangle-by-rows. 98-08-06
*
*   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
*   lower triangle of the  (NxN)  matrix  A  is assumed available in   *
*   COMMON/HMAT/.                                                      *
*                                                                      *
*   This is an adaptation of  Andreas Stathopulos   routine  SPSBMV,   *
*   and is specific to GRASP2 derivatives.                             *
*                                                                      *
*   Call(s) to: [AUXBLAS]: DINIT/SINIT                                 *
*                                                                      *
*   F A Parpia and A Stathopoulos         Last revision: 19 Dec 1992   *
*   MPI version by Xinghong He            Last revision: 29 Jul 1998   *
*
************************************************************************

      integer*8 nelmnt


      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))

      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT

      DIMENSION B(N,M),C(N,M)

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
*
*   Initialise the result matrix; note that this is specific to the
*   data structure of DVDSON --- no overdimensioning
*
      CALL DINIT (N*M, 0.D0, C, 1)

      ibeg = 1
      DO ICOL = myid + 1, N, nprocs
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

      CALL gdsummpi (C, N*M)

      RETURN
      END
