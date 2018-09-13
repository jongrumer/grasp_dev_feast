************************************************************************
*                                                                      *
      SUBROUTINE DNICMV (N,M,B,C)
*                                                                      *
*   Matrix-matrix product: C = AB.  The lower triangle of the  (NxN)   *
*   matrix is assumed available in packed form in the array EMT. The   *
*   matrices B and C are (NxM).                                        *
*                                                                      *
*   This is an adaptation of  Andreas Stathopulos'  routine  TRSBMV,   *
*   and is specific to GRASP2 derivatives.                             *
*                                                                      *
*   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
*               [BLAS]: DAXPY/SAXPY, DDOT/SDOT.                        *
*                                                                      *
*   F A Parpia and A Stathopoulos         Last revision: 09 Oct 1992   *
*   MPI version by Xinghong He            Last revision: 18 Jun 1998   *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)
      POINTER (PIENDC,ENDCDUMMY)
      POINTER (PNIROW,IROWDUMMY)

      POINTER (PNTEMT,EMT(*))

      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT

      DIMENSION B(N,M),C(N,M)

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

*   Initialise the result matrix; note that this is specific to the
*   data structure of DVDSON --- there is no overdimensioning

      CALL DINIT (N*M,0.0D 00,C,1)

      ibeg = 1
      iend = 0
      DO ICOL = myid + 1, N, nprocs
         iend = iend + icol
         nelc = iend - ibeg + 1
         DO IV = 1, M
            DIAG =  C(ICOL,IV) + EMT(iend)*B(icol,IV)
            CALL DMERGE_dnicmv (NELC-1,B(1,IV),C(1,IV),
     :                     EMT(IBEG),B(ICOL,IV),DL)
            C(ICOL,IV) = DIAG + DL
         ENDDO
         ibeg = iend + 1
      ENDDO

      CALL gdsummpi (C,N*M)

      RETURN
      END

************************************************************************
*                                                                      *
      SUBROUTINE dmerge_dnicmv ( n, db, dc, da, dconst, dl )
C 
C  Used by dnimcv
C  Developed from dmerge. The only diff is: idy not needed here
C
************************************************************************
      IMPLICIT REAL*8          (A-H, O-Z)
      DIMENSION da(n), db(*), dc(*)

      dsum = 0.d0
      DO i = 1, n
         dsum = dsum + da(i) * db(i)
         dc(i) = dc(i) + dconst * da(i)
      ENDDO
      dl = dsum

      RETURN
      END
