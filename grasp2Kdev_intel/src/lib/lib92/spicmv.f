************************************************************************
*                                                                      *
      SUBROUTINE SPICMV (N,M,B,C)
*                                                                      *
*   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
*   lower triangle of the  (NxN)  matrix  A  is assumed available in   *
*   COMMON/HMAT/.                                                      *
*                                                                      *
*   This is an adaptation of  Andreas Stathopulos   routine  SPSBMV,   *
*   and is specific to GRASP2 derivatives.                             *
*                                                                      *
*   There are currently 4 versions of spicmv.
*    lib92/spicmv  - for column-mode rscfblk and rciblk
*    rcimpi/spicmv - for row-mode rcimpi
*    rciblk/spicmv - for row-mode rciblk
*    rcii/spicmv   - for row mode rcii (mpi)
*   The last two are very similar; the only difference is those stuffs
*        myid, nprocs, gdsummpi
*
*   Call(s) to: [AUXBLAS]: DINIT/SINIT                                 *
*               [SPBLAS]: DAXPYI/SAXPYI, DDOTI/SDOTI                   *
*                                                                      *
*   F A Parpia and A Stathopoulos         Last revision: 19 Dec 1992   *
*
*   Block version for rscfblk and rciblk by
*   Xinghong He                           Last revision: 26 Nov 1997
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
*
      integer*8 nelmnt
      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
*
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
*
      DIMENSION B(N,M),C(N,M)
*
*   Initialise the result matrix; note that this is specific to the
*   data structure of DVDSON --- no overdimensioning
*
      CALL DINIT (N*M,0.0D 00,C,1)
*
      DO 2 ICOL = 1,N
            IBEG = IENDC(ICOL-1)+1
            IEND = IENDC(ICOL)
            NELC = IEND-IBEG+1
	    DO 1 IV = 1,M
               DIAG =  C(ICOL,IV)
     :                      +EMT(IBEG)*B(IROW(IBEG),IV)
               CALL DMERGE (NELC-1,B(1,IV),C(1,IV),
     :                     IROW(IBEG+1),EMT(IBEG+1),B(ICOL,IV),DL)
               C(ICOL,IV) = DIAG + DL
    1       CONTINUE
    2 CONTINUE
*
      RETURN
      END

*      subroutine dmerge ( n, db, dc, idy, da, dconst, dl )
C 
C  this merge version has the advantage of loading da(i)
C  and idy(i) only once.
C
*      IMPLICIT REAL*8          (A-H,O-Z)
*     DIMENSION DA(N), DB(*), DC(*), IDY(N)
*
*      dsum = 0.0
*      do 30 i = 1, n
*         dsum = dsum + da(i) * db(idy(i))
*         dc(idy(i)) = dc(idy(i)) + dconst * da(i)
* 30   continue
*      dl = dsum
*      return
*      end
