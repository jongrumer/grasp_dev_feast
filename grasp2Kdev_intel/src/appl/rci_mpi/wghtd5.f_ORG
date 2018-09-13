************************************************************************
*                                                                      *
      SUBROUTINE WGHTD5(iatjpo, iaspar)
*                                                                      *
*   Print  the  weights of the largest five CSF contributors to each   *
*   ASF.                                                               *
*                                                                      *
*   Call(s) to: ALLOC, DALLOC, ISPAR, ITJPO.                           *
*                                                                      *
*                                          Last updated: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*4 JLBL,LABJ,LABP
*
      POINTER (PNTRWT,WT(1))
      POINTER (PNNEXT,NEXT(1))
*
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNIVEC,IVEC(1))
*
      DIMENSION WGHT(5),ICONF(5)
*
      COMMON/EIGVEC/PNEVEC
     :      /JLABL/JLBL(32),LABJ(32),LABP(2)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
*
*   Allocate storage for local arrays
*
      CALL ALLOC (PNTRWT,NCF,8)
      CALL ALLOC (PNNEXT,NCF,4)
*
      WRITE (24,300)
*
      IF (NCF .LT. 5) NELT = NCF
*
      NVEX = NVEC
*
      DO 7 IV = 1,NVEX
*
         ICF = IVEC(IV)
*
*   Set up linked list of weights
*
         NEXT(1) = 0
         WT(1) = EVEC(1+(IV-1)*NCF)**2
         IFIRST = 1
         DO 4 I = 2,NCF
            M = IFIRST
            L = 0
            WT(I) = EVEC(I+(IV-1)*NCF)**2
    1       IF (WT(I) .LE. WT(M)) GOTO 2
            IF (L .NE. 0) GOTO 3
            NEXT(I) = IFIRST
            IFIRST = I
            GOTO 4
    2       L = M
            M = NEXT(L)
            IF (M .NE. 0) GOTO 1
    3       NEXT(I) = NEXT(L)
            NEXT(L) = I
    4    CONTINUE
*
*   Print first five elements of list.
*
         M = IFIRST
         I = 0
    5    continue
!	 IF ((ITJPO (M) .NE. IATJPO(IV)) .OR.
!     :       (ISPAR (M) .NE. IASPAR(IV))) GOTO 6
         I = I+1
         WGHT(I) = WT(M)
         ICONF(I) = M
    6    M = NEXT(M)
         IF (M .NE. 0.AND.I .LT. 5) GOTO 5
         IP = (iaspar + 3 ) / 2
         NELT = MIN (I,5)
         WRITE (24,301) ICF,LABJ(IATJPO),LABP(IP),(WGHT(I),I = 1,NELT)
         WRITE (24,302) (ICONF(I),I = 1,NELT)
    7 CONTINUE
*
*   Deallocate storage for local arrays
*
      CALL DALLOC (PNTRWT)
      CALL DALLOC (PNNEXT)
*
      RETURN
*
  300 FORMAT (/'Weights of major contributors to ASF:'
     :       //'Level J Parity      CSF contributions'/)
  301 FORMAT (I3,2X,2A4,5(E12.4))
  302 FORMAT (13X,5(I12))
*
      END
