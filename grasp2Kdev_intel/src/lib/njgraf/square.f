************************************************************************
*                                                                      *
      SUBROUTINE SQUARE
*                                                                      *
*   Reduces  a  circuit  of  order 4 in the two cases which are left   *
*   over by POLYGN, namely two disconnected groups of two points and   *
*   one group of two points plus  the  two  ends of the axis. In the   *
*   latter, the end of the axis is transferred  to the beginning. In   *
*   this  process,  one summation variable and two  6j  symbols  are   *
*   introduced.                                                        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      INTEGER ARR,TAB1
      LOGICAL SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /NAM/NAMSUB
*
      DATA NAME/'SQUARE'/
*
      NAMSUB = NAME
      MP = MP+1
      SUMVAR(MP) = .TRUE.
      IT1 = NPOINT(1)
      IT2 = NPOINT(2)
*
      IF (ICROSS .EQ. 1) THEN
         IT3 = NPOINT(3)
         IT4 = NPOINT(4)
         K23 = 3
         K32 = 2
      ELSE
         IT3 = NPOINT(4)
         IT4 = NPOINT(3)
         K23 = 2
         K32 = 3
      ENDIF
*
      L4 = JDIAG(IT2,1)
*
      IF (ARR(IT2,1) .LE. 0) THEN
         CALL PHASE2 (L4)
         ARR(IT2,1) = 1
         ARR(IT3,1) = -1
      ENDIF
*
      L2 = JDIAG(IT1,1)
      IF (ARR(IT1,1) .GT. 0) CALL PHASE2 (L2)
      JWC = JWC+1
      KW(1,JWC) = L4
      KW(2,JWC) = L2
      KW(3,JWC) = JDIAG(IT2,2)
      JJ1 = JDIAG(IT1,3)
      KW(4,JWC) = JJ1
      KW(5,JWC) = MP
      KW(6,JWC) = JDIAG(IT1,2)
      IF (ARR(IT1,2) .LT. 0) CALL PHASE2 (JDIAG(IT1,2))
      JWC = JWC+1
      KW(1,JWC) = L4
      KW(2,JWC) = L2
      JJ3 = JDIAG(IT3,K23)
      JJ2 = JDIAG(IT4,K32)
      KW(3,JWC) = JJ3
      KW(4,JWC) = JJ2
      KW(5,JWC) = MP
      KW(6,JWC) = JDIAG(IT3,K32)
      IF (ARR(IT3,K32) .LT. 0) CALL PHASE2 (JDIAG(IT3,K32))
      J6(J6C+1) = MP
      J6C = J6C+2
      J6(J6C) = MP
*
      IF (NPART .EQ. 1) THEN
         ITMIN = IT2
         ITMAX = IT3
      ELSE
         ITMIN = MIN (IT2,IT3)
         ITMAX = MAX (IT2,IT3)
      ENDIF
      ITMN = MIN (IT1,IT4)
      ITMX = MAX (IT1,IT4)
*
      TAB1(MP,1) = ITMIN
      TAB1(MP,2) = ITMAX
      JDIAG(IT2,1) = MP
      JDIAG(IT3,1) = MP
      JDIAG(IT2,3) = JJ1
      ARR(IT2,3) = ARR(IT1,3)
      JDIAG(IT3,K32) = JJ2
      ARR(IT3,K32) = ARR(IT4,K32)
*
      IF (ICROSS .EQ. 1) THEN
         J7(J7C+1) = L2
         J7(J7C+2) = L4
         CALL PHASE2 (L4)
         J7C = J7C+3
         J7(J7C) = MP
      ELSE
         CALL PHASE2 (JJ2)
      ENDIF
*
      ITLL = IL(ITMN)
      ITHL = IL(ITMX)
*
      DO 5 I = ITLL+1,ITHL-1
         IT = IH(I)
         ILP = I-1
         IL(IT) = ILP
         IH(ILP) = IT
    5 CONTINUE
      IF (ITHL .NE. NBNODE) THEN
         DO 6 I = ITHL+1,NBNODE
            IT = IH(I)
            ILP = I-2
            IL(IT) = ILP
            IH(ILP) = IT
    6    CONTINUE
      ENDIF
*
      IF (NPART .NE. 2) THEN
         TAB1(JJ1,1) = IH(1)
         TAB1(JJ1,2) = IH(NBNODE-2)
      ENDIF
*
      RETURN
      END
