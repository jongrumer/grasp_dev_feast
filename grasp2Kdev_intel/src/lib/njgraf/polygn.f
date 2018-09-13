************************************************************************
*                                                                      *
      SUBROUTINE POLYGN (JPOL)
*                                                                      *
*   This routine reduces a circuit of arbitrary order NC. It exchan-   *
*   ges nodes on the flat diagram until the distance on the axis be-   *
*   tween nodes equals one. Each exchange introduces a summation va-   *
*   riable  and  a 6j-symbol. The circuit has a maximum of NPART = 2   *
*   disconnected parts on the axis.                                    *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
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
*
      DATA NAME/'POLYGN'/
*
      NC1 = NC+1
      NC2 = NC
      NBC = IPARTL-2
*
   10 DO 8 I = 1,NBC
         IT2 = NPOINT(NC1-I)
         IT1 = NPOINT(NC2-I)
         JB = JDIAG(IT1,1)
         JC = JDIAG(IT2,1)
         JDIAG(IT1,1) = JC
         JDIAG(IT2,1) = JB
         JAR = ARR(IT1,1)
         ARR(IT1,1) = ARR(IT2,1)
         ARR(IT2,1) = JAR
         JE = JDIAG(IT1,2)
         MP = MP+1
         SUMVAR(MP) = .TRUE.
         JDIAG(IT1,2) = MP
         JDIAG(IT2,3) = MP
*
         IF (TAB1(JB,1) .EQ. IT1) THEN
            TAB1(JB,1) = IT2
         ELSE
            TAB1(JB,2) = IT2
         ENDIF
*
         IF (TAB1(JC,1) .EQ. IT2) THEN
            TAB1(JC,1) = IT1
         ELSE
            TAB1(JC,2) = IT1
         ENDIF
*
         IF (ARR(IT1,2) .LE. 0) THEN
            CALL PHASE2 (JE)
            ARR(IT1,2) = 1
            ARR(IT2,3) = -1
         ENDIF
*
         JWC = JWC+1
         KW(1,JWC) = JB
         KW(2,JWC) = MP
         KW(3,JWC) = JE
         KW(4,JWC) = JC
         KW(5,JWC) = JDIAG(IT2,2)
         KW(6,JWC) = JDIAG(IT1,3)
         J6(J6C+1) = MP
         J6C = J6C+2
         J6(J6C) = MP
    8 CONTINUE
*
      NC = NC-NBC
*
      IF (NC .GT. 4) THEN
         NBC = IPARTS-2
         NC1 = IPARTS+1
         NC2 = IPARTS
         GOTO 10
      ENDIF
*
      IF (NPART .NE. 1) THEN
         NPOINT(3) = NPOINT(NC1)
         NPOINT(4) = NPOINT(NC1+1)
      ENDIF
*
      IF (NC .EQ. 2) JPOL = 1
      CALL PRINTJ (NAME,10)
*
      RETURN
      END
