************************************************************************
*                                                                      *
      SUBROUTINE PRINTJ (NAMES,JP)
*                                                                      *
*   This  SUBROUTINE  prints  intermediate  results in standard form   *
*   from wherever it is called.                                        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*8 IBLANK,IFREE,IFR
      CHARACTER*6 NAMES,NSETTB
      CHARACTER*4 I6,I7,I8,I9,IJ1
      CHARACTER IM,IP,IS(3)
      INTEGER ARR,TAB1,ARROW
      LOGICAL TABS,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10,MTAB = 30,MZERO = 20)
*
      DIMENSION IX(7),JTAB(MTAB,3)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CONST/I6C,I7C,I8C,I9C,IDEL,IWC
     :      /ZER/NZERO,JZERO(MZERO)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
*
      EQUIVALENCE(I6C,IX(1))
*
      DATA IBLANK,IFREE,IP,IM/'        ','FREE END','+','-'/
      DATA NSETTB/'SETTAB'/
      DATA I6,I7,I8,I9,IJ1/'I6=','I7=','I8=','I9=','J1='/
*
      IF (IBUG3 .NE. 1) RETURN
      WRITE (99,1050) NAMES
*
*   Initialise variables
*
      I6C = 1
      I7C = 1
      I8C = 1
      I9C = 1
      IDEL = 1
      IWC = 1
*
      JUMP = JP
      IF (JUMP .EQ. 0) THEN
*
         DO 9 I = 1,7
            IX(I) = 1
    9    CONTINUE
*
         WRITE (99,1020) IJ1,(J1(I),I = 1,M)
      ENDIF
*
      IF (JUMP .LT. 8) GOTO 20
      WRITE (99,1000) NBNODE,NBTR,NFIN,IFIRST,ILAST,NFREE
      JUMP = JUMP-8
      WRITE (99,1001)
      K = 0
*
      DO 1 I = 1,NBNODE
         IT = IH(I)
         IFR = IBLANK
         JT = JDIAG(IT,1)
*
         IF ((TAB1(JT,2) .NE. IT) .OR. (JT .EQ. JDIAG(IFIRST,3))) THEN
            K = K+1
            IF (K .GT. MTAB) THEN
               WRITE (*,100) K,MTAB
               STOP
            ENDIF
            JTAB(K,1) = JT
            JTAB(K,2) = TAB1(JT,1)
            JTAB(K,3) = TAB1(JT,2)
         ENDIF
*
         IF (TAB1(JT,2) .GT. ILAST) IFR = IFREE
*
         DO 2 J = 1,3
            IS(J) = IP
            IF (ARR(IT,J) .LT. 1) IS(J) = IM
    2    CONTINUE
*
         WRITE (99,1002) (IS(J),J = 1,3)
         WRITE (99,1003) IL(IT),IT,IFR,(JDIAG(IT,J),J = 1,3)
*
    1 CONTINUE
*
      WRITE (99,1004)
      NTIME = 0
      JT = JDIAG(IFIRST,3)
      IF (JT .NE. JDIAG(ILAST,2)) THEN
         IF (TAB1(JT,2) .LT. 1000) GOTO 5
      ENDIF
    4 K = K+1
      IF (K .GT. MTAB) THEN
         WRITE (*,101) K,MTAB
         STOP
      ENDIF
      JTAB(K,1) = JT
      JTAB(K,2) = TAB1(JT,1)
      JTAB(K,3) = TAB1(JT,2)
    5 NTIME = NTIME+1
*
      IF (NTIME .NE. 2) THEN
         JT = JDIAG(ILAST,2)
         IF (TAB1(JT,2) .EQ. 1000) GOTO 4
      ENDIF
*
      WRITE (99,1005) ((JTAB(I,J),J = 1,3),I = 1,K)
      WRITE (99,1006) (I,SUMVAR(I),I = 1,MP)
   20 IF (JUMP .LT. 4) GOTO 30
      JUMP = JUMP-4
      NBTR1 = 2*N-2
      WRITE (99,1010) NBTR1
      K = 0
*
      DO 11 I = 1,NBTR1
         IF (TABS(I)) GOTO 11
         K = K+1
*
         DO 12 J = 1,3
            IS(J) = IP
            IF (ARROW(I,J) .LT. 1) IS(J) = IM
   12    CONTINUE
*
         WRITE (99,1012) (IS(J),J = 1,3)
         WRITE (99,1013) K,I,(J23(I,J),J = 1,3)
*
   11 CONTINUE
*
      WRITE (99,1014)
      MM = M
      IF (NAMES .NE. NSETTB) MM = M-1
      WRITE (99,1015) (I,(LINE(I,J),LCOL(I,J),J = 1,2),I = 1,MM)
*
   30 IF (JUMP .GE. 2) THEN
         JUMP = JUMP-2
         WRITE (99,1030) NC,NPART,IPARTL,IPARTS,ICROSS,
     :                    (NPOINT(I),I = 1,NC)
      ENDIF
*
      IF (JUMP .GE. 1) WRITE (99,1040) NZERO,(I,JZERO(I),I = 1,NZERO)
      IF (J6C .GE. I6C) WRITE (99,1020) I6,(J6(I),I = I6C,J6C)
      IF (J7C .GE. I7C) WRITE (99,1020) I7,(J7(I),I = I7C,J7C)
      IF (J8C .GE. I8C) WRITE (99,1020) I8,(J8(I),I = I8C,J8C)
      IF (J9C .GE. I9C) WRITE (99,1020) I9,(J9(I),I = I9C,J9C)
      IF (JDEL .GE. IDEL) WRITE (99,1021)
     :                    ((LDEL(I,J),J = 1,2),I = IDEL,JDEL)
      IF (JWC .GE. IWC) WRITE (99,1022)
     :                    ((KW(J,I),J = 1,6),I = IWC,JWC)
      I6C = J6C+1
      I7C = J7C+1
      I8C = J8C+1
      I9C = J9C+1
      IDEL = JDEL+1
      IWC = JWC+1
      RETURN
*
  100 FORMAT (' Dimension error in PRINTJ. K = ',I5,' MTAB = ',I5)
  101 FORMAT (' Dimension error IN PRINTJ. K = ',I5,' MTAB = ',I5)
 1000 FORMAT (/10X,'NBNODE = ',I3,10X,'NBTR = ',I3,10X,'NFIN = ',I3,
     :   /10X,'IFIRST = ',I3,10X,'ILAST = ',I3,9X,'NFREE = ',I3)
 1001 FORMAT (//7X,'IL',3X,'IH',14X,'JDIAG'//)
 1002 FORMAT (28X,3(A1,2X))
 1003 FORMAT (7X,I2,3X,I2,2X,A8,2X,3I3/)
 1004 FORMAT (/5X,'TAB1'/)
 1005 FORMAT (4(I3,1H),2X,I3,I5,5X))
 1006 FORMAT (/2X,'SUMVAR = ',15(I3,L1))
 1010 FORMAT (//10X,'J23',10X,'NBTR1 = ',I3//)
 1012 FORMAT (18X,3(A1,2X))
 1013 FORMAT (I9,I5,2X,3I3/)
 1014 FORMAT (/3X,'J  L1 K1  L2 K2')
 1015 FORMAT (4(I4,1H),I3,I3,I4,I3))
 1020 FORMAT (/3X,A4,3X,3(20I3/))
 1021 FORMAT (/3X,'DELTA = ',7(I5,I3))
 1022 FORMAT (/3X,'KW(ARG. OF 6J)',6I3)
 1030 FORMAT (//2X,'NC = ',I2,4X,'NPART = ',I2,4X,'IPARTL = ',I2,4X,
     :   'IPARTS = ',I2,4X,'ICROSS = ',I2,4X,/2X,'NPOINT = ',20I3)
 1040 FORMAT (//2X,'NZERO = ',I2,5X,12(I4,1H),I3))
 1050 FORMAT (///3X,'Printout after calling SUBROUTINE ',A7)
*
      END
