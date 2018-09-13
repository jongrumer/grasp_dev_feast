************************************************************************
*                                                                      *
      SUBROUTINE SEARCH (FIND)
*                                                                      *
*   This  routine locates circuits or loops of order  NC. NPOINT(NC)   *
*   are the  indices  of the points (triads) pertaining to the first   *
*   such loop found.  NPART  is the number of separate parts (groups   *
*   of contiguous points) on  the  axis of the flat graph. IPARTS is   *
*   the number of points in the smallest  part. IPARTL is the number   *
*   of points in  the  largest  part.  The  SUBROUTINE finds all the   *
*   possible loops  of  order  3  and 4. For NC .GE. 5, it looks for   *
*   only those who are partitionned in NPART .LE. 2. which can even-   *
*   tually reduce to a loop of  order  4  without breaking the basic   *
*   structure of the flat graph. ICROSS = -1, if lines cross.          *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1
      LOGICAL FIND
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
*
      DATA NAME/'SEARCH'/
*
*   Initialization
*
      FIND = .FALSE.
      NCM1 = NC-1
      NCM = NC-2
      ICROSS = 0
*
*   First treat two cases that do not involve do loops:
*
*   1. One isolated point, either the first or the last
*
      NPART = 1
      IPARTL = NC-1
      IPARTS = 1
*
*   A. First
*
      I1 = IFIRST
      K3 = 3
      K2 = 2
  200 JA = JDIAG(I1,1)
      JC = JDIAG(I1,K3)
*
      IF (JA .EQ. JC) THEN
         IF (NC .GT. 1) THEN
            WRITE (*,300) I1,K3,JA,JC,NC
            STOP
         ENDIF
         NPOINT(1) = I1
         GOTO 900
      ENDIF
*
      I2 = TAB1(JA,K2)
      I3 = TAB1(JC,K2)
*
      IF (ABS(IL(I3)-IL(I2))-NCM .LT. 0) THEN
         WRITE (*,301) I2,I3,JA,JC,K2,NC
         STOP
      ENDIF
*
      IF (ABS(IL(I3)-IL(I2))-NCM .GT. 0) THEN
*
*   B. Last
*
         IF (I1 .NE. IFIRST) GOTO 250
         I1 = ILAST
         K3 = 2
         K2 = 1
         GOTO 200
      ENDIF
*
      IC = 1
      NPOINT(IC) = I1
      I20 = MIN (I2,I3)
      I21 = IL(I20)
      I31 = I21+NCM1
*
      DO 203 II = I21,I31
         IC = IC+1
         NPOINT(IC) = IH(II)
  203 CONTINUE
*
      IF (NC .LE. 2) THEN
         IF (JDIAG(IFIRST,1) .NE. JDIAG(ILAST,1))
     :       CALL PHASE (I1,JDIAG,M4TRD)
         GOTO 900
      ENDIF
*
      IF (I1 .NE. ILAST) THEN
         IT = I2
         JT = JDIAG(ILAST,2)
         K4 = 2
         I4 = ILAST
      ELSE
         IT = I3
         JT = JDIAG(IFIRST,3)
         K4 = 3
         I4 = IFIRST
      ENDIF
*
      IF (IT .EQ. I20) CALL PHASE (I1,JDIAG,M4TRD)
      IF ((JT .EQ.JA) .OR. (JT.EQ.JC)) CALL CHANGE (I4,K4)
      GOTO 900
*
*   2. Two isolated points,first and last
*
  250 IF (NC .EQ. 1) RETURN
      IF (NC .LE. 3) GOTO 100
      IPARTL = NC-2
      IPARTS = 1
      I1 = IFIRST
      I2 = ILAST
      JA = JDIAG(I1,1)
      JB = JDIAG(I1,3)
*
      IF (TAB1(JA,2) .NE. I2) THEN
         JA = JDIAG(I1,3)
         JB = JDIAG(I1,1)
         IF (TAB1(JA,2) .NE. I2) GOTO 100
      ENDIF
*
      IF (JA .EQ. JDIAG(I2,1)) THEN
         JC = JDIAG(I2,2)
      ELSE
         JC = JDIAG(ILAST,1)
      ENDIF
*
      I3 = TAB1(JB,2)
      I4 = TAB1(JC,1)
      IDIST = IL(I4)-IL(I3)
*
      IF (ABS(IDIST)-(NCM-1) .LT. 0) THEN
         WRITE (*,302) I3,I4,JB,JC,IDIST,NC
         STOP
      ENDIF
      IF (ABS(IDIST)-(NCM-1) .EQ. 0) THEN
         NPOINT(1) = ILAST
         NPOINT(2) = IFIRST
         ICROSS = SIGN (1,IDIST)
         IC = 2
         I20 = MIN (I3,I4)
         I21 = IL(I20)
         I31 = I21+NCM
*
         DO 261 II = I21,I31
            IC = IC+1
            NPOINT(IC) = IH(II)
  261    CONTINUE
*
         IF (JA .EQ. JDIAG(IFIRST,1)) CALL CHANGE (IFIRST,3)
         IF (JA .EQ. JDIAG(ILAST,1)) CALL CHANGE (ILAST,2)
         GOTO 900
      ENDIF
*
*   First general case: all points in one group
*
  100 NPART = 1
      IPARTS = 0
      IPARTL = NC
      K3 = 1
*
      DO 101 IN = 1,NBNODE
         I = IH(IN)
  108    JA = JDIAG(I,K3)
         IF (I .NE. TAB1(JA,2))THEN
            I2 = TAB1(JA,2)
*
            IF (IL(I2)-IN-NCM1 .LT. 0) THEN
               WRITE (*,303) IN,I,I2,IL(I2),JA,NC
               STOP
            ENDIF
            IF (IL(I2)-IN-NCM1 .EQ. 0) THEN
               I21 = IL(I2)
               IC = 0
*
               DO 103 II = IN,I21
                  IC = IC+1
                  NPOINT(IC) = IH(II)
  103          CONTINUE
*
               IF (JA .EQ. JDIAG(IFIRST,3)) CALL CHANGE (IFIRST,3)
               IF (JA .EQ. JDIAG(ILAST,2)) CALL CHANGE (ILAST,2)
               GOTO 900
            ENDIF
         ENDIF
*
         IF (IN .EQ. 1) THEN
            IF (K3 .NE. 3) THEN
               K3 = 3
               GOTO 108
            ELSE
               K3 = 1
            ENDIF
         ENDIF
*
  101 CONTINUE
*
*   Search did not find loop NC .LE. 3
*
      IF (NC .LE. 3) RETURN
*
*   General case of loop partitionned in 2 groups. DO loop
*   on IPARTS
*
      NPART = 2
      NC2 = NC/2
      K3 = 1
      K2 = 1
*
      DO 400 IPS = 2,NC2
         JPS = IPS-1
         NBN = NBNODE-JPS
*
         DO 401 I1 = 1,NBN
            I = IH(I1)
            I2 = IH(I1+JPS)
  402       JA = JDIAG(I,K3)
            JD = JDIAG(I2,K2)
*
            IF (I .EQ. TAB1(JA,1)) THEN
               II2 = TAB1(JD,2)
               II1 = TAB1(JA,2)
            ELSE
               II1 = TAB1(JA,1)
               II2 = TAB1(JD,1)
            ENDIF
*
            IDIST = IL(II1)-IL(II2)
*
            IF (ABS (IDIST)-(NCM-JPS) .LT. 0) THEN
               WRITE (*,304) JPS,I1,I,I2,JA,JD,II1,II2,IDIST,NC
               STOP
            ENDIF
            IF (ABS (IDIST)-(NCM-JPS) .GT. 0) GOTO 420
            ICROSS = SIGN (1,IDIST)
            IC = 0
            I21 = IL(I2)
*
            DO 410 II = I1,I21
               IC = IC+1
               NPOINT(IC) = IH(II)
  410       CONTINUE
*
            I20 = MIN (II1,II2)
            I30 = MAX (II1,II2)
            I21 = IL(I20)
            I31 = IL(I30)
*
            DO 411 II = I21,I31
               IC = IC+1
               NPOINT(IC) = IH(II)
  411       CONTINUE
*
            IPARTS = IPS
            IPARTL = NC-IPS
            IF ((JDIAG(IFIRST,3) .EQ. JA) .OR.
     :          (JDIAG(IFIRST,3) .EQ. JD)) CALL CHANGE (IFIRST,3)
            IF ((JDIAG(ILAST,2) .EQ. JA) .OR.
     :          (JDIAG(ILAST,2) .EQ. JD)) CALL CHANGE (ILAST,2)
            GOTO 900
*
  420       IF (I1 .EQ. 1) THEN
               IF (K3 .EQ. 3) THEN
                  K3 = 1
                  GOTO 401
               ELSE
                  K3 = 3
                  GOTO 402
               ENDIF
            ENDIF
*
            IF (I2 .EQ. ILAST) THEN
               IF (K2 .NE. 2) THEN
                  K2 = 2
                  GOTO 402
               ENDIF
            ENDIF
*
  401    CONTINUE
  400 CONTINUE
*
*   SEARCH did not find circuit of order NC
*
      RETURN
*
*   Loop found
*
  900 FIND = .TRUE.
      CALL PRINTJ (NAME,10)
*
      RETURN
*
*   Error printout
*
  300 FORMAT (' Error in SEARCH. I1,K3,JA,JC,NC = ',5I5)
  301 FORMAT (' Error in SEARCH. I2,I3,JA,JC,K2,NC = ',6I5)
  302 FORMAT (' Error in SEARCH. I3,I4,JB,JC,IDIST,NC = ',6I5)
  303 FORMAT (' Error in SEARCH. IN,I,I2,IL(I2),JA,NC = ',6I5)
  304 FORMAT (' Error in SEARCH. JPS,I1,I,I2,JA,JD,II1,II2,IDIST,NC = '
     :       ,10I5)
*
      END
