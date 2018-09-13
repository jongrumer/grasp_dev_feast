************************************************************************
*                                                                      *
      SUBROUTINE DIAGRM (JUMP)
*                                                                      *
*   This subroutine builds up a flat diagram from the triads J23 and   *
*   places them in JDIAG . Arrows  are in ARR (INTEGER). The diagram   *
*   is built so as to maximize the number of triads involved, within   *
*   a one-step-forward-check process. If the diagram does not inclu-   *
*   de all the NBTR triads, it will have 'free ends'. JDIAG has dim-   *
*   ension double that of  J23 , because the path may proceed either   *
*   way.                                                               *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1,ARROW
      LOGICAL TABS,FREE,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DATA NAME/'DIAGRM'/
*
      DATA NB/0/
*
*   Initialization
*
      IF (JUMP .GT. 2) GOTO 17
      IF (JUMP .LT. 2) NB = 0
    1 NB = NB+1
      IF (TABS(NB)) GOTO 1
      NODE = NBTR
      ILAST = NBTR
*
      DO 2 J = 1,3
         JDIAG(NODE,J) = J23(NB,J)
         ARR(NODE,J) = ARROW(NB,J)
    2 CONTINUE
*
      TABS(NB) = .TRUE.
*
      DO 15 I = 1,MP
         IAL(I) = 0
   15 CONTINUE
*
      IF1 = JDIAG(NODE,1)
      IF2 = JDIAG(NODE,3)
      IAL(IF1) = 1
      IAL(IF2) = 1
   17 NTIME = 0
      I1 = 1
      K1 = 1
      K2 = 2
      K3 = 3
    3 JB = JDIAG(NODE,K2)
      CALL OTHERJ (0,JB,L,LC,KP)
      CALL NEIBOR (LC,L1,L2)
*
*   Check consistency of triads
*
      IF (TABS(L)) THEN
         WRITE (*,300)
         STOP
      ENDIF
*
      CALL WAY (L,L1,L2,ICH,ND)
      NODE = NODE+I1
      TABS(L) = .TRUE.
      JDIAG(NODE,K3) = J23(L,LC)
      ARR(NODE,K3) = ARROW(L,LC)
      ICT = ICH*I1
*
      IF (ICH .LE. 0) THEN
         LP = L1
         L1 = L2
         L2 = LP
      ENDIF
*
      IF (ICT .LE. 0) CALL PHASE (L,J23,M2TRD)
      JDIAG(NODE,K1) = J23(L,L1)
      ARR(NODE,K1) = ARROW(L,L1)
      JDIAG(NODE,K2) = J23(L,L2)
      ARR(NODE,K2) = ARROW(L,L2)
      J = J23(L,L1)
      IAL(J) = IAL(J)+1
      J = J23(L,L2)
      IAL(J) = IAL(J)+1
      IF (ND .LT. 1) GOTO 3
      NTIME = NTIME+1
      ILAST = MAX (NODE,ILAST)
      IFIRST = MIN (NODE,NBTR)
      NBP = IAL(IF1)+IAL(IF2)
      IF ((NBP .GT. 3) .OR. (NTIME .GT. 1)) THEN
         NBNODE = ILAST-IFIRST+1
         NBTR = NBTR-NBNODE
*
*   Definition of free ends and other quantities.
*
         CALL INTAB
         CALL PRINTJ (NAME,12)
         GOTO 50
      ENDIF
*
      IF (NBP .GT. 2) THEN
         IF (IAL(IF1) .LE. IAL(IF2)) THEN
            JT = JDIAG(NBTR,1)
            JAR = ARR(NBTR,1)
            JDIAG(NBTR,1) = JDIAG(NBTR,3)
            ARR(NBTR,1) = ARR(NBTR,3)
            JDIAG(NBTR,3) = JT
            ARR(NBTR,3) = JAR
            CALL PHASE (NBTR,JDIAG,M4TRD)
         ENDIF
      ENDIF
*
      NODE = NBTR
      I1 = -1
      K2 = 3
      K3 = 2
      GOTO 3
*
   50 RETURN
*
  300 FORMAT ('DIAGRM: Flat graph impossible to build.')
*
      END
