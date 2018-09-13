************************************************************************
*                                                                      *
      SUBROUTINE BESSJ (W)
*                                                                      *
*   This routine evaluates Bessel fuctions J K ( W*R/C ) at the grid   *
*   points for K=L-1,L,L+1 and stores  them in the  arrays BJ(..,1),   *
*   BJ(..,2),BJ(..,3) respectively. It uses a power series expansion   *
*   for  small  r and switches to sin/cos expansion when more than 4   *
*   terms in power series are required.                                *
*                                         Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)

      include 'parameters.def'
CFF     PARAMETER (NNNP = 590) 
CFF     PARAMETER (NNN1 = NNNP+10)

      LOGICAL LDBPR
*
      COMMON/BESS2/BJ(NNNP,3),TC(NNNP),TD(NNNP)
     :      /DEBUGR/LDBPR(30)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /OSC2/L,KK
      COMMON/DEF2/C
     :      /DEF9/CVAC,PI

*
      EPSI = 1.0D-05
*
      IW = 1
      NN = L-1
      IF (KK .EQ. 0) GOTO 1
      IW = 2
      NN = L
    1 CONTINUE
      IEND = 2*NN+1
      IPROD = 1
      DO 2 I = 1,IEND,2
         IPROD = IPROD*I
    2 CONTINUE
      DFN = IPROD
      JCHAN = N
      BJ(1,IW) = 1.0D 00
      DO 5 J = 2,N
         WR = W*R(J)
         WA = -WR*WR*0.5D 00
         XBESS1 = 1.0D 00
         S1 = 1.0D 00
         DO 3 I = 1,4
            XBESS1 = XBESS1*WA/DBLE(I*(2*(NN+I)+1))
            S1 = S1+XBESS1
            IF (ABS (XBESS1) .LT. ABS (S1)*EPSI) GOTO 4
    3    CONTINUE
         JCHAN = J
         GOTO 6
    4    CONTINUE
         BJ(J,IW) = S1*(WR**NN)/DFN
    5 CONTINUE
*
*   Use sin/cos expansion when power series takes longer
*   than 4 terms to converge
*
    6 IF (JCHAN .GE. N) GOTO 22
      ISWAP = 0
      MODNN4 = MOD (NN-1,4)+1
      GOTO (7,8,9,10), MODNN4
*
*   NN = 1, 5, 9, ...
*
    7 SSN = -1.0D 00
      SCN =  1.0D 00
      ISWAP = 1
      GOTO 13
*
*   N = 2, 6, 10,....
*
    8 SSN = -1.0D 00
      SCN = -1.0D 00
      GOTO 13
*
*   NN = 3, 7, 11,...
*
    9 SSN =  1.0D 00
      SCN = -1.0D 00
      ISWAP = 1
      GOTO 13
*
*   NN = 0, 4, 8,...
*
   10 SSN = 1.0D 00
      SCN = 1.0D 00
*
   13 CONTINUE
      DO 21 J = JCHAN,N
         WA = W*R(J)
         IF (ISWAP .GT. 0) GOTO 14
         SN = SSN*SIN (WA)
         CN = SCN*COS (WA)
         GOTO 15
   14    CONTINUE
         SN = SSN*COS (WA)
         CN = SCN*SIN (WA)
   15    CONTINUE
         I = -1
         S1 = 0.0D 00
   16    I = I+1
         IF (I .GT. NN) GOTO 20
         IF (I) 18,17,18
   17    B = 1.0D 00/WA
         GOTO 19
   18    CONTINUE
         B = B*DBLE((NN+I)*(NN-I+1))/DBLE(2*I)/WA
   19    S1 = S1+B*SN
         SKEEP = SN
         SN = CN
         CN = -SKEEP
         GOTO 16
   20    CONTINUE
         BJ(J,IW) = S1
   21 CONTINUE
   22 CONTINUE
*
      IF ((NN .GE. L+1) .OR. (KK .EQ. 1)) THEN
*
*   Print out Bessel functions if (debug) option set
*
         IF (LDBPR(16)) THEN
            DO 23 JJ = 1,3
               JJJ = L-2+JJ
               WRITE (99,300) JJJ,(BJ(II,JJ),II = 1,N)
   23       CONTINUE
         ENDIF
*
*   All done
*
c zou
         DO  JJ = 1,3
               JJJ = L-2+JJ
           DO II = 1,N
            BJ(II,JJ)=BJ(II,JJ)*(C/CVAC)**JJJ
           END DO
         END DO    
c zou
         RETURN
      ELSE
         NN = NN+1
         IW = IW+1
         GOTO 1
      ENDIF
*
  300 FORMAT (/' Bessel function of order ',I3,/(1P,7D18.10))
*
      END
