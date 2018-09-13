************************************************************************
*                                                                      *
      FUNCTION FUNL (X,K)
*                                                                      *
*   This  function  evaluates the LK(X) functions using the analytic   *
*   functions defined  in table 5  and equations  (20) and  (21)  of   *
*   Fullerton and Rinker.                                              *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      DIMENSION F(6,2),G(3,2),H(2,2)
*
      DATA F/2.008188D 00, -2.397605D 00,
     :       1.046471D 00, -3.670660D-01,
     :       6.374000D-02, -3.705800D-02,
     :       1.646407D 00, -2.092942D 00,
     :       9.623100D-01, -2.549600D-01,
     :       1.644040D-01,  0.0D 00/
*
      DATA G/7.51198D-01,  1.38889D-01,
     :       2.0886D-02,   1.37691D-01,
     :      -4.16667D-01, -9.7486D-02/
*
      DATA H/-4.44444D-01, -3.472D-03,
     :        4.44444D-01,  1.7361D-02/
*
      DATA A,B/2.2D 00, -1.72D 00/
*
      IF ((K .LT. 0) .OR. (K .GT. 1)) GOTO 99
      IF (X .GT. 2.0D 00) GOTO 3
      IF (X .EQ. 0.0D 00) GOTO 6
*
*   Use rational approximation for X < 2
*
      K1 = K+1
      SUM = 0.0D 00
      XN = 1.0D 00
      DO 1 I = 1,6
         SUM = SUM+XN*F(I,K1)
         XN = XN*X
    1 CONTINUE
      X2 = X*X
      SUMG = G(1,K1)+X2*(G(2,K1)+X2*G(3,K1))
      SUMH = H(1,K1)+X2*X2*H(2,K1)
      XN = LOG(X)
      SUMG = XN*(SUMG+XN*SUMH)
      IF (K .EQ. 0) GOTO 2
      SUM = SUM+SUMG
      GOTO 7
    2 SUM = SUM+X*SUMG
      GOTO 7
*
    3 CONTINUE
      SUM = A+B/X
      IF (K .EQ. 0) GOTO 4
      SUM = SUM+(SUM+B/X)/X
    4 SUM = SUM/X
      XM = -X
      SUM = SUM*EXP (XM)
      GOTO 7
    6 IF (K .EQ. 1) GOTO 98
      SUM = F(1,1)
    7 FUNL = SUM
      RETURN
*
*   Error section
*
   98 WRITE (*,302)
      STOP
   99 WRITE (*,301)
      STOP
*
  301 FORMAT (/'FUNL: K must be either 0 or 1')
  302 FORMAT (/'FUNL: Attempt to calculate function for'
     :        /' zero argument and K value of 1')
*
      END
