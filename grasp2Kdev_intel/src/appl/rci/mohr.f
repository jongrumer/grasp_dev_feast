************************************************************************
*                                                                      *
      SUBROUTINE MOHR (N,KAPPA,Z,FZALFA)
*                                                                      *
*   The  function  F (Z*alpha)  for the  1s  2s  2p-  2p  symmetries   *
*   is computed here.    A value is obtained by interpolating in, or   *
*   extrapolating from, the table due to  P J Mohr.   See  P J Mohr,   *
*   At Data Nucl Data Tables 29 (1983) 453.                            *
*                                                                      *
*   Call(s) to: [LIB92]: INTERP.                                       *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
*----------------------------------------------------------------------*
*
*   Number of data points
*
      PARAMETER (NUMVAL = 12)
*
      DIMENSION VAL1S(NUMVAL),
     :          VAL2S(NUMVAL),VAL2P1(NUMVAL),VAL2P3(NUMVAL),
     :          ARG(NUMVAL)
*
*   1s data :
*
      DATA VAL1S/10.3168D 00, 4.6540D 00,
     :            3.2460D 00, 2.5519D 00,
     :            2.1351D 00, 1.8644D 00,
     :            1.6838D 00, 1.5675D 00,
     :            1.5032D 00, 1.4880D 00,
     :            1.5317D 00, 1.6614D 00/
*
*   2s data:
*
      DATA VAL2S/10.5468D 00, 4.8930D 00,
     :            3.5063D 00, 2.8391D 00,
     :            2.4550D 00, 2.2244D 00,
     :            2.0948D 00, 2.0435D 00,
     :            2.0650D 00, 2.1690D 00,
     :            2.3870D 00, 2.7980D 00/
*
*   2p- data:
*
      DATA VAL2P1/-0.1264D 00,-0.1145D 00,
     :            -0.0922D 00,-0.0641D 00,
     :            -0.0308D 00, 0.0082D 00,
     :             0.0549D 00, 0.1129D 00,
     :             0.1884D 00, 0.2934D 00,
     :             0.4530D 00, 0.7250D 00/
*
*   2p data:
*
      DATA VAL2P3/0.1235D 00, 0.1303D 00,
     :            0.1436D 00, 0.1604D 00,
     :            0.1794D 00, 0.1999D 00,
     :            0.2215D 00, 0.2440D 00,
     :            0.2671D 00, 0.2906D 00,
     :            0.3141D 00, 0.3367D 00/
*
*   Z data:
*
      DATA ARG/  1.0D 00, 10.0D 00, 20.0D 00,
     :          30.0D 00, 40.0D 00, 50.0D 00,
     :          60.0D 00, 70.0D 00, 80.0D 00,
     :          90.0D 00,100.0D 00,110.0D 00/
*
*----------------------------------------------------------------------*
*
*   Convergence criterion for interpolation
*
      DATA ACCY/1.0D-03/
*
*   Interpolate or issue error message as appropriate
*
      IF     (N .EQ. 1) THEN
         IF (KAPPA .EQ. -1) THEN
            CALL INTERP (ARG,VAL1S,NUMVAL,Z,VALUE,ACCY)
         ELSE
            WRITE (*,300)
            WRITE (*,301) N,KAPPA
            STOP
         ENDIF
      ELSEIF (N .EQ. 2) THEN
         IF     (KAPPA .EQ. -1) THEN
            CALL INTERP (ARG,VAL2S,NUMVAL,Z,VALUE,ACCY)
         ELSEIF (KAPPA .EQ.  1) THEN
            CALL INTERP (ARG,VAL2P1,NUMVAL,Z,VALUE,ACCY)
         ELSEIF (KAPPA .EQ. -2) THEN
            CALL INTERP (ARG,VAL2P3,NUMVAL,Z,VALUE,ACCY)
         ELSE
            WRITE (*,300)
            WRITE (*,301) N,KAPPA
            STOP
         ENDIF
      ELSE
         WRITE (*,300)
         WRITE (*,302) N
         STOP
      ENDIF
*
      FZALFA = VALUE
*
      RETURN
*
  300 FORMAT ('MOHR:')
  301 FORMAT (' Principal quantum number, ',I12,', kappa, ',1I3,'.')
  302 FORMAT (' Principal quantum number, ',1I2,
     :        ', Should be either 1 or 2.')
*
      END
