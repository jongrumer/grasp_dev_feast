************************************************************************
*                                                                      *
      SUBROUTINE CFP7 (NEL,IJD,IVD,IJP,IVP,COEFP)
*                                                                      *
*   Table look-up for fractional parentage coefficients of  equival-   *
*   ent electrons with j = 7/2. See listing of CFP for argument list.  *
*                                                                      *
*                                           Last update: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
*
      DIMENSION IJ(4,8),IV(4,8),NUM3(6,4),NORM3(6),NUM4(8,6),NORM4(8)
*
*   1.0 Tables of data
*
*   1.1 Floating point constants
*
      PARAMETER (D0 = 0.0D 00,
     :           D1 = 1.0D 00,
     :           D9 = 9.0D 00)
*
*   1.2 Table for CFPs.
*
      DATA IJ/7,0,7,0,0,4,3,4,0,8,5,8,0,12,9,12,0,0,11,4,0,0,15,8,3*0,10
     :         ,3*0,16/
      DATA IV/1,0,1,0,10,2,3,2,10,2,3,2,10,2,3,2,10,10,3,4,10,10,3,4,3*1
     :0,4,3*10,4/
      DATA NUM3/9,5*0,-5,3,121,143,-55,0,-9,-11,12,-900,39,5,-13,0,-65,3
     :43,104,-17/
      DATA NORM3/36,14,198,1386,198,22/
      DATA NUM4/1,280,308,1144,4*0,0,54,-121,0,-968,169,462,0,0,-231,-14
     :,195,-77,2366,-343,0,0,-65,250,-245,-1755,90,-945,140,0,-210,91,62
     :4,280,2275,650,234,0,0,140,-1224,0,-560,680,627/
      DATA NORM4/1,840,924,3432,3080,5460,3080,1001/
*
*   2.0 Locate entry in CFP table.
*
      IF (NEL .LE. 0) GOTO 14
      IF (NEL .GE. 5) GOTO 1
*
      N = NEL
      IJ1 = IJD
      IV1 = IVD
      IJ2 = IJP
      IV2 = IVP
      GOTO 2
    1 IF (NEL .GT. 8) GOTO 14
      N = 9-NEL
      IJ1 = IJP
      IV1 = IVP
      IJ2 = IJD
      IV2 = IVD
*
*   2.1 Find 'daughter' index.
*
    2 K = 0
    3 K = K+1
      IF (K .GT. 8) GOTO 14
      IF (IJ(N,K) .NE. IJ1) GOTO 3
      IF (IV(N,K) .NE. IV1) GOTO 3
      KD = K
*
*   2.2 Find 'parent' index.
*
      IF (N .NE. 1) GOTO 4
      IF (IV2 .NE. 0) GOTO 14
      IF (IJ2 .EQ. 0) GOTO 6
      GOTO 14
    4 K = 0
    5 K = K+1
      IF (K .GT. 8) GOTO 14
      IF (IJ(N-1,K) .NE. IJ2) GOTO 5
      IF (IV(N-1,K) .NE. IV2) GOTO 5
      KP = K
*
*   3.0 Compute coefficients.
*
*   3.1 Table look-up
*
      GOTO (6,6,7,11),N
    6 COEFP = D1
      GOTO 12
    7 CONTINUE
      COEFP = DBLE (NUM3(KD,KP))
      DENOM = DBLE (NORM3(KD))
    8 IF (COEFP) 10,13,9
    9 CONTINUE
      COEFP = SQRT (COEFP/DENOM)
      GOTO 12
   10 CONTINUE
      COEFP = -SQRT (-COEFP/DENOM)
      GOTO 12
   11 CONTINUE
      COEFP = DBLE (NUM4(KD,KP))
      DENOM = DBLE (NORM4(KD))
      GOTO 8
*
*   3.2 Insert additional factors for hole states
*
   12 IF (NEL .LE. 4) GOTO 13
      DNEL = DBLE (NEL)
      FACT = ((D9-DNEL)/DNEL)
     :       *(D1+DBLE(IJP))/(D1+DBLE(IJD))
      COEFP = COEFP* SQRT (FACT)
C GG tikras      IS = IABS ((IJD-IJP-IVD+IVP)/2-3)
      IS = IABS ((IJD-IJP+IVD-IVP)/2)
      IF (MOD (IS,2) .EQ. 0) GOTO 13
      COEFP = -COEFP
   13 RETURN
*
*   4.0 Fault mode section.
*
   14 WRITE (*,300) NEL,IJD,IVD,IJP,IVP
      STOP
*
  300 FORMAT ('CFP7: Error in trying to compute a CFP',
     :        ' for a state with ',1I2,' electrons with j = 7/2;'
     :       /' IJD = ',1I2,', IVD = ',1I2,
     :       ', IJP = ',1I2,', IVP = ',1I2,'.')
*
      END
