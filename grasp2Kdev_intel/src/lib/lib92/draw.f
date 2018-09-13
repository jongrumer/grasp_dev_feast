************************************************************************
*                                                                      *
      SUBROUTINE DRAW (P,SP,Q,SQ,MF)
*                                                                      *
*   This  subroutine  generates  a  printer plot. P and Q are radial   *
*   functions with the maximum tabulation point MF. SP is the factor   *
*   by which P is to be scaled, SQ is the factor by which Q is to be   *
*   scaled.                                                            *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 10 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
      CHARACTER*132 CBLANK,CDASH,CPLOT
      CHARACTER*1 CP,CQ
      LOGICAL FIRST
*
      DIMENSION P(NNNP),Q(NNNP)
      DIMENSION CPLOT(60)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
*
      DATA SIXTY / 60.0D 00/,
     :     FIFT4 / 54.0D 00/,
     :     OHTHT /131.0D 00/
*
      DATA FIRST /.TRUE./
*
*   This initialization is carried out once per run
*
      IF (FIRST) THEN
         CBLANK(1:1) = '|'
         CDASH(1:1) = '-'
         DO 1 I = 2,131
            CBLANK(I:I) = ' '
            CDASH(I:I) = '-'
    1    CONTINUE
         CBLANK(132:132) = '|'
         CDASH(132:132) = '-'
         FIRST = .FALSE.
      ENDIF
*
*   Initialization
*
      IF (SQ .EQ. 0.0D 00) THEN
         CP = 'X'
         CQ = '-'
      ELSE
         CP = 'P'
         CQ = 'Q'
      ENDIF
*
*   Determine the range of amplitude
*
      DMX = MAX (SP*P(1),SQ*Q(1))
      DMN = MIN (SP*P(1),SQ*Q(1))
      DO 2 I = 2,MF
         SPI = SP*P(I)
         SQI = SQ*Q(I)
         DMX = MAX (DMX,SPI,SQI)
         DMN = MIN (DMN,SPI,SQI)
    2 CONTINUE
*
*   Determine the radial extent of the function
*
      RMX = R(MF)
*
*   Determine the scale factors
*
      IF ((DMX .EQ. 0.0D 00) .AND.
     :    (DMN .EQ. 0.0D 00)) THEN
         WRITE (99,300)
         RETURN
      ELSE
         XSCAL = FIFT4/ABS (DMX-DMN)
      ENDIF
      YSCAL = OHTHT/DBLE (N)
*
*   Locate x = 0 if this is in the range
*
      IF ((DMX .GT. 0.0D 00) .AND.
     :    (DMN .LT. 0.0D 00)) THEN
         IOFFST = 4+XSCAL*ABS (DMN)
      ELSE
         IOFFST = 1
      ENDIF
*
*   Initialize the array CPLOT
*
      DO 3 I = 1,60
         IF ((I .EQ. 1) .OR. (I .EQ. IOFFST) .OR. (I .EQ. 60)) THEN
            CPLOT(I) = CDASH
         ELSE
            CPLOT(I) = CBLANK
         ENDIF
    3 CONTINUE
*
*   Generate the plot
*
*   Note that 'P' is the dominant character
*
      SQX = SQ*XSCAL
      SPX = SP*XSCAL
      DO 4 I = 1,MF
         IXLOCQ = NINT (SQX*Q(I))+IOFFST
         IXLOCP = NINT (SPX*P(I))+IOFFST
         IYLOC  = MAX (1,NINT (DBLE (I)*YSCAL))
         CPLOT(IXLOCQ)(IYLOC:IYLOC) = CQ
         CPLOT(IXLOCP)(IYLOC:IYLOC) = CP
    4 CONTINUE
*
*   Print plot
*
      WRITE (99,301)
      DO 5 I = 60,1,-1
         WRITE (99,302) CPLOT(I)
    5 CONTINUE
*
*   Plot-information line
*
      WRITE (99,303) RMX,DMX
*
      RETURN
*
  300 FORMAT (' No plot: function is identically zero')
  301 FORMAT ('1')
  302 FORMAT (1X,1A132)
  303 FORMAT (/50X,1P,' r (max) = ',1D10.3,'Bohr radii,'
     :            ,' Maximum of functions = ',1D10.3)
*
      END
