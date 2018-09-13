************************************************************************
*                                                                      *
      SUBROUTINE NCHARG
*                                                                      *
*   This routine evaluates the nuclear charge density, and stores it   *
*   in the  common  array  ZDIST .                                     *
*                                                                      *
*   Call(s) to: [LIB92]: ES.                                           *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
      LOGICAL FORM1,FORM2
*
      COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF9/CVAC,PI
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NCDIST/ZDIST(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
*   Initialize array to zero
*
      DO 1 I = 1,N
         ZDIST(I) = 0.0D 00
    1 CONTINUE
*
*   Fermi charge distribution
*
      IF (NPARM .EQ. 2) THEN
         C = PARM(1)
         A = PARM(2)
         CBA = C/A
         PI2 = PI*PI
         ABC = A/C
         ABC2 = ABC*ABC
         ABC3 = ABC2*ABC
         CALL ES (-CBA,S2MCBA,S3MCBA)
         EN = 1.0D 00+PI2*ABC2-6.0D 00*ABC3*S3MCBA
         ZNORM = 3.0D 00*Z/(4.0D 00*PI*EN*C**3)
         FORM1 = .TRUE.
         FORM2 = .FALSE.
         DO 2 I = 1,N
            IF     (FORM1) THEN
               EXTRM = EXP ( (R(I)-C)/A )
               ZDIST(I) = ZNORM/( 1.0D 00+EXTRM )
               IF (1.0D 00/EXTRM .LE. PRECIS) THEN
                  FORM1 = .FALSE.
                  FORM2 = .TRUE.
               ENDIF
            ELSEIF (FORM2) THEN
               ZDISTI = ZNORM*EXP ( -(R(I)-C)/A )
               IF (ABS (ZDISTI) .GT. 0.0D 00) THEN
                  ZDIST(I) = ZDISTI
               ELSE
                  MTP = I
                  GOTO 3
               ENDIF
            ENDIF
    2    CONTINUE
      ENDIF
*
    3 RETURN
*
      END
