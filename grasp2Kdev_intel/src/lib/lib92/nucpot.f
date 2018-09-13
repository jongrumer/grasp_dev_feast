************************************************************************
*                                                                      *
      SUBROUTINE NUCPOT
*                                                                      *
*   Evaluate the nuclear potential for point and Fermi models.         *
*                                                                      *
*   Call(s) to: [LIB92] ES.                                            *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 05 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
      LOGICAL LDBPR,SET
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF9/CVAC,PI
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
*
*   Point nucleus
*
      IF (NPARM .EQ. 0) THEN
*
         DO 1 I = 1,N
            ZZ(I) = Z
    1    CONTINUE
*
*   Fermi distribution
*
      ELSEIF (NPARM. EQ. 2) THEN
*
         C = PARM(1)
         A = PARM(2)
         ABC = A/C
         TABC = 2.0D 00*ABC
         ABC2 = ABC*ABC
         THABC2 = 3.0D 00*ABC2
         ABC3 = ABC2*ABC
         CBA = C/A
         PI2 = PI*PI
         HPIAC2 = 0.5D 00*PI2*ABC2
         H3PHP = 1.5D 00+HPIAC2
         CALL ES (-CBA,S2MCBA,S3MCBA)
         SABC3 = 6.0D 00*ABC3
         DMSAS = -SABC3*S3MCBA
         EN = 1.0D 00 + ABC2*PI2 + DMSAS
         ZBN = Z/EN
*
         SET = .FALSE.
         DO 2 I = 1,N
            RI = R(I)
            RMC = RI-C
            RMCBA = RMC/A
            RBC = RI/C
            IF (RBC .LE. 1.0D 00) THEN
               CALL ES (RMCBA,S2RCBA,S3RCBA)
               ZZ(I) = ZBN*( DMSAS + SABC3*S3RCBA
     :                      +RBC*( H3PHP-THABC2*S2RCBA
     :                            -0.5D 00*RBC*RBC) )
            ELSE
               IF (.NOT. SET) THEN
                  NNUC = I
                  SET = .TRUE.
               ENDIF
               CALL ES (-RMCBA,S2RCBA,S3RCBA)
               ZZ(I) = Z * ( 1.0D 00
     :                        +THABC2 * ( RBC *S2RCBA
     :                                   +TABC*S3RCBA ) / EN )
            ENDIF
    2    CONTINUE
      ENDIF
*
      IF (LDBPR(2)) THEN
         WRITE (99,300)
         NB3 = N/3
         IF (3*NB3 .EQ. N) THEN
            NROWS = NB3
         ELSE
            NROWS = NB3+1
         ENDIF
         DO 4 II = 1,NROWS
            II1 = II
            II2 = II1+NROWS
            II3 = II2+NROWS
            IF (II3 .LE. N) THEN
               WRITE (99,301) R(II1),ZZ(II1),R(II2),ZZ(II2),
     :                        R(II3),ZZ(II3)
            ELSEIF (II2 .LE. N) THEN
               WRITE (99,301) R(II1),ZZ(II1),R(II2),ZZ(II2)
            ELSE
               WRITE (99,301) R(II1),ZZ(II1)
            ENDIF
    4    CONTINUE
      ENDIF
*
      RETURN
*
  300 FORMAT (/'From SUBROUTINE NUCPOT:'
     :        /3(' -------- r -------- ----- -r*V(r) -----'))
  301 FORMAT (1P,6(1X,1D19.12))
*
      END
