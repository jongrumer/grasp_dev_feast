************************************************************************
*                                                                      *
      SUBROUTINE TIINIG(CIIN,NCSF,NCIV,I,L,CONST,CIOUT,NTESTG)
*                                                                      *
*   Calculates the action of the operator                              *
*   Const ** E(li,li) on a set of vectors                              *
*                                                                      *
*   Adapted for GRASP, daughter of TIINI, born February 1996           *
*                                                                      *
*   =====                                                              *
*   Input                                                              *
*   =====                                                              *
*   CIIN      : Input CI vectors                                       *
*   NCSCF     : Length of CI expansion                                 *
*   NCIV      : Number of CI vectors                                   *
*   I         : Shell number                                           *
*   L         : symmetry                                               *
*   NSHLP(L,K): Gives the shell number as defined in getcsl for        *
*               the K:th shell with symmetry L                         *
*   CONST     : The constant                                           *
*                                                                      *
*   ======                                                             *
*   Output                                                             *
*   ======                                                             *
*                                                                      *
*   CIOUT : List of output CI vectors                                  *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'parameters.def'
CGG      PARAMETER (KEY = 121)
CGG begin
      PARAMETER (KEY = KEYORB)
CGG end
      PARAMETER (NLMAX = 20)                                            
CGG      PARAMETER (NNNW = 120)

      DIMENSION CIIN(NCSF,NCIV)
      DIMENSION CIOUT(NCSF,NCIV)

      POINTER (PJANN,JANN(*)),(PPPINT,INTGRL(*)),
     :        (PJBNN,JBNN(*)),(PCNN,CNN(*)),(PINTPTR,INTPTR(*))

      COMMON/MCPDATA/PJANN,PJBNN,PPPINT,PCNN,PINTPTR,NCOEFF,NINTG
      COMMON/SBDAT1/NSHLP(NLMAX,NLMAX),NSHLPP(NLMAX,NNNW)
*
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      NTEST = 000

      IF(NTEST.GE.10) WRITE(6,*) ' Entering TIINI'

      CALL SETVEC(CIOUT,0.0D0,NCSF*NCIV)
*
*.  Obtain address of first coupling coefficient for h(il,il) :  IFIRST
*.  Obtain number of        coupling coefficient for h(il,il) :  NFOUND
*.
*     NFOUND : Number of coefficients obtained
*     IVAL   : actual RACAH coefficient <CSF(L)!E(il,il)!CSF(R)>
*     ILEFT  = CSF(L) ?
*
      NFOUND = 0
      DO K = 1,NINTG
        IA = INTGRL(K)/KEY
        IB = MOD(INTGRL(K),KEY)
        IF (NSHLP(L,I).EQ.IA.AND.IA.EQ.IB) THEN
          IF (K.EQ.1) THEN
            IFIRST = 1
          ELSE
            IFIRST = INTPTR(K-1) + 1
          ENDIF
          NFOUND = INTPTR(K) - IFIRST + 1
          GOTO 15
        ENDIF
      ENDDO

   15 DO IELMNT = 1, NFOUND
        IVAL = CNN(IFIRST-1+IELMNT)
        CONSTN = CONST ** IVAL
        ILEFT = JANN(IFIRST-1+IELMNT)
        DO IVEC=1, NCIV
          CIOUT(ILEFT,IVEC) = CONSTN * CIIN(ILEFT,IVEC)
        ENDDO
      ENDDO
*
*.  The previous provided us with all
*   terms with nonvanishing occupation.
*   For terms with vanishing occupation of il,
*   just copy coefficients, since (x) ** 0 = 1
*
      DO IVEC = 1, NCIV
        DO IELMNT = 1, NCSF
          IF(CIOUT(IELMNT,IVEC).EQ.0.0D0)
     &     CIOUT(IELMNT,IVEC) = CIIN(IELMNT,IVEC)
        ENDDO
      ENDDO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input and output vectors from TIINI I,L',I,L
        CALL WRTMAT(CIIN,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
        CALL WRTMAT(CIOUT,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
      END IF
*
      IF(NTEST.GE.10) WRITE(6,*) ' Leaving  TIINI'

      RETURN
      END
