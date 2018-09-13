************************************************************************
*                                                                      *
      SUBROUTINE FRMTFP (INDEX,NSUBS)
*                                                                      *
*   This  subroutine  is  used  to  produce  estimates  of  the wave   *
*   functions by use of the Thomas-Fermi approximation to the direct   *
*   potential. SUBROUTINE  SOLVH  is used  to obtain the radial wave   *
*   functions.                                                         *
*                                                                      *
*   Call(s) to: [ERWF]: SOLVH.                                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL FAIL,SET
      CHARACTER*256 SOURCE
      CHARACTER*2 NH
*
      DIMENSION INDEX (NNNW)
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/LEFT/SET(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WHFROM/SOURCE(NNNW)
*
      DO 1 J = 1,NSUBS
*
         LOC = INDEX(J)
*
         IF (.NOT. SET(LOC)) THEN
*
*   Estimate the leading coefficient of the series expansion for
*   the orbital
*
            PZ(LOC) = 10.0D 00
*
*   Calculate radial wavefunctions
*
            CALL SOLVH (LOC,FAIL)
*
*   Message if SOLVH did not converge; reset SET(LOC) otherwise
*
            IF (FAIL) THEN
               WRITE (*,300) NP(LOC),NH(LOC)
            ELSE
               SET(LOC) = .TRUE.
               !SOURCE(LOC) = 'Thomas-Fermi estimate'
               SOURCE(LOC) = 'T-F'
            ENDIF
*
         ENDIF
*
    1 CONTINUE
*
      RETURN
*
  300 FORMAT (/'TFWAVE: Unable to compute radial'
     :       //' wavefunction for ',I2,A2,' subshell;')
*
      END
