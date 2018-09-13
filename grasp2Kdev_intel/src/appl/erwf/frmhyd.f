************************************************************************
*                                                                      *
      SUBROUTINE FRMHYD (INDEX,NSUBS,modify)
*                                                                      *
*   This  subroutine  is  used  to  produce  estimates  of  the wave   *
*   functions from the hydrogenic approximation. The screening cons-   *
*   tant is taken to be 0.4 (Changed, see comments below,XHH)          *
*                                                                      *
*   Call(s) to: [ERWF]: DCBSRW.                                        *
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

      LOGICAL YES, GETYN
      CHARACTER*256 SOURCE
*
      DIMENSION INDEX (NNNW)
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      LOGICAL SET
      CHARACTER*2 NH
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /LEFT/SET(NNNW)
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WHFROM/SOURCE(NNNW)
     &      /hydpar/sigma(NNNW)

      COMMON/iounit/istdi,istdo,istde
      COMMON/DEFAULT/NDEF

      LOGICAL modify
*
      WRITE(istde,*)
      WRITE(istde,*) '***** Screening parameters ******'
      DO 1 J = 1, NSUBS
         LOC = INDEX(J)
         IF (.NOT. SET(LOC)) THEN
            
            WRITE (istde, '(I2,A2,F10.2)') np(loc), nh(loc), sigma(loc)
            IF (modify) THEN
               WRITE (istde, *) 'Input new value >'
               READ *, sigma(loc)
            ENDIF

            ZEFF = Z - sigma(loc)
*P fix for negative zeff
            if(ZEFF <= 0.0) ZEFF = 1.0d0

*           ...Calculate radial wavefunctions
            CALL DCBSRW (NP(LOC),NAK(LOC),ZEFF,E(LOC),PZ(LOC),
     :                   PF(1,LOC),QF(1,LOC),MF(LOC))
            SET(LOC) = .TRUE.
            !SOURCE(LOC) = 'Screened hydrogenic estimate'
            SOURCE(LOC) = 'Hyd'
         ENDIF
    1 CONTINUE
      RETURN
      END
