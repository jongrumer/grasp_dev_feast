************************************************************************
*                                                                      *
      SUBROUTINE LODRWF (IiERR)
*                                                                      *
*   This subroutine loads  radial wavefunctions from the  .rwf  file   *
*   and performs some related setup.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
*   Block version by Xinghong He          Last revision: 27 May 1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      POINTER (PNTRPA,PA(*))
      POINTER (PNTRQA,QA(*))
      POINTER (PNTRRA,RA(*))
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)

!-----------------------------------------------------------------------
*
*   Write entry message
*
      PRINT *, 'Loading Radial WaveFunction File ...'
*
*   Allocate storage to orbital arrays
*
      CALL ALLOC (PNTRPF,NNNP*NW,8)
      CALL ALLOC (PNTRQF,NNNP*NW,8)
*
*   Setup: (1) Orbital arrays to zero
*          (2) Array E to -1 (no orbitals estimated)
*          (3) Parameters GAMMA for each orbital
*
      CON = Z/C
      CON = CON*CON
*
      DO 2 J = 1,NW
*
         DO 1 I = 1,N
            PF(I,J) = 0.0D 00
            QF(I,J) = 0.0D 00
    1    CONTINUE
*
         E(J) = -1.0D 00
*
         K = ABS (NAK(J))
         IF (NPARM .GT. 0) THEN
            GAMA(J) = DBLE (K)
         ELSEIF (NPARM .EQ. 0) THEN
            FKK = DBLE (K*K)
            IF (FKK .GE. CON) THEN
               GAMA(J) = SQRT (FKK-CON)
            ELSE
               !WRITE (istde,*) 'LODRWF: Imaginary gamma parameter'
               !WRITE (istde,*) ' for ',NP(J),NH(J),' orbital; the'
               !WRITE (istde,*) ' point model for the nucleus'
               !WRITE (istde,*) ' is inappropriate for Z > ',C,'.'
               STOP 'lodrwf: Inappropriate gamma'
            ENDIF
         ENDIF
*
    2 CONTINUE
*
*   Read orbital information from Read Orbitals File; write summary
*   to  .dbg  file if option set
*
      IF (LDBPR(3)) WRITE (99,300)
      NWIN = 0
    3 CONTINUE

      READ (23,IOSTAT = IOS) NPY,NAKY,EY,MY

      IF (IOS .EQ. 0) THEN
         CALL ALLOC (PNTRPA,MY,8)
         CALL ALLOC (PNTRQA,MY,8)
         CALL ALLOC (PNTRRA,MY,8)

         READ (23) PZY,(PA(I),I = 1,MY),(QA(I),I =1 ,MY)
         READ (23) (RA(I),I = 1,MY)

         DO 4 J = 1,NW
            IF ( (E(J) .LT. 0.0D 00) .AND.
     :           (NPY .EQ. NP(J)) .AND. (NAKY .EQ. NAK(J)) ) THEN
               PZ(J) = PZY
               E(J) = EY
               CALL INTRPQ (PA,QA,MY,RA,J,DNORM)
               IF (LDBPR(3)) WRITE (99,301) NP(J),NH(J),E(J),DNORM
               NWIN = NWIN+1
            ENDIF
    4    CONTINUE
         CALL DALLOC (PNTRPA)
         CALL DALLOC (PNTRQA)
         CALL DALLOC (PNTRRA)
         GOTO 3
      ENDIF
      IF (LDBPR(3)) WRITE (99,*) ' orbitals renormalised;'
*
*   Stop with an error message if all orbitals are not known
*
      IF (NWIN .LT. NW) THEN
         IiERR = 1
         GOTO 5
      ENDIF
*
*   Schmidt orthogonalise the orbitals
*
      CALL ORTHSC
      IF (LDBPR(3)) WRITE (99,*) ' orbitals orthogonalised'
     :                         //' and renormalised;'
*
      IiERR = 0
    5 RETURN
*
  300 FORMAT (/'From SUBROUTINE LODRWF:'
     :        /' Orbital',8X,'Eigenvalue',19X,'Norm')
  301 FORMAT (2X,I2,A2,4X,1P,1D22.15,4X,1D22.15)
*
      END
