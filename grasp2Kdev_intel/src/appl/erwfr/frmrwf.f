************************************************************************
*                                                                      *
      SUBROUTINE FRMRWF (INDEX,NSUBS,FILNAM)
*                                                                      *
*   This subroutine loads  radial wavefunctions from the  .rwf  file   *
*   and performs some related setup.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, OPENFL.                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL FOUND,LDBPR,SET
      CHARACTER*256 SOURCE, FILNAM*(*)
      !CHARACTER*11 FORM
      !CHARACTER*3 STATUS
      CHARACTER*6 G92RWF
      CHARACTER*2 NH
*
      DIMENSION INDEX(NNNW)
*
      POINTER (PNTRPA,PA(*))
      POINTER (PNTRQA,QA(*))
      POINTER (PNTRRA,RA(*))
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /LEFT/SET(NNNW)
     :      /ORB1/E(NNNW),GAMA(NNNW)
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WHFROM/SOURCE(NNNW)

      COMMON/iounit/istdi,istdo,istde


      !FORM = 'UNFORMATTED'
      !STATUS = 'OLD'

      CALL OPENFL (23,FILNAM,'UNFORMATTED','OLD',IERR)
      IF (IERR .EQ. 1) THEN
         WRITE (istde,*) 'Error openning file "'
     &                   ,FILNAM(1:LEN_TRIM (FILNAM)), '"'
         CLOSE (23)
         STOP
      ENDIF
*
*   Check the file; if not as expected, try again
* 
      READ (23,IOSTAT = IOS) G92RWF
      IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
         WRITE (istde,*) 'This is not a Radial WaveFunction File;'
         CLOSE (23)
         STOP
      ENDIF
*
*   Read orbital information from Read Orbitals File; write summary
*   to  .dbg  file if option set
*
      IF (LDBPR(3)) WRITE (99,300)
    2 FOUND = .FALSE.
      READ (23,IOSTAT = IOS) NPY,NAKY,EY,MY
      IF (IOS .EQ. 0) THEN
         DO 3 J = 1,NSUBS
            LOC = INDEX(J)
            IF ((.NOT. SET(LOC)) .AND.
     :          (NP(LOC) .EQ. NPY) .AND.
     :          (NAK(LOC) .EQ. NAKY)) THEN
               FOUND = .TRUE.
               E(LOC) = EY
               CALL ALLOC (PNTRPA,MY,8)
               CALL ALLOC (PNTRQA,MY,8)
               CALL ALLOC (PNTRRA,MY,8)
               READ (23) PZ(LOC),(PA(I),I = 1,MY),(QA(I),I = 1,MY)
               READ (23) (RA(I),I = 1,MY)
               CALL INTRPQ (PA,QA,MY,RA,LOC,DNORM)
               IF (LDBPR(3)) WRITE (99,301)
     :            NP(LOC),NH(LOC),E(LOC),DNORM
               CALL DALLOC (PNTRPA)
               CALL DALLOC (PNTRQA)
               CALL DALLOC (PNTRRA)
               LENTH = LEN_TRIM(FILNAM)
               SET(LOC) = .TRUE.
               !SOURCE(LOC)(1:LENTH) = FILNAM(1:LENTH)
               SOURCE(LOC) = filnam(1:3)
               GOTO 2
            ENDIF
    3    CONTINUE
         IF (.NOT. FOUND) THEN
            READ (23)
            READ (23)
            GOTO 2
         ENDIF
      ENDIF
      IF (LDBPR(3)) WRITE (99,*) ' orbitals renormalised;'
*
      CLOSE (23)
*
      RETURN
*
  300 FORMAT (/'From SUBROUTINE FRMRWF:'
     :        /' Orbital',8X,'Eigenvalue',19X,'Norm')
  301 FORMAT (2X,I2,A2,4X,1P,1D22.15,4X,1D22.15)
*
      END
