************************************************************************
*                                                                      *
      SUBROUTINE GENRWF
*                                                                      *
*   Controls the computation of the subshell radial wavefunctions.     *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, GETRSL, GETYN.                         *
*               [ERWF]: FRMHYD, FRMRWF, FRMTFP, PRTREM, SUMMRY,        *
*                       TFPOT.                                         *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*   Updated by Xinghong He                Last revision: 11 Nov 1997   *
*    for menu-driven and less questions
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL ALL,GETYN,SET,YES
      CHARACTER*256 SOURCE, infile*128
      CHARACTER*2 NH
*
      DIMENSION INDEX (NNNW)
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
cff
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/COUN/THRESH
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEFAULT/NDEF
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /LEFT/SET(NNNW)
     :      /NPAR/PARM(2),NPARM
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WHFROM/SOURCE(NNNW)
      COMMON/iounit/istdi,istdo,istde
      LOGICAL modify, existed
      DATA modify/.FALSE./
*
*   Set the threshold for node counting
*
      THRESH = 0.05
*
*   Set up the Thomas-Fermi potential
*
      CALL TFPOT
*
*   Allocate storage to the arrays that store the subshell
*   radial wavefunction arrays; initialise these and all
*   associated arrays
*
      CALL ALLOC (PNTRPF,NNNP*NW,8)
      CALL ALLOC (PNTRQF,NNNP*NW,8)

      CON = Z/C
      CON = CON*CON

      DO J = 1, NW
         SET(J) = .FALSE.
         SOURCE(J) = ' '
         DO I = 1,N
            PF(I,J) = 0.d0
            QF(I,J) = 0.d0
         ENDDO

         K = ABS (NAK(J))
         IF (NPARM .GT. 0) THEN
            GAMA(J) = DBLE (K)
         ELSEIF (NPARM .EQ. 0) THEN
            FKK = DBLE (K*K)
            IF (FKK .GE. CON) THEN
               GAMA(J) = SQRT (FKK-CON)
            ELSE
               WRITE (istde,*) 'LODRWF: Imaginary gamma parameter '
     &                      , 'for ',NP(J),NH(J),' orbital;'
               WRITE (istde,*) 'the point model for the nucleus '
     &                      , 'is inappropriate for Z > ',C,'.'
               STOP
            ENDIF
         ENDIF
       ENDDO
*
*   Write out the complete list of subshell radial wave functions
*
      CALL PRTREM (ALL)
*
* Direct to read radial functions till finish
*
  123 CONTINUE
      IF (.NOT. ALL) THEN
  234    WRITE (istde,*) 
         WRITE (istde,*) 'Read subshell radial wavefunctions. ',
     &                   'Choose one below'
         WRITE (istde,*) '    1 -- GRASP2K File'
         WRITE (istde,*) '    2 -- Thomas-Fermi'
         WRITE (istde,*) '    3 -- Screened Hydrogenic'

         READ (istdi,*) nradial
         IF (nradial .LT. 1 .OR. nradial .GT.3) THEN
            WRITE (istde,*) nradial, 'is not a valid choice, redo'
            GOTO 234
         ENDIF

         IF (nradial .EQ. 1) THEN
  345       WRITE (istde,*) 'Enter the file name (Null then "rwfn.out")'
            READ (istdi, '(A)') infile
            IF (LEN_TRIM (infile) .EQ. 0) infile = 'rwfn.out'

            INQUIRE (FILE = infile, EXIST = existed)
            IF (.NOT. existed) THEN 
               WRITE (istde,*) ' File "', infile(1:LEN_TRIM (infile))
     &                         , '" does not exist, redo'
               GOTO 345
            ENDIF
         ENDIF

         WRITE(istde,*) 'Enter the list of relativistic subshells:'
         CALL GETRSL (INDEX,NSUBS)

         IF (nradial .EQ. 1) THEN
            CALL FRMRWF (INDEX,NSUBS, infile)
         ELSEIF (nradial .EQ. 2) THEN
            CALL FRMTFP (INDEX,NSUBS)
         ELSE
            CALL FRMHYD (INDEX,NSUBS,modify)
         ENDIF

         CALL PRTREM (ALL)
         IF (.NOT. ALL) THEN
            !WRITE (istde,*) 'Radial functions incomplete, need more...'
            ! PRTREM has more informative prompt
            GOTO 234
         ENDIF

      ENDIF
*
* All read. Let know, and allow modifying if non default
*
      WRITE (istde,*) 'All required subshell radial wavefunctions '
     :       , ' have been estimated:'
      CALL SUMMRY (istde)

      IF (ndef .EQ. 0) THEN
         modify = .FALSE.
      ELSE
         WRITE (istde,*) 'Revise any of these estimates?'
         modify = GETYN ()
      ENDIF

      IF (modify) THEN
  456    WRITE (istde,*) 'Enter the list of subshells whose radial ',
     &                  'wavefunctions are to be revised:'
         CALL GETRSL (INDEX,NSUBS)
         IF (NSUBS .EQ. 0) GOTO 456
         DO J = 1,NSUBS
            LOC = INDEX(J)
            SET(LOC) = .FALSE.
            DO I = 1,N
               PF(I,LOC) = 0.d0
               QF(I,LOC) = 0.d0
            ENDDO
            SOURCE(LOC) = ' '
         ENDDO
         ALL = .FALSE.
         GOTO 123
      ENDIF

      IF (NDEF.NE.0) CALL SUMMRY (24)

      RETURN
      END
