************************************************************************
*                                                                      *
      SUBROUTINE LODRWFI(NAME)
*                                                                      *
*   This subroutine loads  radial wavefunctions from the  .rwf  file   *
*   and performs some related setup.                                   *
*                                                                      *
*   Written by Per Jonsson                              June 1996      *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)


      include 'parameters.def'
      CHARACTER*2 NHII
      CHARACTER*128 NAME
      CHARACTER*6 G92RWF
*
      POINTER (PNTRPA,PA(*)),(PNTRQA,QA(*)),(PNTRRA,RA(*))
*
      POINTER (PNTRPFII,PFII(NNNP,*)),(PNTRQFII,QFII(NNNP,*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
*
*   Common relevant for the initial state
*
      COMMON/ORB1II/EII(NNNW),GAMAII(NNNW)
     :      /ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB10II/NHII(NNNW)
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*
*   Write entry message
*
*R    PRINT *, 'Loading Radial WaveFunction File for initial state...'
*
*   Open the radial wave function file
*
      J = INDEX(NAME,' ')
      OPEN (UNIT = 69,FILE=NAME(1:J-1)//'.bw',FORM='UNFORMATTED',
     :     STATUS='OLD')
*
*   Allocate storage to orbital arrays
*
      CALL ALLOC (PNTRPFII,NNNP*NWII,8)
      CALL ALLOC (PNTRQFII,NNNP*NWII,8)
*
      CON = Z/C
      CON = CON*CON
*
*R    write(*,*) 'NWII',nwii
      DO J = 1,NWII
*R      write(*,*) NAKII(J) ,NPII(J),NHII(J)
        DO I = 1,NNNP
          PFII(I,J) = 0.0D 00
          QFII(I,J) = 0.0D 00
        ENDDO
*
        EII(J) = -1.0D 00
*
        K = ABS (NAKII(J))
        IF (NPARM .GT. 0) THEN
          GAMAII(J) = DBLE (K)
        ELSEIF (NPARM .EQ. 0) THEN
          FKK = DBLE (K*K)
          IF (FKK .GE. CON) THEN
            GAMAII(J) = SQRT (FKK-CON)
          ELSE
            PRINT *, 'LODRWF: Imaginary gamma parameter'
            PRINT *, ' for ',NPII(J),NHII(J),' orbital; the'
            PRINT *, ' point model for the nucleus'
            PRINT *, ' is inappropriate for Z > ',C,'.'
            STOP
          ENDIF
        ENDIF
*
      ENDDO
*
*   Read orbital information from Read Orbitals File;
*
      NWIN = 0
      READ (69,IOSTAT = IOS) G92RWF
      IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
        PRINT *, 'This is not a Radial WaveFunction File;'
        CLOSE (69)
      ENDIF

    3 READ (69,IOSTAT = IOS) NPYII,NAKYII,EYII,MYII
      IF (IOS .EQ. 0) THEN
        CALL ALLOC (PNTRPA,MYII,8)
        CALL ALLOC (PNTRQA,MYII,8)
        CALL ALLOC (PNTRRA,MYII,8)
        READ (69) PZY,(PA(I),I = 1,MYII),(QA(I),I = 1,MYII)
        READ (69) (RA(I),I = 1,MYII)

        DO J = 1,NWII
          IF ( (EII(J) .LT. 0.0D 00) .AND.
     :         (NPYII .EQ. NPII(J)) .AND.
     :         (NAKYII .EQ. NAKII(J)) ) THEN
            PZII(J) = PZY
            EII(J) = EYII
            MFII(J) = MYII
            DO JJ = 1,MFII(J)
               PFII(JJ,J) = PA(JJ)
               QFII(JJ,J) = QA(JJ)
            ENDDO
            NWIN = NWIN+1
          ENDIF
        ENDDO
        CALL DALLOC (PNTRPA)
        CALL DALLOC (PNTRQA)
        CALL DALLOC (PNTRRA)
        GOTO 3
      ENDIF
*
*   Stop with an error message if all orbitals are not known
*
      IF (NWIN .LT. NWII) THEN
        PRINT *, 'LODRWF: All required orbitals not'
        PRINT *, ' found.'
        STOP
      ENDIF
*
*R    PRINT *, ' ... load complete;'
*
      CLOSE(69)

      RETURN
      END
