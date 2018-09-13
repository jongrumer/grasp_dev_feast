************************************************************************
*                                                                      *
      SUBROUTINE LODRWFF(NAME)
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
CFF      PARAMETER (NNNP = 590) 
CFF      PARAMETER (NNN1 = NNNP+10)
CFF      PARAMETER (NNNW = 120)

      CHARACTER*2 NHFF
      CHARACTER*24 NAME
      CHARACTER*6 G92RWF
*
      POINTER (PNTRPA,PA(1)),(PNTRQA,QA(1)),(PNTRRA,RA(1))
*
      POINTER (PNTRPFFF,PFFF(NNNP,1)),(PNTRQFFF,QFFF(NNNP,1))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
*
*   Common relevant for the final state
*
      COMMON/ORB1FF/EFF(NNNW),GAMAFF(NNNW)
     :      /ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)
*
*   Write entry message
*
      PRINT *, 'Loading Radial WaveFunction File for final state...'
*
*   Open the radial wave function file
*
      J = INDEX(NAME,' ')
      OPEN (UNIT = 69,FILE=NAME(1:J-1)//'.bw',FORM='UNFORMATTED',
     :     STATUS='OLD')
*
*   Allocate storage to orbital arrays
*
      CALL ALLOC (PNTRPFFF,NNNP*NWFF,8)
      CALL ALLOC (PNTRQFFF,NNNP*NWFF,8)
*
*   Setup: (1) Orbital arrays to zero
*          (2) Array E to -1 (no orbitals estimated)
*          (3) Parameters GAMMA for each orbital
*
      CON = Z/C
      CON = CON*CON
*
      DO J = 1,NWFF
        DO I = 1,NNNP
          PFFF(I,J) = 0.0D 00
          QFFF(I,J) = 0.0D 00
        ENDDO
*
        EFF(J) = -1.0D 00
*
        K = ABS (NAKFF(J))
        IF (NPARM .GT. 0) THEN
          GAMAFF(J) = DBLE (K)
        ELSEIF (NPARM .EQ. 0) THEN
          FKK = DBLE (K*K)
          IF (FKK .GE. CON) THEN
            GAMAFF(J) = SQRT (FKK-CON)
          ELSE
            PRINT *, 'LODRWF: Imaginary gamma parameter'
            PRINT *, ' for ',NPFF(J),NHFF(J),' orbital; the'
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

    3 READ (69,IOSTAT = IOS) NPYFF,NAKYFF,EYFF,MYFF
      IF (IOS .EQ. 0) THEN
        CALL ALLOC (PNTRPA,MYFF,8)
        CALL ALLOC (PNTRQA,MYFF,8)
        CALL ALLOC (PNTRRA,MYFF,8)
        READ (69) PZY,(PA(I),I = 1,MYFF),(QA(I),I = 1,MYFF)
        READ (69) (RA(I),I = 1,MYFF)

        DO J = 1,NWFF
          IF ( (EFF(J) .LT. 0.0D 00) .AND.
     :         (NPYFF .EQ. NPFF(J)) .AND.
     :         (NAKYFF .EQ. NAKFF(J)) ) THEN
            PZFF(J) = PZY
            EFF(J) = EYFF
            MFFF(J) = MYFF
            DO JJ = 1,MFFF(J)
               PFFF(JJ,J) = PA(JJ)
               QFFF(JJ,J) = QA(JJ)
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
      IF (NWIN .LT. NWFF) THEN
        PRINT *, 'LODRWF: All required orbitals not'
        PRINT *, ' found.'
        STOP
      ENDIF
*
      PRINT *, ' ... load complete;'
*
      CLOSE(69)

      RETURN
      END
