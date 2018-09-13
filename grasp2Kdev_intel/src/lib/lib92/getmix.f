************************************************************************
*                                                                      *
      SUBROUTINE GETMIX
*                                                                      *
*   Open, check, load data from and close the  .mix  file. This file   *
*   is always attached to stream 25.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, OPENFL.                                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 25 Dec 1992   *
*                                                                      *
************************************************************************
*

! LENGTH replaced by LEN_TRIM
! Some PRINT *, changed to WRITE(istde,*); lines joined
!XHH 1997.02.05

      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL FOUND
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*6 G92MIX
      CHARACTER*3 STATUS
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*))
      POINTER (PIASPA,IASPAR(*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
      COMMON/iounit/istdi,istdo,istde
*
*   The  .mix  file is UNFORMATTED; it must exist
*
      DEFNAM = 'grasp92.mix'
      FORM = 'UNFORMATTED'
      STATUS = 'OLD'
*
*   Look for  grasp92.mix
*
      INQUIRE (FILE = DEFNAM, EXIST = FOUND)
*
      IF (FOUND) THEN
*
*   File  grasp92.csl  exists; ascertain that it is to be used;
*   if it is not to be used, determine another filename
*
         WRITE(istde,*) 'File  grasp92.mix  found; enter another file '
     &,                 'name if this file is not to be'
         WRITE(istde,*) 'used as the GRASP92 MIXing Coefficients File;'
     &,                 ' <cr> otherwise:'
         READ (*,'(A)') FILNAM
*
         IF ( LEN_TRIM(FILNAM) .EQ. 0) FILNAM = DEFNAM
*
      ELSE
*
*   File  grasp92.mix  does not exist; determine the name of the
*   .mix  file
*
    1    WRITE(istde,*) 'Enter the name of the GRASP92 MIXing '
     &,                 'Coefficients File:'
         READ (*,'(A)') FILNAM
*
         IF ( LEN_TRIM(FILNAM) .EQ. 0) GOTO 1
*
      ENDIF
!xhb
! To prevent jumping into an IF block, added the following 6
! statements and changed the subsequent "GOTO 1" to "GOTO 2"
! Mon Jan 13 15:52:01 CST 1997
      GOTO 3
    2 WRITE(istde,*) 'Enter the name of the GRASP92 MIXing '
     &,              'Coefficients File:'
      READ (*,'(A)') FILNAM
      IF ( LEN_TRIM(FILNAM) .EQ. 0) GOTO 2
    3 CONTINUE

!xhe
*
      CALL OPENFL (25,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) GOTO 2
*
*   Check the header of the file; if not as expected, try again
*
      READ (25,IOSTAT = IOS) G92MIX
      IF ((IOS .NE. 0) .OR.
     :    (G92MIX .NE. 'G92MIX')) THEN
         WRITE(istde,*) 'Not a GRASP92 MIXing Coefficients File;'
         CLOSE (25)
         GOTO 2
      ENDIF
*
      READ (25) NELECT,NCFT,NWT
      IF ((NELEC .NE. NELECT) .OR.
     :    (NCF .NE. NCFT) .OR.
     :    (NW .NE. NWT)) THEN
         WRITE(istde,*) 'This MIXing Coefficients File is not '
     &,                 'appropriate to Coefficients file is'
         WRITE(istde,*) 'not the Configuration Symmetry List File.'
         CLOSE (25)
         GOTO 2
      ENDIF
*
*   Load data from the  .mix  file
*
      PRINT *, 'Loading MIXing Coefficients File ...'
*
      READ (25) NVEC
      CALL ALLOC (PNEVAL,NVEC,8)
      CALL ALLOC (PNEVEC,NCF*NVEC,8)
      CALL ALLOC (PNIVEC,NVEC,4)
      CALL ALLOC (PIATJP,NVEC,4)
      CALL ALLOC (PIASPA,NVEC,4)
      READ (25) (IVEC(I),I = 1,NVEC)
      READ (25) (IATJPO(I),IASPAR(I),I = 1,NVEC)
      READ (25) EAV,(EVAL(I),I = 1,NVEC)
      READ (25) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
*
      PRINT *, ' ... load complete;'
*
*   Close the  .mix  file
*
      CLOSE (25)
*
      RETURN
      END
