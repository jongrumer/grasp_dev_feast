************************************************************************
*                                                                      *
      SUBROUTINE GETMIXCZ(NAME)
*                                                                      *
*   Open, check, load data from and close the  .mix  file. This file   *
*   is always attached to stream 25.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC,  OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 25 Dec 1992   *
*                                                                      *
************************************************************************
*

! Some PRINT *, changed to WRITE(istde,*); lines joined
!XHH 1997.02.05

      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL FOUND
      CHARACTER*256 FILNAM
      CHARACTER*24 NAME
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
      K = INDEX(NAME,' ')
      FILNAM = NAME(1:K-1)//'.cm'
      FORM = 'UNFORMATTED'
      STATUS = 'OLD'
*
      CALL OPENFL (25,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) THEN
         WRITE(istde,*) 'Error when opening',FILNAM
         STOP
      ENDIF
*
*   Check the header of the file; if not as expected, try again
*
      READ (25,IOSTAT = IOS) G92MIX
      IF ((IOS .NE. 0) .OR.
     :    (G92MIX .NE. 'G92MIX')) THEN
         WRITE(istde,*) 'Not a GRASP92 MIXing Coefficients File;'
         CLOSE (25)
         STOP
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
         STOP
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
      nncf=ncf/10000+1
      do j = 1,nvec
      do ii=1,nncf
      ncf0=1+(ii-1)*10000
      ncf1=ii*10000
      if(ncf1.gt.ncf) ncf1=ncf
      READ (25) (EVEC(I+(J-1)*NCF),I = NCF0,NCF1)
c     READ (25,IOSTAT = IOS) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
      enddo
      enddo

*
      PRINT *, 'IOS= ',ios,' ... load complete;'
*
*   Close the  .mix  file
*
      CLOSE (25)
*
      RETURN
      END
