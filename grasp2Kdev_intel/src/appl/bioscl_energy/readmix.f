************************************************************************
*                                                                      *
      SUBROUTINE READMIX(NAME,INPCI,INIT)
*                                                                      *
*   Open and read the mixing coefficent files                          *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

Cww      INTEGER PNTRIQ, PNTIQR
      POINTER(PNTRIQ,RIQDUMMY)
      POINTER(PNTIQR,IQRDUMMY)

      CHARACTER*24 NAME(2)
      CHARACTER*6 G92MIX

      POINTER(PNEVALII,EVALII(1))
      POINTER(PNEVECII,EVECII(1))
      POINTER(PNIVECII,IVECII(1))
      POINTER(PIATJPII,IATJPOII(1)),(PIASPAII,IASPARII(1))

      POINTER(PNEVALFF,EVALFF(1))
      POINTER(PNEVECFF,EVECFF(1))
      POINTER(PNIVECFF,IVECFF(1))
      POINTER(PIATJPFF,IATJPOFF(1)),(PIASPAFF,IASPARFF(1))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF1R/NELECR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB2R/NCFR,NWR,PNTIQR

      COMMON/DEF1II/EMNII,IONCTYII,NELECII,ZII
     :      /EIGVALII/EAVII,PNEVALII
     :      /EIGVECII/PNEVECII
     :      /ORB2II/NCFII,NWII
     :      /PRNTII/NVECII,PNIVECII,NVECMXII
     :      /SYMAII/PIATJPII,PIASPAII

      COMMON/DEF1FF/EMNFF,IONCTYFF,NELECFF,ZFF
     :      /EIGVALFF/EAVFF,PNEVALFF
     :      /EIGVECFF/PNEVECFF
     :      /ORB2FF/NCFFF,NWFF
     :      /PRNTFF/NVECFF,PNIVECFF,NVECMXFF
     :      /SYMAFF/PIATJPFF,PIASPAFF

      if (INIT.eq.1) then
*
*   Read the initial state mixing file    
*
C      J = INDEX(NAME(1),' ')
C      IF (INPCI.EQ.0) THEN
C        OPEN (UNIT = 68,FILE=NAME(1)(1:J-1)//'.cbm',FORM='UNFORMATTED',
C     :        STATUS='OLD')
C      ELSE
C        OPEN (UNIT = 68,FILE=NAME(1)(1:J-1)//'.bm',FORM='UNFORMATTED',
C     :        STATUS='OLD')
C      ENDIF
C      READ(68,IOSTAT=IOS) G92MIX
C      IF ((IOS.NE.0).OR.(G92MIX.NE.'G92MIX')) THEN
C        WRITE(*,*) 'Not a GRASP mixing file'
C        STOP
C      ENDIF

C      READ(68) N11,N12,N13
C      IF ((N11.NE.NELECII).OR.(N12.NE.NCFII).OR.(N13.NE.NWII)) THEN
C        PRINT *, 'This MIXing Coefficients File is not'
C        PRINT *, 'appropriate for the initial state'
C        STOP
C      ENDIF
*
      READ(68) ib,ncfii,NVECII,iatjp,iaspa
      CALL ALLOC(PNEVALII,NVECII,8)
      CALL ALLOC(PNEVECII,NCFII*NVECII,8)
      CALL ALLOC(PNIVECII,NVECII,4)
      CALL ALLOC(PIATJPII,NVECII,4)
      CALL ALLOC(PIASPAII,NVECII,4)
      READ(68) (IVECII(I),I=1,NVECII)
c     READ(68) (IATJPOII(I),IASPARII(I),I=1,NVECII)
      do i=1,nvecii
        IATJPOII(I) = iatjp
        IASPARII(I) = iaspa
      enddo
      READ(68) EAVII,(EVALII(I),I=1,NVECII)

C      DO 10 I=1,NVECII
C         EVALII(I)=EAVII+EVALII(I)
C10    CONTINUE
      READ(68) ((EVECII(I+(J-1)*NCFII),I=1,NCFII),J=1,NVECII)
*
*   Close the initial state mixing  file
*
c     CLOSE(68)
*
*   Read the final state mixing file    
*
      else
*
*   Read the final state mixing file    
*

c      J = INDEX(NAME(2),' ')
c      IF (INPCI.EQ.0) THEN
c        OPEN (UNIT = 78,FILE=NAME(2)(1:J-1)//'.cbm',FORM='UNFORMATTED',
c     :        STATUS='OLD')
c      ELSE 
c        OPEN (UNIT = 78,FILE=NAME(2)(1:J-1)//'.bm',FORM='UNFORMATTED',
c     :        STATUS='OLD')
c      ENDIF

c      READ(78,IOSTAT=IOS) G92MIX
c      IF ((IOS.NE.0).OR.(G92MIX.NE.'G92MIX')) THEN
c        WRITE(*,*) 'Not a GRASP mixing file'
c        STOP
c      ENDIF

c      READ(78) N11,N12,N13
c      IF ((N11.NE.NELECFF).OR.(N12.NE.NCFFF).OR.(N13.NE.NWFF)) THEn
c        PRINT *, 'This MIXing Coefficients File is not'
c        PRINT *, 'appropriate for the final state'
c        STOP
c      ENDIF
*
c     READ(78) NVECFF
      READ(78) ib,ncfff,NVECFF,iatjp,iaspa
      CALL ALLOC(PNEVALFF,NVECFF,8)
      CALL ALLOC(PNEVECFF,NCFFF*NVECFF,8)
      CALL ALLOC(PNIVECFF,NVECFF,4)
      CALL ALLOC(PIATJPFF,NVECFF,4)
      CALL ALLOC(PIASPAFF,NVECFF,4)
      READ(78) (IVECFF(I),I=1,NVECFF)
c     READ(78) (IATJPOFF(I),IASPARFF(I),I=1,NVECFF)
      do i=1,nvecff
        IATJPOFF(I) = iatjp
        IASPARFF(I) = iaspa
      enddo
      READ(78) EAVFF,(EVALFF(I),I=1,NVECFF)

C      DO 20 I=1,NVECFF
C         EVALFF(I)=EAVFF+EVALFF(I)
C20    CONTINUE
      READ(78) ((EVECFF(I+(J-1)*NCFFF),I=1,NCFFF),J=1,NVECFF)
*
*   Close the initial state mixing  file
*
c      CLOSE(78)
      endif

      RETURN
      END
