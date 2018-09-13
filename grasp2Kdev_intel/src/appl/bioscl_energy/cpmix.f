      SUBROUTINE CPMIX(NAME,INPCI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*24 NAME(2)
      CHARACTER*6 G92MIX

      POINTER(PNEVAL,EVAL(1))
      POINTER(PNEVEC,EVEC(1))
      POINTER(PNIVEC,IVEC(1))


      NAME(2)= TRIM(NAME(1))//'_CP'
      IF (INPCI.EQ.0) THEN
        OPEN (UNIT = 78,FILE=TRIM(NAME(2))//'.cbm',FORM='UNFORMATTED',
     :        STATUS='UNKNOWN')
      ELSE
        OPEN (UNIT = 78,FILE=TRIM(NAME(2))//'.bm',FORM='UNFORMATTED',
     :        STATUS='UNKNOWN')
      ENDIF
      IF (INPCI.EQ.0) THEN
        OPEN (UNIT = 68,FILE=TRIM(NAME(1))//'.cbm',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ELSE
        OPEN (UNIT = 68,FILE=TRIM(NAME(1))//'.bm',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ENDIF
      READ(68,IOSTAT=IOS) G92MIX
      IF ((IOS.NE.0).OR.(G92MIX.NE.'G92MIX')) THEN
        WRITE(*,*) 'Not a GRASP mixing file'
        STOP
      ENDIF
      WRITE(78) G92MIX
      READ (68) nelec, ncftot, nw, nvectot, nvecsiz, nblock
      WRITE(78) nelec, ncftot, nw, nvectot, nvecsiz, nblock

      DO IBLK =1,NBLOCK
      READ(68) ib,ncf,NVEC,iatjp,iaspa
      WRITE(78) ib,ncf,NVEC,iatjp,iaspa
      CALL ALLOC(PNEVAL,NVEC,8)
      CALL ALLOC(PNEVEC,NCF*NVEC,8)
      CALL ALLOC(PNIVEC,NVEC,4)
      READ(68) (IVEC(I),I=1,NVEC)
      READ(68) EAV,(EVAL(I),I=1,NVEC)
      READ(68) ((EVEC(I+(J-1)*NCF),I=1,NCF),J=1,NVEC)

      WRITE(78) (IVEC(I),I=1,NVEC)
      WRITE(78) EAV,(EVAL(I),I=1,NVEC)
      WRITE(78) ((EVEC(I+(J-1)*NCF),I=1,NCF),J=1,NVEC)
      CALL DALLOC(PNEVAL)
      CALL DALLOC(PNEVEC)
      CALL DALLOC(PNIVEC)
      ENDDO
      NAME(2)= NAME(1)
      CLOSE(68)
      CLOSE(78)
      RETURN
      END
