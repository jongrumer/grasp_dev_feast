************************************************************************
*                                                                      *
      SUBROUTINE GETMIX (NAME,INPCI,IBLK)
*                                                                      *
*   Open, check, load data from and close the rscf.mix file.           *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, LENGTH, OPENFL.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 25 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*24 NAME
      CHARACTER*6 G92MIX
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
      common/BLKTOT/NELECTOT,NCFTOT,NWTOT,NVECTOT,NVECSIZTOT,NBLOCK1
*

      IK = 30
      IF(IBLK.EQ.1) THEN
        J = INDEX(NAME,' ')                                            
        IF (INPCI.EQ.0) THEN
          OPEN (UNIT = IK,FILE=NAME(1:J-1)//'.cm',FORM='UNFORMATTED',
     :         STATUS='OLD')
        ELSE
          OPEN (UNIT = IK,FILE=NAME(1:J-1)//'.m',FORM='UNFORMATTED',
     :         STATUS='OLD')
        ENDIF

        READ (IK,IOSTAT = IOS) G92MIX
        IF ((IOS .NE. 0) .OR.
     :      (G92MIX .NE. 'G92MIX')) THEN
           PRINT *, 'File',IK,'Not a GRASP MIXing Coefficients File;'
           CLOSE (IK)
           STOP 
        ENDIF
*
        READ (IK) NELECTOT,NCFTOT,NWTOT,NVECTOT,NVECSIZTOT,NBLOCK1
        IF ((NELEC .NE. NELECTOT) .OR.
     :      (NW .NE. NWTOT)) THEN
C     :    (NCF .NE. NCFT) .OR.
           PRINT *, 'File',IK,'is not an'
           PRINT *, ' appropriate to Coefficients file'
           CLOSE (IK)
           STOP
        ENDIF
      ENDIF
*
*   Load data from the  rscf.mix  file
*
      PRINT *, 'Loading MIXing Coefficients File ...'
*
      READ (IK) NB,NCF,NVEC,IATJP,IASPA
      CALL ALLOC (PNEVAL,NVEC,8)
      CALL ALLOC (PNEVEC,NCF*NVEC,8)
      CALL ALLOC (PNIVEC,NVEC,4)
      CALL ALLOC (PIATJP,NVEC,4)
      CALL ALLOC (PIASPA,NVEC,4)
*
*   These arrays are deallocated in mcp
*
      READ (IK) (IVEC(I),I = 1,NVEC)
c     READ (IK) (IATJPO(I),IASPAR(I),I = 1,NVEC)
      DO I = 1,NVEC
      IATJPO(I)=IATJP
      IASPAR(I)=IASPA
      ENDDO
      READ (IK) EAV,(EVAL(I),I = 1,NVEC)
      READ (IK) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
*
      PRINT *, ' ... load complete;'
*
*   Close the  rscf.mix  file
*
c     CLOSE (IK)
*
      RETURN
      END
