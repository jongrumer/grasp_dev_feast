************************************************************************
*                                                                      *
      SUBROUTINE ANGDATA(NAME,AVAIL)
*                                                                      *
*   Checks if the angular file name.TH is available and appropriate    *
*                                                                      *
*   Written by Per Jonsson                      6 March 1997           *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL AVAIL,FOUND
Cww      INTEGER PNTRIQ                                         
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*24 NAME
*
      COMMON/PRNT/NVEC,PNIVEC,NVECMX
      COMMON/ORB2/NCF,NW,PNTRIQ

      NF = 129
      J = INDEX(NAME,' ')
      INQUIRE (FILE = NAME(1:J-1)//'.TH', EXIST = FOUND)
      IF (.NOT.FOUND) THEN
        PRINT *, ' Angular file not available'
        AVAIL = .FALSE.
        OPEN (UNIT = NF,FILE = NAME(1:J-1)//'.TH',STATUS='NEW',
     :        FORM='UNFORMATTED')
        WRITE(NF) NCF,NW,NVEC
      ELSE
*
*  Open the file and check if it is appropriate for the present case
*
        OPEN (UNIT = NF,FILE = NAME(1:J-1)//'.TH',STATUS='OLD',
     :        FORM='UNFORMATTED')
        REWIND(NF)
        READ(NF) NCFD,NWD,NVECD
        IF (.NOT.(NCFD.EQ.NCF.AND.NWD.EQ.NW.AND.NVECD.EQ.NVEC)) THEN
          PRINT *, ' Angular file not appropriate.'
          PRINT *, ' Remove ',NAME(1:J-1)//'.TH',' and rerun program.'
          STOP 
        ELSE
          PRINT *, ' Angular data read from file'
          AVAIL = .TRUE.
        ENDIF
      ENDIF

      RETURN
      END
*
*

