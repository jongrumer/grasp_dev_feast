************************************************************************
*                                                                      *
      SUBROUTINE ANGDATA(NAME,AVAIL,KAMAX)
*                                                                      *
*   Checks if the angular file name.T is available and appropriate     *
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
      COMMON/ORB2/NCF,NW,PNTRIQ

      J = INDEX(NAME,' ')
      INQUIRE (FILE = NAME(1:J-1)//'.TB', EXIST = FOUND)
      IF (.NOT.FOUND) THEN
        PRINT *, ' Angular file not available. Angular data computed'
        AVAIL = .FALSE.
        RETURN
      ELSE
*
*  Open the file and check if it is appropriate for the present case
*
        NF = 200
        OPEN (UNIT = NF,FILE = NAME(1:J-1)//'.TB',STATUS='OLD',
     :        FORM='UNFORMATTED')
        REWIND(NF)
        READ(NF) NCFD,NWD,KAMAXD
        IF (.NOT.(NCFD.EQ.NCF.AND.NWD.EQ.NW.AND.KAMAXD.EQ.KAMAX)) THEN
          PRINT *, ' Angular file not appropriate'
          AVAIL = .FALSE.
          RETURN
        ELSE
          PRINT *, ' Angular data read from file'
          AVAIL = .TRUE.
        ENDIF
      ENDIF
      RETURN
      END
*
*

