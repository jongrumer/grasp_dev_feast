************************************************************************
*                                                                      *
      SUBROUTINE ANGDATA(NAME,AVAIL,JKP,NFILE2)
*                                                                      *
*   Checks if the angular file name(1).name(2).T is available          *
*   and appropriate                                                    *
*                                                                      *
*   Written by Per Jonsson                      6 March 1997           *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL AVAIL,FOUND
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)

      POINTER (PNTRKP,KP(1))

      CHARACTER*24 NAME(2)
      CHARACTER*2 S(-9:9)
*
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /OSC6/NKP,PNTRKP

      S(-9)  = '-9'
      S(-8)  = '-8'                                                     
      S(-7)  = '-7'                                                     
      S(-6)  = '-6'                                                     
      S(-5)  = '-5'                                                     
      S(-4)  = '-4'                                                     
      S(-3)  = '-3'                                                     
      S(-2)  = '-2'                                                     
      S(-1)  = '-1'                                                     
      S(0)   = '+0'                                                     
      S(1)   = '+1'   
      S(2)   = '+2'                
      S(3)   = '+3'       
      S(4)   = '+4'   
      S(5)   = '+5'       
      S(6)   = '+6'    
      S(7)   = '+7'
      S(8)   = '+8'
      S(9)   = '+9'

      J1 = INDEX(NAME(1),' ')
      J2 = INDEX(NAME(2),' ')
      INQUIRE (FILE =
     :         NAME(1)(1:J1-1)//'.'//NAME(2)(1:J2-1)//'.'
     :         //S(KP(JKP))//'T',
     :         EXIST = FOUND)
      IF (.NOT.FOUND) THEN
        PRINT *
        PRINT *, ' Angular file not available'
        AVAIL = .FALSE.
        RETURN
      ELSE
*
*  Open the file and check if it is appropriate for the present case
*
        OPEN (UNIT = NFILE2,FILE = 
     :        NAME(1)(1:J1-1)//'.'//NAME(2)(1:J2-1)//'.'
     :        //S(KP(JKP))//'T',
     :        STATUS='OLD',FORM='UNFORMATTED')
        REWIND(NFILE2)
        READ(NFILE2) iblki,iblkf,NWD,NKPD
        IF (.NOT.(NWD.EQ.NW.AND.NKPD.EQ.NKP)) THEN
          PRINT *, ' Angular file not appropriate'
          AVAIL = .FALSE.
          CLOSE (NFILE2,STATUS='DELETE')
          RETURN
        ELSE
          rewind (nfile2)
          PRINT *, ' Angular data read from file'
          AVAIL = .TRUE.
        ENDIF
      ENDIF
      RETURN
      END
*
*

