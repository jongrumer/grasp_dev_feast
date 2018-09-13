************************************************************************
*                                                                      *
      SUBROUTINE LDLBL2 (NAME)
*                                                                      *
*   Open, check and load data from the  .lsj.lbl   file of the         *
*   final state.                                                       *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC.                                        *
*                                                                      *
*   Written by G. Gaigalas,                                            *
*   NIST                                                  May 2011     *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*1 Lev_par_2
      CHARACTER*4 Lev_J_2
      CHARACTER*24 NAME
      CHARACTER*64 string_CSF2
      POINTER(PNELev_POS_2,Lev_POS_2(1))
      POINTER(PNELev_J_2,Lev_J_2(1))
      POINTER(PNELev_Par_2,Lev_Par_2(1))
      POINTER(PNELev_ENER_2,RLev_ENER_2(1))
      POINTER(PNEstring_CSF2,string_CSF2(1))
      COMMON/JJ2LSJ2/ NVECTOTF,IOPEN_STATUS2,PNELev_POS_2,PNELev_J_2,
     :       PNELev_Par_2,PNELev_ENER_2,PNEstring_CSF2
*
*   Check the first record of the file; if not as expected, try again
*
      J = INDEX(NAME,' ')
      OPEN (UNIT = 31,FILE=NAME(1:J-1)//'.lsj.lbl',FORM='FORMATTED',
     :     STATUS='OLD',IOSTAT = IOPEN_STATUS2)
      IF(IOPEN_STATUS2 .EQ. 0) THEN
         READ (31,'(1A15)',IOSTAT = IOS) RECORD
         IF (IOS .NE. 0) THEN
            PRINT *, 'Not a i *.lsj.lbl  File;'
            CLOSE (31)
            STOP
         ELSE
            CALL ALLOC (PNELev_POS_2,NVECTOTF,4)
            CALL ALLOC (PNELev_J_2,NVECTOTF,4)
            CALL ALLOC (PNELev_Par_2,NVECTOTF,4)
            CALL ALLOC (PNELev_ENER_2,NVECTOTF,8)
            CALL ALLOC (PNEstring_CSF2,NVECTOTF,64)
*
            ICount = 1
            READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)
     :        Lev_Pos_2(ICount),Lev_J_2(ICount),Lev_Par_2(ICount),
     :        RLev_ENER_2(ICount)
            IF (IOS .NE. 0) GO TO 1
*
            READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF2(ICount)
*
    2       READ (31,'(1X,I2)',IOSTAT = IOS),ITEST
            IF (IOS .NE. 0) GO TO 1
            IF (ITEST .EQ. 0) GO TO 2
            BACKSPACE 31
            ICount = ICount + 1
            READ (31,'(1X,I2,1X,A4,5X,A1,8X,F16.9)',IOSTAT = IOS)
     :        Lev_Pos_2(ICount),Lev_J_2(ICount),Lev_Par_2(ICount),
     :        RLev_ENER_2(ICount)
            IF (IOS .NE. 0) GO TO 1
            READ (31,'(7X,F12.8,17X,A)') WEIGHTS,string_CSF2(ICount)
            GO TO 2
         ENDIF
      ENDIF
    1 CONTINUE
      CLOSE (31)
      RETURN
      END
