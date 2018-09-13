************************************************************************
*                                                                      *
      SUBROUTINE CONVRT2 (INTNUM,CNUM,LENTH, from)
*                                                                      *
*   Converts the  INTEGER number  INTNUM  into the  CHARACTER string   *
*   CNUM of length LENTH. INTEGER lengths of up to 64 bits are acco-   *
*   modated.                                                           *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 22 Sep 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*(*) CNUM, from
      CHARACTER*6 FORM
      CHARACTER*2 C1020
      CHARACTER*1 C19
*
      DIMENSION C19(1:9)
      DIMENSION C1020(0:10)
*
      DATA C19   /      '1', '2', '3', '4', '5', '6', '7', '8', '9'/
      DATA C1020 /'10','11','12','13','14','15','16','17','18','19',
     :            '20'/
      common/iounit/istdi, istdo, istde
*
      IF     (INTNUM .LT. 0) THEN
         LENTH = LOG10 (DBLE (-INTNUM))+2
      ELSEIF (INTNUM .EQ. 0) THEN
         LENTH = 1
      ELSE
         LENTH = LOG10 (DBLE ( INTNUM))+1
      ENDIF
*
*   Ensure that the length of CNUM as dimensioned is adequate;
*   stop with an error message if it isn't
*
      IF (LENTH .GT. LEN (CNUM)) THEN
         write(istde,*) 'CONVRT: Length of CNUM inadeuate. (from:'
     &                 , from, ')'
         STOP
      ELSE
         IF (LENTH .LE. 9) THEN
            FORM = '(1I'//C19(LENTH)//')'
            WRITE (CNUM(1:LENTH),FORM(1:5)) INTNUM
         ELSE
            FORM = '(1I'//C1020(LENTH-10)//')'
            WRITE (CNUM(1:LENTH),FORM(1:6)) INTNUM
         ENDIF
      ENDIF
*
      RETURN
      END
