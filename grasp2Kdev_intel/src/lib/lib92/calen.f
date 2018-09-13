************************************************************************
*                                                                      *
      SUBROUTINE CALEN (JTIME,JDATE)
*                                                                      *
*   Loads the character strings JTIME and JDATE with the time of day   *
*   and the date when called.                                          *
*                                                                      *
*                                         Last revision: 25 Sep 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*8 JTIME,JDATE
Cww      CHARACTER*8 clock_,date
*
Cww      JTIME = clock_ ()
Cww      JDATE = date ()
*
      RETURN
      END
