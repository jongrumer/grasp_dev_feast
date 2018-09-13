!***********************************************************************
      subroutine starttime (ncount1, progname)
      implicit none
      integer ncount1
      character*(*) progname

! Calls DATE_AND_TIME to get date, time, zone;


      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

      integer ncount_rate, ncount_max

! For printing
		character msg*80

C      print *, '===================================================='
C      print *, '       ', progname, ': Execution Begins ...'
C      print *, '===================================================='

*=======================================================================
*  Get date, time, zone and print
*=======================================================================

CGG      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)
CGG      print *, 'Date and Time:'
CGG      msg = '  Date: ' // chdate //
CGG     &      '  Time: ' // chtime //
CGG     &      '  Zone: ' // chzone
CGG      print *, msg(1:len_trim(msg))
CGG      Print*, '  Date (Yr/Mon/Day): ',
CGG     :        chdate(1:4),'/',chdate(5:6),'/',chdate(7:8)
CGG      Print*, '  Time (Hr/Min/Sec): ',
CGG     :        chtime(1:2),'/',chtime(3:4),'/',chtime(5:10)
CGG      Print*, '  Zone: ',chzone

*=======================================================================
*  Start timing - Record the wall clock
*=======================================================================

      Call system_clock (ncount1, ncount_rate, ncount_max)

      return
      end
