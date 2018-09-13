!***********************************************************************
      subroutine stoptime (ncount1, progname)
      implicit none
      integer ncount1
      character*(*) progname

! Calls DATE_AND_TIME to get date, time, zone;

! Things for timing
      INTEGER   ncount2, ncount_rate, ncount_max, nseconds
      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

! For printing
		character str2nds*8, msg*80

*=======================================================================
*  Get processor info: myid, nprocs, host name; and print
*=======================================================================

C      print *, '       ', progname, ': Execution Finished ...'
      print *
      print *,        'Wall time:'

      call system_clock (ncount2, ncount_rate, ncount_max)
		ncount2 = ncount2 - ncount1
		nseconds = ncount2 / ncount_rate
		write (str2nds, '(i8)') nseconds
		msg = str2nds // ' seconds'
		print *, msg(1:len_trim(msg))

      print *
      print *, 'Finish Date and Time:'

      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

CGG      msg = '  Date: ' // chdate //
CGG     &      '  Time: ' // chtime //
CGG     &      '  Zone: ' // chzone
CGG      print *, msg(1:len_trim(msg))
      Print*, '  Date (Yr/Mon/Day): ',
     :        chdate(1:4),'/',chdate(5:6),'/',chdate(7:8)
      Print*, '  Time (Hr/Min/Sec): ',
     :        chtime(1:2),'/',chtime(3:4),'/',chtime(5:10)
      Print*, '  Zone: ',chzone

      print *
      print *, progname // ': Execution complete.'

      return
      end
