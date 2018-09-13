! Collection of FORTRAN subroutines for general purpose MPI operations.
! No COMMONs
! Mostly IMPLICIT NONE
! Only INCLUDE 'mpif.h'

!***********************************************************************
      subroutine startmpi (myid, nprocs, host, lenhost)
      implicit none
      integer myid, nprocs, lenhost
      character*(*) host

! This subroutine starts MPI to get the process id (myid), the number
! of processes (nprocs) and the status of enrollment (ierr)

      include 'mpif.h'
      integer ierr

      call MPI_Init (ierr)
      call MPI_Comm_rank (MPI_COMM_WORLD, myid, ierr)
      call MPI_Comm_size (MPI_COMM_WORLD, nprocs, ierr)
      call MPI_Get_processor_name (host, lenhost, ierr)
      if (ierr .ne. MPI_SUCCESS) print *, ' ierr=', ierr, ' id=', myid

      return
      end

!***********************************************************************
      subroutine startmpi2 (myid, nprocs, host, lenhost, ncount1,
     &                     startdir, permdir, tmpdir, progname)
      implicit none
      integer myid, nprocs, lenhost, ncount1
      character*(*) host, startdir, permdir, tmpdir, progname

! Calls startmpi to get mpi environmen;
! Calls cpath to get various paths;
! Calls DATE_AND_TIME to get date, time, zone;
! Calls mpi_printmsg to print some of these info to screen.

      character idstring*3

      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

      integer lenstart, lenperm, lentmp, ncount_rate, ncount_max
! For printing
		character msg*80

*=======================================================================
*  Get processor info: myid, nprocs, host name; and print
*=======================================================================

      CALL startmpi (myid, nprocs, host, lenhost)
      WRITE (idstring, '(I3.3)') myid
      IF (myid .EQ. 0) THEN
         print *, '===================================================='
         print *, '       ', progname, ': Execution Begins ...'
         print *, '===================================================='
         print *,        'Participating nodes:'
       ENDIF
      msg = '  Host: ' // host(1:lenhost) // '    ID: ' // idstring
      CALL mpix_printmsg (msg, myid, nprocs)

*=======================================================================
*  Get date, time, zone and print
*=======================================================================

      CALL DATE_AND_TIME (chdate, chtime, chzone, nYMDUHMSM)
      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Date and Time:'
      ENDIF
      msg = '  ' // host(1:lenhost) // ': ' //
     &      '  Date: ' // chdate //
     &      '  Time: ' // chtime //
     &      '  Zone: ' // chzone
      CALL mpix_printmsg (msg, myid, nprocs)

*=======================================================================
*  Set up local working dir and go there
*     tmpdir  - local working dir of the node. mcpXXX files are there
*     permdir - for I/O specific to node-0.
*=======================================================================

!     CALL cpath (startdir, permdir, tmpdir)

      lenstart = LEN_TRIM (startdir)
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Start Dir:'
      ENDIF
      msg = '  ' // host(1:lenhost) // ': ' // startdir(1:lenstart)
      CALL mpix_printmsg (msg, myid, nprocs)

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Serial I/O Dir (node-0 only):'
         PRINT *, '  ' // host(1:lenhost) // ': ' // permdir(1:lenperm)
      ENDIF

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *, 'Work Dir (Parallel I/O):'
      ENDIF
      msg = '  ' // host(1:lenhost) // ': ' // tmpdir(1:lentmp)
      CALL mpix_printmsg (msg, myid, nprocs)

*=======================================================================
*  Start timing - Record the wall clock
*=======================================================================

      CALL SYSTEM_CLOCK (ncount1, ncount_rate, ncount_max)

      return
      end
!***********************************************************************
      subroutine stopmpi (what, myid)
      implicit none
      character*(*) what
      integer myid

! This subroutine stops MPI, but before doing so it issues an error
! message to node-0 about what went wrong and where it happened.
!   what - string, type of error (most likely a subroutine name)
!   myid - The id of the PE where error happened

      include 'mpif.h'
      integer ierr

      print *, 'mpi stopped by node-', myid, ' from ', what
      call MPI_Barrier (MPI_COMM_WORLD,ierr)
      call MPI_Finalize (ierr)

      stop
      !return
      end

!***********************************************************************
      subroutine stopmpi2 (myid, nprocs, host, lenhost, ncount1, 
     &                      progname)
      implicit none
      integer myid, nprocs, lenhost, ncount1
      character*(*) host, progname

! Calls DATE_AND_TIME to get date, time, zone;
! Calls mpi_printmsg to print some of these info to screen.

! Things for timing
      INTEGER   ncount2, ncount_rate, ncount_max, nseconds
      character chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      integer  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

! For printing
		character str2nds*8, msg*80

      include 'mpif.h'
      integer ierr

*=======================================================================
*  Get processor info: myid, nprocs, host name; and print
*=======================================================================

      if (myid .eq. 0) then
         print *, '===================================================='
         print *, '       ', progname, ': Execution Finished ...'
         print *, '===================================================='
         print *,        'Wall time:'
       endif

      call system_clock (ncount2, ncount_rate, ncount_max)
		ncount2 = ncount2 - ncount1
		nseconds = ncount2 / ncount_rate
		write (str2nds, '(i8)') nseconds
		msg = str2nds // ' seconds on ' // host(1:lenhost)
		call mpix_printmsg (msg, myid, nprocs)

      if (myid .eq. 0) then
         print *
         print *, 'Finish Date and Time:'
      endif

      call date_and_time (chdate, chtime, chzone, nYMDUHMSM)

      msg = '  ' // host(1:lenhost) // ': ' //
     &      '  Date: ' // chdate //
     &      '  Time: ' // chtime //
     &      '  Zone: ' // chzone
      CALL mpix_printmsg (msg, myid, nprocs)

      if (myid .eq. 0) print *

      call MPI_Barrier (MPI_COMM_WORLD,ierr)
      call stopmpi (progname // ': Execution complete.', myid)

      return
      end
!***********************************************************************
      subroutine mpix_printmsg (msg, myid, nprocs)
      implicit none
      character*(*) msg
      integer myid, nprocs

! Displays on node-0's screen info from all nodes including node 0.

      include 'mpif.h'
      integer inID, istat(MPI_STATUS_SIZE), ierr, msgLength

      msgLength = len_trim (msg)
      
      if (myid .ne. 0) then
         call MPI_Send (msgLength, 1, MPI_INTEGER, 0, myid,
     &                  MPI_COMM_WORLD, ierr)   ! Send nsgLength
         call MPI_Send (msg, msgLength, MPI_CHARACTER, 0, myid+nprocs,
     &                  MPI_COMM_WORLD, ierr)   ! Send msg
      else
         print *, msg(1:msgLength)      ! msg from node 0 itself
         do inID = 1, nprocs - 1
            call MPI_Recv (msgLength, 1, MPI_INTEGER, inID,
     &                     inID, MPI_COMM_WORLD, istat, ierr)
            call MPI_Recv (msg, msgLength, MPI_CHARACTER, inID,
     &                     inID+nprocs, MPI_COMM_WORLD, istat, ierr)
            print *, msg(1:msgLength)
         enddo
      endif

      return
      end

!***********************************************************************
      subroutine mpix_chkpt (myid, what)
      implicit none
      integer myid
      character *(*) what

!  To set a check-point in mpi program

      include 'mpif.h'
      integer ierr
      character chtmp*1

      if (myid .eq. 0) then
         print *, what, 'hit any key to continue'
         read (*, '(a)') chtmp
      endif

      call mpi_barrier (MPI_COMM_WORLD, ierr)

      return
      end

!***********************************************************************
      subroutine mpix_bytes (n, newType, ierr)
      implicit none
      integer n, newType, ierr

!  Constructs new mpi data type newType of n-bytes long

      include 'mpif.h'
      integer ier0

      call MPI_Type_contiguous (n, MPI_BYTE, newType, ier0)
      call MPI_Type_commit (newType, ierr)

      ierr = ierr + ier0

      return
      end
!***********************************************************************
      SUBROUTINE gdsummpi_ori (x, n)
      IMPLICIT REAL*8           (A-H, O-Z)
	
      INCLUDE 'mpif.h'
      PARAMETER (ITUNE = 4096)
      DIMENSION x(*)

      DIMENSION y(ITUNE)
      POINTER (P,tmp1(*))
      INTEGER itmp	

*	In order to make life easier and also have a tuning possibility
*	let's split x a little bit
*	should not create as much swapping as before ? 
*
      DO I = 1, N, ITUNE
         p = LOC (x(I))
         itmp = MIN (ITUNE, N-I+1)
         CALL MPI_Allreduce (tmp1, y, itmp, MPI_DOUBLE_PRECISION,
     $                             MPI_SUM, MPI_COMM_WORLD, ierr)
         CALL dcopy (itmp, y, 1, tmp1, 1)
      ENDDO

     	RETURN
     	END
!***********************************************************************
! From misha (gdsummpii.F), slightly modified
      SUBROUTINE gisummpi_ori (ix, n)
      IMPLICIT REAL*8           (A-H, O-Z)
	
      INCLUDE 'mpif.h'
      PARAMETER (ITUNE = 4096)
      DIMENSION ix(*)
      INTEGER iy(ITUNE)
      POINTER (P, jtmp(*))	
      INTEGER jtmp
	
*	In order to make life easier and also have a tuning possibility
*	let's split x a little bit
*	should not create as much swapping as before ? 

      DO I =1, N, ITUNE
         p = LOC (ix(I))
         itmp = MIN (ITUNE, N-I+1)
         CALL MPI_Allreduce (jtmp, iy, itmp, MPI_INTEGER, MPI_SUM,
     $                        MPI_COMM_WORLD, ierr)
         CALL icopy (itmp, iy, 1, jtmp, 1)
      ENDDO

     	RETURN
     	END
!***********************************************************************
      SUBROUTINE gdsummpi (x, n)
      IMPLICIT NONE
	
      INCLUDE 'mpif.h'
      INTEGER ierr, n
      REAL*8           x(n), y(n)

      CALL dinit (n, 0.d0, y, 1)
      CALL MPI_Allreduce (x, y, n, MPI_DOUBLE_PRECISION,
     &                             MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL dcopy (n, y, 1, x, 1)   ! copy y to x

      RETURN
      END

!***********************************************************************
      SUBROUTINE gisummpi (ix, n)
      IMPLICIT NONE
	
      INCLUDE 'mpif.h'
      INTEGER ix(n), iy(n), n, ierr, i
	
      DO i = 1, n
         iy(i) = 0
      ENDDO

      CALL MPI_Allreduce (ix, iy, n, MPI_INTEGER, MPI_SUM,
     &                               MPI_COMM_WORLD, ierr)
      CALL icopy (n, iy, 1, ix, 1)

      RETURN
      END
