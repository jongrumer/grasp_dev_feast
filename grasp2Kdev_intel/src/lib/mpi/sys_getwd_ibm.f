!***********************************************************************
c      subroutine sys_getwd (dir, ierr, machine)
      subroutine sys_getwd (dir, ierr)
      implicit none
      character(len=*), intent(out):: dir
      integer, intent(out):: ierr
      character(len=*), parameter::iam = 'IBM'  ! soly for machine
c      character(len=len (iam)), optional, intent(out):: machine

!   This routine gets current working directory and assigns it to dir.
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-24
!
!***********************************************************************

      integer getwd

      integer myid, nprocs, ir
      COMMON /mpi/ myid, nprocs, ir 
c      CALL mpix_printmsg ('within sys_getwd_IBM ', myid, nprocs)
c      machine = iam

      ierr = getwd (dir)
      if (ierr .gt. 0) then   ! successful
         ierr = 0
      else if (ierr .eq. 0) then
         ierr = 1
      else
         print *, 'Can an address be negative ?'
         ierr = 2
      endif

c      CALL mpix_printmsg (' finishing sys_getwd_IBM ', myid, nprocs)
      return
      end
