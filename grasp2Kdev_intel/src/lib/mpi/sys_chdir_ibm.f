!***********************************************************************
c      subroutine sys_chdir (dir, lendir, ierr, machine)
      subroutine sys_chdir (dir, lendir, ierr)
      implicit none
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr
      character(len=*), parameter::iam = 'IBM'  ! soly for machine
c      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine changes current working directory to dir.
!   lendir is the length of character string dir; 
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-21
!
!***********************************************************************

      integer chdir

c      integer myid, nprocs, ir
c      COMMON /mpi/ myid, nprocs, ir 
c      CALL mpix_printmsg ('entering sys_chdir_IBM ', myid, nprocs)
c      machine = iam
c      CALL mpix_printmsg ('       dir = ' // dir, myid, nprocs)

      ierr = chdir (dir(1:lendir) // '\0')

      return
      end
