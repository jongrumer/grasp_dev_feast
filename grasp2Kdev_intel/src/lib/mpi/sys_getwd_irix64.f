!***********************************************************************
c      subroutine sys_getwd (dir, ierr, machine)
      subroutine sys_getwd (dir, ierr)
      implicit none
      character(len=*), intent(out):: dir
      integer, intent(out):: ierr
      character(len=*), parameter::iam = 'IRIX64'  ! soly for machine
c      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine gets current working directory and assigns it to dir.
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-24
!
!***********************************************************************

c      machine = iam

      dir = 'N/A'
      ierr = 0

      return
      end
