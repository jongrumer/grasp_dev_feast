!***********************************************************************
      subroutine sys_getwd (dir, ierr, machine)
      implicit none
      character(len=*), intent(out):: dir
      integer, intent(out):: ierr
      character(len=*), parameter::iam = 'DIC'  ! soly for machine
      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine gets current working directory and assigns it to dir.
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-24
!
!***********************************************************************

      machine = iam

      dir = 'N/A'
      ierr = 0

      return
      end
