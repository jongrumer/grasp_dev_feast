!***********************************************************************
      subroutine sys_chdir (dir, lendir, ierr, machine)
      implicit none
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr
      character(len=*), parameter::iam = 'PVP'  ! soly for machine
      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine changes current working directory to dir.
!   lendir is the length of character string dir; 
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-21
!
!***********************************************************************

      integer chdir

      machine = iam

      ierr = chdir (dir(1:lendir))

      return
      end
