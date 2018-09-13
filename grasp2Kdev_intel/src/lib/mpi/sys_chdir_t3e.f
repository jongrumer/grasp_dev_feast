!***********************************************************************
c      subroutine sys_chdir (dir, lendir, ierr, machine)
      subroutine sys_chdir (dir, lendir, ierr)
      implicit none
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr
c      character(len=*), parameter::iam = 'T3E'  ! soly for machine
c      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine changes current working directory to dir.
!   lendir is the length of character string dir; 
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!***********************************************************************

      call pxfchdir (dir, lendir, ierr)  ! T3E
      if (ierr .ne. 0) then
         print *, dir(1:lendir), ' doesn''t exist: will be created'
      endif

      return
      end
