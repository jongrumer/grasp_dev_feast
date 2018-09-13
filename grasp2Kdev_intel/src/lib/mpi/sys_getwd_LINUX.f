!***********************************************************************
c      subroutine sys_getwd (dir, ierr, machine)
      subroutine sys_getwd (dir,lcwd)
      implicit none
      character(len=*), intent(out):: dir
      integer, intent(out):: lcwd
      integer serr
      character(len=*), parameter::iam = 'SUN'  ! soly for machine
c      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine gets current working directory and assigns it to dir.
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-24
!
!***********************************************************************

      integer getcwd
c      machine = iam
      serr = getcwd(dir)
      lcwd = len_trim(dir)

      if (serr.ne.0) then
         print*, 'couldn''t get the current directory, exiting...'
         call exit(23);
      end if


      return
      end
