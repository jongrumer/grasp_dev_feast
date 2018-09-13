!***********************************************************************
      subroutine cpath (startdir, permdir, tmpdir)
      implicit none
      character(len=*), intent(out):: startdir, permdir, tmpdir

!   This subroutine returns path names of local data disks.
!
!      startdir - path where the current node started from.
!      permdir  - path where node-0 performs serial i/o.
!      tmpdir   - path where the current node goes to.
!
!   It is developed by modifying Misha's code create_paths
!   where some C/system functions are called to get the current working
!   directory (permdir) and to create tmpdir if it does not exist.
!
!   This version reads (by node-0) the paths from a disk file under the 
!   starting directory of the node-0, determine the length and do 
!   sending/receiving. Only if the paths defined here do not exist will
!   C functions be called to create them.
!
!   Xinghong He  98-10-30
!
!***********************************************************************
      !...mpi
      include 'mpif.h'
      integer  myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      COMMON /mpi/ myid, nprocs, ierr
      !...locals - be careful when change lenidstring !
      integer, parameter:: lendisk0 = 128, lenidstring = 3
      character(len=lendisk0) disk
      character(len=lenidstring) idstring
      integer   lendisk, i
!=======================================================================
!  Open file, read paths and send/receive them. Each node will have 
!  its preliminary path stored in variable disk. In addition, node-0 
!  will have the current working dir stored in permdir.
!=======================================================================
      if (myid .eq. 0) then
         open (unit=1001, file='disks', status='old')
         !...paths for serial i/o, node-0 only
         read (1001,*) permdir  ! paths for serial i/o, node-0 only
         read (1001,*) tmpdir   ! tempory for local disk of node-0 
         !...paths for slaves, read and send; 
         do i = 1, nprocs - 1
            read (1001,*) disk
            call MPI_Send (disk, lendisk0, MPI_CHARACTER, i, i, 
     &                                MPI_COMM_WORLD, ierr)
         enddo
         disk = tmpdir           ! local disk of node-0 
         close (1001)
      else 
         !...slaves receive their local dirs
         if (nprocs .gt. 1)
     &   call MPI_Recv (disk, lendisk0, MPI_CHARACTER, 0, myid,
     &                             MPI_COMM_WORLD, istat, ierr)
      endif
      lendisk = len_trim (disk)
!=======================================================================
!  Get the current dir name. This part can be removed - but why remove
!=======================================================================
      call sys_getwd (startdir, ierr)
      if (ierr .ne. 0) then
         print *, 'Failed to get current working dir, myid = ', myid
         ! goto 999 ! Not a serious error, go on
      endif
!=======================================================================
!  Go to local disk - They must have been there.
!=======================================================================
      call sys_chdir (disk, lendisk, ierr)
      if (ierr .ne. 0) then
         print *, 'Failed to go to pre-defined dir ' // disk(1:lendisk)
         goto 999
      endif
!=======================================================================
!  Go to sub-dir defined by its identification number. Create it if
!  not there.
!=======================================================================
      write (idstring,'(I3.3)') myid
      call sys_chdir (idstring, lenidstring, ierr)
      if (ierr .ne. 0) then
         call sys_mkdir (idstring, lenidstring, ierr)
         if (ierr .ne. 0) then
            print *, 'Failed to make sub-dir ' // idstring
            goto 999
         endif
         call sys_chdir (idstring, lenidstring, ierr)  ! try again
         if (ierr .ne. 0) then
            print *, 'Failed to go to sub-dir ' // idstring
            goto 999
         endif
      endif
  999 continue
!=======================================================================
!  Handle error cases. Since it has to succeed, we invoke stop.
!=======================================================================
      if (ierr .ne. 0) then
         print *, 'cpath failed, myid = ', myid
         stop
      endif
      tmpdir = disk
      return
      end
