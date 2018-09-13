************************************************************************
*                                                                      *
      SUBROUTINE genmat (atwinv, jblock, myid, nprocs, elsto, irestart,
     : slf_en)
      IMPLICIT REAL*8          (A-H, O-Z)
*
*        Generate Hamiltonian matrix for all blocks
*        This routine calls setham to do the computation. It makes 
*        sure that the hamiltonian matrix is complete.
*  Only node-0 has the correct, non-zero elsto. And in any case, elsto
*  is not added to the matrix elements in files *.res
*  This routine has been re-written to work in both single- and
*  multi- processors. 
*
* Xinghong He 1998-06-23
*
************************************************************************
      DIMENSION SLF_EN(*)
      POINTER (pntriq,dummy)
      POINTER (PNEVAL,EVALdum)
      COMMON/EIGVAL/EAV,PNEVAL
     :      /FOPARM/ICCUT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /WHERE/IMCDF

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT

      COMMON/setham_to_genmat2/CUTOFFtmp,
     &  NCOEItmp, NCOECtmp, NCTEItmp, NCTECtmp, NTPItmp(6), NMCBPtmp, 
     &  NCOREtmp, NVPItmp, NKEItmp, NVINTItmp, NELMNTtmp, NCFtmp

      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
      nelmnt = 0     ! Counting continues in setham
      eav    = 0.d0
      elsto  = 0.d0

! See how much had been done (Hamiltonian matrix)
! irestart is set;
! iread accumulated;
! nelmnt, eav, elsto obtained (to be further modified in setham)

      iread = 0      ! # of rows read, initialization necessary

      READ (imcdf, IOSTAT=ios) ncfdum, ICCUTdum, myiddum, nprocsdum
      irestart = ios

      IF (ios .EQ. 0) THEN

         IF (ncf .NE. ncfdum .OR. iccut .NE. iccutdum .OR.
     &       myid .NE. myiddum .OR. nprocs .NE. nprocsdum) THEN
            WRITE (istde,*) ncf, ncfdum, iccut, iccutdum, 
     &                      myid, myiddum, nprocs, nprocsdum, 'check'
            STOP 'genmat:1'
         ENDIF

         DO i = myid + 1, ncf, nprocs
            READ (imcdf, IOSTAT=ios2) NELC,STOEL,(dum,IR=2,NELC),eav0,
     &                              (IROWdum,IR = 1,NELC)
                            ! Lower triangle row-mode, diagonal last
            IF (ios2 .EQ. 0) THEN
               iread = iread + 1
               nelmnt = nelmnt + nelc
               eav    = eav + eav0
               elsto  = stoel
            ELSE
               EXIT
            ENDIF
         ENDDO
         ipos = 7+nw+nw + iread + 1
      ELSE
         ipos = 7+nw+nw
      ENDIF
      
! Find the maximum number of rows

      nrows = (ncf - myid - 1 + nprocs) / nprocs
      IF (ncf .LT. nprocs) nrows = ncf / (myid+1)

! Report the number of rows read. 
! A more suitable report on all nodes can be done here, but this will
! set a synchronization point.

      WRITE (istde,*) iread,' (total ',nrows,') rows read from .res'
      IF (myid .EQ. 0) THEN
         WRITE (24,*) iread,' (total ',nrows,') rows read from .res'
      ENDIF

! Position the file for the next record from setham

      DO i = 1, jblock - 1
         j = (ncfblk(i) - myid - 1 + nprocs) / nprocs
         IF (ncfblk(i) .LT. nprocs) j = ncfblk(i) / (myid+1)
         ipos = ipos + j + 1
      ENDDO
      CALL posfile (0, imcdf, ipos)

      IF (ios .NE. 0) THEN
         WRITE (imcdf) ncf, iccut, myid, nprocs
         ! ncf, iccut are specific to the current block
      ENDIF

      IF (iread .LT. nrows) THEN
         icstrt = iread * nprocs + myid + 1
*     ...Generate the rest of the Hamiltonian matrix
         CALL setham (myid, nprocs, jblock, elsto, icstrt, nelmnt
     &                , atwinv,slf_en)
      ELSE
         NELMNTtmp = NELMNT
         NCFtmp = NCF
      ENDIF

      RETURN
      END
