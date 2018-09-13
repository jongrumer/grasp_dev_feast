************************************************************************
*                                                                      *
      SUBROUTINE LODMIX (idblk)
*                                                                      *
*   Determines the eigenpairs required;  this information is written   *
*   to the head of the  .mix  file.                                    *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, lodstate                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Nov 1992   *
*   Block version by Xinghong He          Last revision:  9 Jun 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER idblk(*)*8
*
      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ

      POINTER (PCCMIN,ICCMIN(*))
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX   ! NCMAX not used throughout

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (pidxblk, idxblk(*))
      COMMON/blkidx/pidxblk

      COMMON/iounit/istdi,istdo,istde

!-----------------------------------------------------------------------

! lodstate generates
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
! Thus we let node-0 do it and then broadcast here

      CALL alloc (pncmaxblk, nblock, 4)
      CALL alloc (pnevblk, nblock, 4)

      CALL lodstate (nblock, ncfblk(1), idblk, nevblk, ncmaxblk)
*
*   Determine other auxiliary arrays, parameters
*
*    idxblk() is the block number of an eigenstate
*    nvecsiz  is the total size of the eigenvector array
*
      CALL alloc (pidxblk, ncmin, 4)
      ncftot = 0
      noffset = 0
      nvecsiz = 0
      DO jb = 1, nblock
         DO j = 1, nevblk(jb)
            idxblk(j + noffset) = jb
         ENDDO
         ncftot = ncftot + ncfblk(jb)
         noffset = noffset + nevblk(jb)
         nvecsiz = nvecsiz + ncfblk(jb) * nevblk(jb)
      ENDDO
      IF (noffset .NE. ncmin) STOP 'lodmix: ncmin trouble'
*
*   Write header data to the  .mix  file
*
      WRITE (25) NELEC, ncftot, NW, ncmin, nvecsiz, nblock
*
      RETURN
      END
