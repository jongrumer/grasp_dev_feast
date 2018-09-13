************************************************************************
      SUBROUTINE SETMIX (name)
      IMPLICIT REAL*8          (A-H,O-Z)
      CHARACTER name*(*)

*   Opens the  .mix  file on stream 25; writes a header to this file.  *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*   Modified by Xinghong He               Last revision: 13 Jul 1998   *
*
************************************************************************

      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNTRIQ,RIQDUM)

      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /ORB2/NCF,NW,PNTRIQ

      !...Wants nblock only
      POINTER (pncfblk,ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
      CALL OPENFL (25, name,'UNFORMATTED', 'NEW', IERR)
      IF (IERR .NE. 0) THEN
         WRITE (istde,*) 'Error when opening ', name(1:LEN_TRIM (name))
         STOP
      ENDIF
*
*   Write the file header
*
      WRITE (25) 'G92MIX'
      WRITE (25) NELEC, NCF, NW, 0, 0, nblock
*     ...The above record will be overidden in matrix.f
*        with the final form of
*     WRITE (25)  NELEC, NCF, NW, nvectot, nvecsiz, nblock

      RETURN
      END
