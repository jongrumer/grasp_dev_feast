************************************************************************
*                                                                      *
      SUBROUTINE SETTMP (nb, kmax, filehead)
      IMPLICIT NONE
*
*   Opens temporary files tmp.XX in units XX. XX=30 to 32+kmax
*
*   File shared (hard link) by mcpmpi, mcpblk
*
*   Modified by Xinghong He               Last revision:  1 Jul 1998   *
*                                                                      *
************************************************************************
      INTEGER   nb, kmax
      CHARACTER filehead*(*)

      INTEGER   istdi, istdo, istde, k,lck,ierr,i, lng
      CHARACTER CK*2

      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
*
* All files  filehead.XX  are UNFORMATTED;
*
      lng = LEN_TRIM (filehead)
      DO K = 30, 32 + kmax
         CALL CONVRT (K, CK, LCK)
         CALL OPENFL (K, filehead(1:lng)//'.'//CK(1:2), 
     &               'UNFORMATTED', 'UNKNOWN', IERR)
         IF (IERR .NE. 0) THEN
            DO I = 30, K
               CLOSE (I)
            ENDDO
            WRITE (istde,*) 'Error when opening the tmp files'
            STOP
          ENDIF
      ENDDO

      DO K = 30, 32 + kmax
         WRITE (K) 'MCP', nb
      ENDDO

      RETURN
      END
