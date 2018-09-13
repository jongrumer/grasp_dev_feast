************************************************************************

      SUBROUTINE setsum (name)

*   Open the  .csum  file on stream 24.
*   Xinghong He                                          10 Jun 1998

************************************************************************

      CHARACTER*(*) name

!-----------------------------------------------------------------------

      k = INDEX (name,' ')

      CALL openfl (24,name(1:k-1)//'.csum','FORMATTED','UNKNOWN',IERR)

      RETURN
      END
