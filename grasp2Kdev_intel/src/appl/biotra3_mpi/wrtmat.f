************************************************************************
*                                                                      *
      SUBROUTINE WRTMAT(A,NROW,NCOL,NMROW,NMCOL)
*                                                                      *
************************************************************************
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NMROW,NMCOL)
 
      DO I=1,NROW
        WRITE(6,1010) I,(A(I,J),J=1,NCOL)
      ENDDO

 1010 FORMAT(1H0,I5,2X,4(1X,E14.8),/,(1H ,7X,4(1X,E14.8)))

      RETURN
      END
