************************************************************************
*                                                                      *
      SUBROUTINE HMOUT
*                                                                      *
*   Routine for printing the Hamiltonian matrix.                       *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 21 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
*
      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
*
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
*
*
      ICI = 0
      DO 1 I = 1,NELMNT
         IRI = IROW(I)
         IF (I .GT. IENDC(ICI)) ICI = ICI+1
         WRITE (99,301) IRI,ICI,EMT(I)
    1 CONTINUE
*
  300 FORMAT (//)
  301 FORMAT (' H(DC) (',1I2,',',1I2,') = ',1PD22.15)
*
      END
