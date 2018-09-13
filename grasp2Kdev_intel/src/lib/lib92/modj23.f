************************************************************************
*                                                                      *
      SUBROUTINE MODJ23
*                                                                      *
*   Restores  COMMON  block  /COUPLE/ from saved values for exchange   *
*   case.                                                              *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      PARAMETER (MANGM = 60,MTRIAD = 12)
*
      LOGICAL FREE
*
      COMMON/L2/J2S(MTRIAD,3),J3S(MTRIAD,3)
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
*
      NS2 = NJA-1
      DO 2 J = 1,3
         DO 1 I = 1,NS2
            J2(I,J) = J2S(I,J)
            J3(I,J) = J3S(I,J)
    1    CONTINUE
    2 CONTINUE
*
      I = J3(1,3)
      J3(1,3) = J2(1,1)
      J2(1,1) = I
*
      RETURN
      END
