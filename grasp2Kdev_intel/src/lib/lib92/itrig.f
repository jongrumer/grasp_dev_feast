************************************************************************
*                                                                      *
      FUNCTION ITRIG (I1,I2,I3)
*                                                                      *
*   The  triangular delta. Input: Values of 2*J+1; Output: 1, IF J'S   *
*   form a triangle; 0, otherwise.                                     *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      I4 = I2-I3
      IF ((I1 .GE. (ABS (I4)+1)) .AND. ((I1 .LE. (I2+I3-1)))) THEN
         ITRIG = 1
      ELSE
         ITRIG = 0
      ENDIF
*
      RETURN
      END
