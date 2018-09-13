********************************************************************
*                                                                  *
      FUNCTION JFAZE(I1,I2,I3,I4)
*                                                                  *
*     ------------  SECTION METWO    SUBPROGRAM 20  -------------  *
*                                                                  *
*     DETERMINATE THE PHASE FACTOR WHICH APPEAR FROM PERMUTATION   *
*     OF OPERATORS OF SECOND QUANTIZATION                          *
*                                                                  *
*     NO FUNCTION CALLED                                           *
*                                                                  *
********************************************************************
*
      JFAZE=1
      IF(I1.GT.I2)JFAZE=-JFAZE
      IF(I1.GT.I3)JFAZE=-JFAZE
      IF(I1.GT.I4)JFAZE=-JFAZE
      IF(I2.GT.I3)JFAZE=-JFAZE
      IF(I2.GT.I4)JFAZE=-JFAZE
      IF(I3.GT.I4)JFAZE=-JFAZE
      RETURN
      END
