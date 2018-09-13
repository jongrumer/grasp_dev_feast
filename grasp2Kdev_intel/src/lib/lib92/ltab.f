************************************************************************
*                                                                      *
      SUBROUTINE LTAB (IS,NQS,KS,IROWS)
*                                                                      *
*   locates rows of possible parents of active shell states for acc-   *
*   essing  NTAB. It is assumed that empty shells have been elimina-   *
*   ted from consideration by SUBROUTINE RKCO.                         *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      DIMENSION IS(4),NQS(4),KS(4),IROWS(4),KQ(4)
*
      COMMON/TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      IF (IS(1) .EQ. IS(2)) NQS(1) = NQS(2)-1
      IF (IS(3) .EQ. IS(4)) NQS(3) = NQS(4)-1
*
      DO 1 I = 1,4
*
*   Check that input data are consistent
*
         IF ((NQS(I) .LE. 0) .OR. (NQS(I) .GT. KS(I))) THEN
            WRITE (*,300) NQS(I),IS(I),KS(I)
            STOP
         ENDIF
*
         KQ1 = NQS(I)-1
         KQ2 = KS(I)-KQ1
         KQ(I) = MIN (KQ1,KQ2)+1
         IF (KQ(I) .NE. 1) THEN
            IROWS(I) = (KS(I)*(KS(I)-2))/8+KQ(I)
         ELSE
            IROWS(I) = 1
         ENDIF
*
         IF (IROWS(I) .GT. NROWS) THEN
            WRITE (*,301)
            STOP
         ENDIF
*
    1 CONTINUE
*
      RETURN
*
  300 FORMAT ('LTAB: ',1I3,' Electrons in shell ',1I3,
     :        ' with 2j+1 = ',1I3)
  301 FORMAT ('LTAB: Extend COMMON block TERMS')
*
      END
