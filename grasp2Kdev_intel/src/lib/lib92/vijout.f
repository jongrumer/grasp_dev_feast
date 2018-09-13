************************************************************************
*                                                                      *
      SUBROUTINE VIJOUT (JA,JB)
*                                                                      *
*   Prints  out tables of configurational quantum numbers defined by   *
*   SETQNA for current matrix element.                                 *
*                                                                      *
*                                           Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*2 IC,IJ
*
      DIMENSION IJ(2),JC(2),IC(2)
*
CGG      INTEGER NNNW
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      COMMON/M0/JJC1(NNNW),JJC2(NNNW)
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
*
      DATA IJ/'/2','  '/
*
      IF (NPEEL .GT. 0) THEN
*
*   Identify CSFs
*
         WRITE (99,300) JA,JB
*
*   Print active shell quantum numbers from JLIST table
*
         WRITE (99,301)
*
         DO 1 J = 1,NPEEL
*
            JW = JLIST(J)
*
            JC(1) = JJQ1(3,JW)-1
            IF (MOD (JC(1),2) .EQ. 1) THEN
               IC(1) = IJ(1)
            ELSE
               JC(1) = JC(1)/2
               IC(1) = IJ(2)
            ENDIF
*
            JC(2) = JJQ2(3,JW)-1
*
            IF (MOD (JC(2),2) .EQ. 1) THEN
               IC(2) = IJ(1)
            ELSE
               JC(2) = JC(2)/2
               IC(2) = IJ(2)
            ENDIF
*
            WRITE (99,302) JW,NQ1(JW),JJQ1(1,JW),JJQ1(2,JW),JC(1),IC(1)
     :                       ,NQ2(JW),JJQ2(1,JW),JJQ2(2,JW),JC(2),IC(2)
    1    CONTINUE
*
*   Print coupling angular momenta if NPEEL .GE. 2
*
         IF (NPEEL .GT. 2) THEN
*
            WRITE (99,303)
            DO 2 J = 2,NPEEL
*
               JC(1) = JJC1(J-1)-1
*
               IF (MOD (JC(1),2) .EQ. 1) THEN
                  IC(1) = IJ(1)
               ELSE
                  JC(1) = JC(1)/2
                  IC(1) = IJ(2)
               ENDIF
*
               JC(2) = JJC2(J-1)-1
               IF (MOD (JC(2),2) .EQ. 1) THEN
                  IC(2) = IJ(1)
               ELSE
                  JC(2) = JC(2)/2
                  IC(2) = IJ(2)
               ENDIF
*
               WRITE (99,304) (JC(I),IC(I),I = 1,2)
    2       CONTINUE
         ENDIF
*
      ENDIF
*
      WRITE (99,305) NCORE
*
      RETURN
*
  300 FORMAT (/'From VIJOUT: CSF ',1I2,35X,'CSF ',1I2)
  301 FORMAT (3X,'subshell',4X,'q',4X,'v',2X,'w',2X,'J',
     :                     19X,'q',4X,'v',2X,'w',2X,'J')
  302 FORMAT (7X,I3,I6,I5,2I3,A2,15X,I3,I5,2I3,A2)
  303 FORMAT (' coupling schemes:')
  304 FORMAT (14X,I2,A2,27X,I2,A2)
  305 FORMAT (' there are ',I3,' inactive closed shells.'/)
*
      END
