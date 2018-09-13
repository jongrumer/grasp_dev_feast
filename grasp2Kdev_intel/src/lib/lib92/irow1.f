************************************************************************
*                                                                      *
      FUNCTION IROW1 (NELC,KSI)
*                                                                      *
*   Locate the row position of configuration j(**n) in table NTAB.     *
*                                                                      *
*                                           Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IF ((NELC .LE. 0) .OR. (NELC .GT. KSI)) THEN
         WRITE (*,300) NELC,KSI
         STOP
      ENDIF
*
      KQ1 = NELC-1
      KQ2 = KSI-KQ1
      KQL = MIN (KQ1,KQ2)+1
      IF (KQL .EQ. 1) THEN
         IROW1 = 1
      ELSE
         IROW1 = (KSI*(KSI-2))/8+KQL
      ENDIF
*
      RETURN
*
  300 FORMAT ('IROW1: ',I3,' electrons in shell with 2j+1 = ',I3)
*
      END
