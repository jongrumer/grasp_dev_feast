************************************************************************
*                                                                      *
      SUBROUTINE FNDBEG (JASTRT,JBSTRT,INDEX,LLISTT,LLISTV)
*                                                                      *
*   Determines the pair of CSFs with which to restart the generation   *
*   of MCP coefficients; positions all files appropriately.            *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 07 Dec 1992   *
*                                                                      *
*   Modified by Per Jonsson to be used with the SMS92 program          * 
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*8 SRTLAB
      CHARACTER*3 MCPLAB
*
      DIMENSION LLISTV(0:*)
*
      COMMON/MCPA/KMAX
     :      /ORB2/NCF,NW,PNTRIQ
*
*   Begin by examining file 30; this is the last to be updated by
*   SUBROUTINE MCP
*
*   Read and check the character part of the file header
*
      REWIND (30)
      READ (30) MCPLAB,SRTLAB
      IF (SRTLAB .EQ. '  SORTED') THEN
         JASTRT = NCF+1
         JBSTRT = NCF+1
         LLISTT = 0
         DO 1 K = 0,KMAX
            LLISTV(K) = 0
    1    CONTINUE
         GOTO 8
      ENDIF
*
      READ (30)
      READ (30)
*
*   Read as many records as possible
*
    2 READ (30,IOSTAT = IOS) ICREAD,IRREAD,INDX
*
      IF (IOS .EQ. 0) THEN
*
*   No errors or end-of-file; keep reading
*
         JASTRT = ICREAD
         JBSTRT = IRREAD
         INDEX = INDX
         GOTO 2
*
      ELSE
*
         IF ((JASTRT .EQ. NCF) .AND. (JBSTRT .EQ. NCF)) THEN
*
*   All coefficients have been generated; sorting may still
*   be necessary; force this option
*
            JASTRT = NCF+1
            JBSTRT = NCF+1
*
         ELSE
*
*   Some coefficients remain to be generated; reposition all files
*   for augmentation of lists by SUBROUTINE MCP; update JBSTRT and,
*   if appropriate, JASTRT
*
            DO 6 K = 31,32+KMAX
               REWIND (K)
               NREC = 3
               DO 3 I = 1,NREC
                  READ (K)
    3          CONTINUE
    4          READ (K,IOSTAT = IOS) INDX,LABEL,COEFF
               IF ((IOS .EQ. 0) .AND. (INDX .LE. INDEX)) THEN
                  NREC = NREC+1
                  GOTO 4
               ELSE
                  REWIND (K)
                  DO 5 I = 1,NREC
                     READ (K)
    5             CONTINUE
                  IF (K .GT. 31) THEN
                     LLISTV(K-32) = NREC-3
                  ELSE
                     LLISTT = NREC-3
                  ENDIF
               ENDIF
    6       CONTINUE
*
*  Now, reposition the sms file. This file should contain the
*  same number of data records as file 33.

            REWIND (20)
            DO 7 I = 1, LLISTV(1)
               READ (20)
    7       CONTINUE

            JBSTRT = JBSTRT+1
            IF (JBSTRT .GT. NCF) THEN
               JASTRT = JASTRT+1
               JBSTRT = JASTRT
            ENDIF
*
         ENDIF
*
      ENDIF
*
    8 RETURN
      END
