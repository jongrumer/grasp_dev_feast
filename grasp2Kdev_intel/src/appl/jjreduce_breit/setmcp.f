************************************************************************
*                                                                      *
      SUBROUTINE SETMCP (RESTRT)
*                                                                      *
*   Open and check the  .mcp  files. File 30 stores the structure of   *
*   H(DC) ; file 31 stores the  T  coefficients;  files 32, 33, ...,   *
*   store V(0), V(1), ... .                                            *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, GETYN, LENGTH, LODRES, OPENFL.        *
*               [GENMCP]: GETINF.                                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
! LENGTH --> LEN_TRIM
! XHH 1997.02.13
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL FOUND,FOUND1,GETYN,RESTRT,YES
Cww      INTEGER PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(*))
      CHARACTER*256 FILNAM,FULNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 MCPRES,STATUS
      CHARACTER*2 CK
*
      COMMON/MCPA/KMAX
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /DEFAULT/NDEF
*
*   Determine KMAX; this is the number of  .mcp  files for the
*   two-electron integrals
*
      KMAX = 0
      DO 1 K = 1,NW
         KMAX = MAX (KMAX,NKJ(K))
    1 CONTINUE
*
*   All files  grasp92.mcp.xx  are UNFORMATTED;
*
      FORM = 'UNFORMATTED'
      DEFNAM = 'mcp'
*
*   Determine if this is a restart
*
      IF (NDEF.EQ.1) THEN
         PRINT *, 'Restarting GENMCP ?'
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF
*
      IF (YES) THEN
*
*   Assume that the computation can be restarted
*
         RESTRT = .TRUE.
         STATUS = 'OLD'
*
*   Look for  mcp.30 , ...
*
         FOUND = .TRUE.
         DO 2 K = 30,32+KMAX
            CALL CONVRT (K,CK,LCK)
            INQUIRE (FILE = DEFNAM//'.'//CK(1:2), EXIST = FOUND1)
            FOUND = FOUND .AND. FOUND1
    2    CONTINUE
*
         IF (FOUND) THEN
*
*   All files  mcp.xx  exist; ascertain that these are to
*   be used; if not, determine another root filename
*
            PRINT *, 'Files  mcp  found; enter another'
            PRINT *, ' root file name if these files are not'
            PRINT *, ' to be used as the GRASP92 MCP Files;'
            PRINT *, ' null otherwise:'
            READ (*,'(A)') FILNAM
*
            IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = DEFNAM
*
         ELSE
*
    3       PRINT *, 'Enter the root name for the GRASP92'
            PRINT *, ' MCP Files:'
            READ (*,'(A)') FILNAM
*
            IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 3
*
         ENDIF
*
      ELSE
*
*   Not a restart
*
         RESTRT = .FALSE.
         STATUS = 'NEW'
         FILNAM = DEFNAM
*
      ENDIF
*
      LFN = LEN_TRIM (FILNAM)
      DO 6 K = 30,32+KMAX
         CALL CONVRT (K,CK,LCK)
         FULNAM = FILNAM(1:LFN)//'.'//CK(1:2)
         CALL OPENFL (K,FULNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
            DO 5 I = 30,K
               CLOSE (I)
    5       CONTINUE
            IF (RESTRT) THEN
               GOTO 3
            ELSE
               PRINT *, 'Error when opening the mcp files'
               STOP
            ENDIF
         ENDIF
    6 CONTINUE
*
      IF (RESTRT) THEN
*
*   Check file headers if restarting
*
         DO 8 K = 30,32+KMAX
            READ (K,IOSTAT = IOS) MCPRES
            IF ((IOS .NE. 0) .OR. (MCPRES .NE. 'MCP')) THEN
               PRINT *, 'Not a GRASP92 MCP File;'
               DO 7 I = 30,32+KMAX
                  CLOSE (I)
    7          CONTINUE
               GOTO 3
            ENDIF
    8    CONTINUE
*
*   Determine the data that control the calculation
*
         CALL LODRES (IERR)
*
         IF (IERR .NE. 0) THEN
            DO 9 I = 30,32+KMAX
               CLOSE (I)
    9       CONTINUE
            GOTO 3
         ENDIF
*
*   Write file headers if not restarting
*
      ELSE
*
         DO 10 K = 30,32+KMAX
            WRITE (K) 'MCP','UNSORTED'
   10    CONTINUE
*
*   Determine the data that control the calculation; add to
*   the file headers
*
         CALL GETINF
*
      ENDIF
*
      RETURN
      END
