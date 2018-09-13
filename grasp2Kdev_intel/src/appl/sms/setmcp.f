************************************************************************
*                                                                      *
      SUBROUTINE SETMCP(AVAIL)
*                                                                      *
*   Open and check the  .mcp  files. File 30 stores the structure of   *
*   H(DC) ; file 31 stores the  T  coefficients;  files 32, 33, ...,   *
*   store V(0), V(1), ... .                                            *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, LENGTH, OPENFL.                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL DIAG,DIAGT,FOUND,FOUND1,LFORDR,LFORDT,AVAIL
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 FILNAM,FULNAM
      CHARACTER*3 DEFNAM
      CHARACTER*11 FORM
      CHARACTER*8 SRTLAB
      CHARACTER*3 MCPLAB,STATUS
      CHARACTER*2 CK
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /FOPARM/ICCUT
     :      /MCPA/KMAXF
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
*
*   Determine KMAXF; this is one less than the number of  .mcp
*   files for the two-electron integrals
*
      AVAIL = .TRUE.
      KMAXF = 0
      DO 1 K = 1,NW
         KMAXF = MAX (KMAXF,NKJ(K))
    1 CONTINUE
*
*   All files  grasp92.mcp.xx  are UNFORMATTED; they must exist
*
      FORM = 'UNFORMATTED'
      DEFNAM = 'mcp'
      STATUS = 'OLD'
*
*   Look for  grasp92.mcp.30 , ...
*
      FOUND = .TRUE.
      DO 2 K = 30,32+KMAXF
         CALL CONVRT (K,CK,LCK)
         INQUIRE (FILE = DEFNAM//'.'//CK(1:2), EXIST = FOUND1)
         FOUND = FOUND .AND. FOUND1
    2 CONTINUE
*
      IF (FOUND) THEN
         FILNAM = DEFNAM
      ELSE
         PRINT *, 'The mcp files do not exist'
         AVAIL = .FALSE.
         RETURN
      ENDIF
*
*   Open the files; check file headers
*
      LFN =  LEN_TRIM (FILNAM)
      DO 5 K = 30,32+KMAXF
         CALL CONVRT (K,CK,LCK)
         FULNAM = FILNAM(1:LFN)//'.'//CK(1:2)
         CALL OPENFL (K,FULNAM,FORM,STATUS,IERR)
         IF (IERR .EQ. 0) THEN
            READ (K,IOSTAT = IOS) MCPLAB,SRTLAB
            IF ((IOS .NE. 0) .OR.
     :          (MCPLAB .NE. 'MCP') .OR.
     :          (SRTLAB .NE. '  SORTED')) THEN
               PRINT *, 'Not a sorted GRASP92 MCP File;'
               IERR = IERR+1
            ENDIF
         ENDIF
         IF (IERR .EQ. 0) THEN
            READ (K) NELECT,NCFT,NWT
            IF ((NELECT .NE. NELEC) .OR.
     :          (NCFT   .NE. NCF  ) .OR.
     :          (NWT    .NE. NW   )) THEN
               PRINT *, 'Sorted GRASP92 MCP File not appropriate'
               PRINT *, '  to Configuration Symmetry List;'
               IERR = IERR+1
            ENDIF
            IF (K .EQ. 30) THEN
               READ (K) DIAG,ICCUT,LFORDR
            ELSE
               READ (K) DIAGT,ICCUTT,LFORDT
               IF ((DIAGT .NEQV. DIAG) .OR.
     :             (ICCUTT .NE. ICCUT) .OR.
     :             (LFORDT .NEQV. LFORDR)) THEN
                  PRINT *, 'Sorted GRASP92 MCP Files are not'
                  PRINT *, ' consistent;'
                  IERR = IERR+1
               ENDIF
            ENDIF
         ENDIF
         IF (IERR .NE. 0) THEN
            DO 4 I = 30,K
               CLOSE (I)
    4       CONTINUE
            AVAIL = .FALSE.
            RETURN
         ENDIF
    5 CONTINUE
*
      RETURN
      END
