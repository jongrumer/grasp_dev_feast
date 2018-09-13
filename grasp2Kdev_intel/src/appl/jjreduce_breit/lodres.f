************************************************************************
*                                                                      *
      SUBROUTINE LODRES (IERR)
*                                                                      *
*   Loads and checks certain data from the  .res  files.               *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, GETINF, LENGTH, LODRES, OPENFL.       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 07 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL DIAG,DIAGT,LFORDR,LFRDRT
Cww      INTEGER PNTRIQ
      POINTER(PNTRIQ,RIQDUMMY(*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /MCPA/KMAX
     :      /FOPARM/ICCUT
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
*
      IERR = 0
*
*   The second and third record in each file should be identical
*
      READ (30) NELECT,NCFT,NWT
      IF (NELECT .NE. NELEC) THEN
         PRINT *, 'LODRES: The number of electrons does not'
         PRINT *, ' match that from the Configuration Symmetry'
         PRINT *, ' List File;'
         IERR = IERR+1
      ENDIF
      IF (NCFT .NE. NCF) THEN
         PRINT *, 'LODRES: The number of CSFs does not'
         PRINT *, ' match that from the Configuration'
         PRINT *, ' Symmetry List File;'
         IERR = IERR+1
      ENDIF
      IF (NWT .NE. NW) THEN
         PRINT *, 'LODRES: The number of subshells does not'
         PRINT *, ' match that from the Configuration'
         PRINT *, ' Symmetry List File;'
         IERR = IERR+1
      ENDIF
*
      READ (30) DIAG,ICCUT,LFORDR
      IF (LFORDR) THEN
         IF ((ICCUT .LT. 1) .OR. (ICCUT .GT. NCF)) THEN
            PRINT *, 'LODRES: The first-order calculation'
            PRINT *, ' determined from the GENMCP REStart'
            PRINT *, ' files is not appropriate to the'
            PRINT *, ' Configuration Symmetry List File;'
            IERR = IERR+1
         ENDIF
      ENDIF
*
      DO 1 K = 31,32+KMAX
         READ (K) NELECT,NCFT,NWT
         IF ((NELEC .NE. NELECT) .OR.
     :       (NCF .NE. NCFT) .OR.
     :       (NWT .NE. NW)) THEN
            PRINT *, 'LODRES: Inconsistent electron cloud data'
            PRINT *, ' detected on comparing GENMCP REStart'
            PRINT *, ' Files;'
            IERR = IERR+1
         ENDIF
         READ (K) DIAGT,ICCUTT,LFRDRT
         IF ((DIAG .NEQV. DIAGT) .OR.
     :       (ICCUT .NE. ICCUTT) .OR.
     :       (LFORDR .NEQV. LFRDRT)) THEN
            PRINT *, 'LODRES: Inconsistent calculation type'
            PRINT *, ' detected on comparing GENMCP REStart'
            PRINT *, ' Files;'
            IERR = IERR+1
         ENDIF
    1 CONTINUE
*
      RETURN
      END
