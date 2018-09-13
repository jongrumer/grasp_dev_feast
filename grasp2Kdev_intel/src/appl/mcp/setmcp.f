************************************************************************
*                                                                      *
      SUBROUTINE SETMCP (myid, nprocs, ncore, idblk, filehead)
      IMPLICIT REAL*8          (A-H,O-Z)
*                                                                      *
*   Open and check the  .mcp  files. File 30 stores the structure of   *
*   H(DC) ; file 31 stores the  T  coefficients;  files 32, 33, ...,   *
*   store V(0), V(1), ... .                                            *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, GETYN, OPENFL.                        *
*               [GENMCP]: GETINF.                                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 08 Dec 1992   *
*   MPI version by Xinghong He            Last revision: 30 Jun 1998   *
*
*   Used by mcpvu, mcpmpivu
*                                                                      *
************************************************************************

      INTEGER   myid, nprocs, ncore
      CHARACTER idblk(*)*8, filehead*(*)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      LOGICAL FOUND,FOUND1,GETYN,YES
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER CK*2

! Current ncf in the common is ncftot <-- obtained from setcsl/setcsll
      COMMON/MCPA/KMAX
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /DEFAULT/NDEF

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
*   Determine KMAX; this is the number of  .mcp  files for the
*   two-electron integrals

      KMAX = 0
      DO K = 1, NW
         KMAX = MAX (KMAX,NKJ(K))
      ENDDO

*   All files  mcp.xx  are UNFORMATTED;

      lng = LEN_TRIM (filehead)
      DO K = 30, 32 + KMAX
         CALL CONVRT (K,CK,LCK)
         CALL OPENFL (K, filehead(1:lng) // '.' // CK(1:2), 
     &                  'UNFORMATTED', 'UNKNOWN', IERR)
         IF (IERR .NE. 0) THEN
             DO I = 30, K
                CLOSE (I)
             ENDDO
             WRITE (istde,*) 'Error when opening the mcp files'
             STOP
          ENDIF
      ENDDO
*
*   We want to know kmax before openning other mcp files (not mcp.30) 
*   in rscf
*
      WRITE (30) ncore, nblock, kmax
      WRITE (30) (ncfblk(i), i=1, nblock)
      WRITE (30) (idblk(i), i=1, nblock)

      DO K = 30, 32+KMAX
         WRITE (K) 'MCP', nblock, myid, nprocs
      ENDDO

      RETURN
      END
