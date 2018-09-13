************************************************************************
*                                                                      *
      FUNCTION IUNPCK (IPACKD,ISUBSH)
*                                                                      *
*   It is assumed that IPACKD is an integer vector of minimum length   *
*   NINT (NW/4) on 32-bit machines, and NINT (NW/8) on 64-bit machi-   *
*   nes; this cannot be checked from the argument list.   ISUBSH  is   *
*   the index of a subshell.  Array  IPACKD is loaded by  SUBROUTINE   *
*   PACK; the header of the latter should be consulted to understand   *
*   the operation of the present subprogram.                           *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
*                                                                      *
************************************************************************
*
! To save space on CRAY.
! 1997.02.12
      INTEGER*4 IPACKD
      DIMENSION IPACKD(1:*)
*
      ISUBM1 = ISUBSH-1
      LOC1 = ISUBM1/4+1
      LOC2 = MOD (ISUBM1,4)
*
      IUNPCK = IPACKD (LOC1)
*
      DO 1 I = 1,LOC2
         IUNPCK = IUNPCK/256
    1 CONTINUE
*
      IUNPCK = IUNPCK-(IUNPCK/256)*256
*
      RETURN
      END
