************************************************************************
*                                                                      *
      SUBROUTINE DECNSL (I,IERR)
*                                                                      *
*   This SUBROUTINE  determines twice the  orbital angular momentum,   *
*   N2LNR, given the nonrelativistic symmetry,  NRSYM(I), where I is   *
*   the index of an orbital.                                           *
*                                                                      *
*   Written by Farid A Parpia, at Oxford.   Last update: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      CHARACTER*1 NRSYM,SYMI,SYMLST
*
      DIMENSION SYMLST(9)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/NRORBN/NPNR(NNNW),N2LNR(NNNW)
     :      /NRORBS/NRSYM(NNNW)
*
      DATA SYMLST/'s','p','d','f','g','h','i','k','l'/
*
*   Initialize
*
      SYMI = NRSYM(I)
*
*   Perform comparison
*
      DO 1 IND = 1,9
*
         IF (SYMI .EQ. SYMLST(IND)) THEN
*
*   Match found:
*
*      determine twice orbital angular momentum quantum number
*
            N2LNR(I) = 2*(IND-1)
*
*      normal exit
*
            IERR = 0
            GOTO 2
*
         ENDIF
*
    1 CONTINUE
*
*   No match, error exit
*
      IERR = 1
      PRINT *, 'DECNSL: Symmetry ',SYMI,' could not be decoded.'
*
*   All done
*
    2 RETURN
      END
