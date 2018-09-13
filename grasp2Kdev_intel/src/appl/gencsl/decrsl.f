************************************************************************
*                                                                      *
      SUBROUTINE DECRSL (I,IERR)
*                                                                      *
*   This subroutine determines twice the total angular momentum, N2J,  *
*   and the orbital angular momentum, NL, given the symmetry, SYM(I),  *
*   where I is the index of an orbital.                                *
*                                                                      *
*   Written by Farid A Parpia, at Oxford.   Last update: 24 Aug 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*2 SYM,SYMI,SYMLST
*
      DIMENSION SYMLST(16)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/ORBNUM/NP(NNNW),N2J(NNNW),NL(NNNW)
     :      /ORBSYM/SYM(NNNW)
*
      DATA SYMLST/'s ','p-','p ','d-','d ','f-','f ','g-','g ',
     :                 'h-','h ','i-','i ','k-','k ','l-'/
*
*   Initialize
*
      SYMI = SYM(I)
*
*   Perform comparison
*
      DO 1 IND = 1,16
*
         IF (SYMI .EQ. SYMLST(IND)) THEN
*
*   Match found:
*
*      determine orbital angular momentum quantum number
*
            NL(I) = IND/2
*
*      determine twice total angular momentum quantum number
*
            IF (MOD (IND,2) .EQ. 0) THEN
               N2J(I) = IND-1
            ELSE
               N2J(I) = IND
            ENDIF
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
      PRINT *, 'DECRSL: Symmetry ',SYMI,' could not be decoded.'
*
*   All done
*
    2 RETURN
      END
