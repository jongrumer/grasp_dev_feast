************************************************************************
*                                                                      *
      SUBROUTINE LODCSH (nfile, NCORE)
*                                                                      *
*   Loads the data from the  .csl  file. A number of checks are made   *
*   to ensure correctness and consistency.                             *
*                                                                      *
*   Call(s) to: [LIB92]:  PRSRCN, PRSRSL       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*   Modified by C. F. Fischer to read only the header information
*
* Input:
*   nfile
*
* Output:
*   ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      CHARACTER*256 str
      CHARACTER*2 NH
*
      DIMENSION IOCC(NNNW)
*
      INTEGER*4 IQAdum

      POINTER (PNTRIQ,IQAdum)
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)

      COMMON/iounit/istdi,istdo,istde
*
*   Entry message
*
      PRINT *, 'Loading CSF file ... Header only'
*
*   Get the list of subshells
*
      NW = 0
*
*   Read the list of core subshells; set up the arrays NP, NAK,
*   NKL, NKJ, NH for these subshells
*
      CALL PRSRSL (nfile,1)
      NCORE = NW
*
*   Skip the peel subshell identification header; read the list of
*   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
*   these subshells
*
      READ (nfile,*)
      CALL PRSRSL (nfile,2)
      NPEEL = NW-NCORE
*
*   Ensure that the sets of core and peel subshell are disjoint
*
      DO 2 J = NCORE+1,NW
         NPJ = NP(J)
         NAKJ = NAK(J)
         DO 1 I = 1,NCORE
            IF ((NP(I) .EQ. NPJ) .AND. (NAK(I) .EQ. NAKJ)) THEN
               WRITE(istde,*) 'lodcsh: The lists of core and'
     &,                ' peel subshells must form disjoint sets.'
               STOP
            ENDIF
    1    CONTINUE
    2 CONTINUE

      PRINT *, 'There are/is ', nw, ' relativistic subshells;'
*
*   Skip the header for the list of CSFs
*
      READ (nfile,*)
*
*   To determine the number of electrons. This was done very much later
*   in the non-block mcp program( near the end of lodcsh). In the block
*   version, that will be used as a check to the value obtained below.
*
      READ (nfile,'(A)',IOSTAT = IOS) str
      CALL PRSRCN (str,NCORE,IOCC,IERR)
      BACKSPACE (nfile)	! return to the first CSF item

      nelec = 0
*             Number of electrons in the peel shells
      DO i = ncore + 1, nw
         nelec = nelec + iocc(i)
      ENDDO

*             Add the number of electrons in the core shells
      DO i = 1, ncore
         nelec = nelec + nkj(i) + 1
      ENDDO
      
      RETURN
      END
