************************************************************************
*                                                                      *
      SUBROUTINE PRTRSL
*                                                                      *
*   This subroutine is now called only in GETOLD to print out all 
*   orbitals so that the user can do cut-and-paste in supplying the
*   orbitals to be varied. I.e., LFIX(LOC) will always be FALSE and
*   are thus commented out.
* Xinghong HE 6 Apr 1998
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LFIX
      CHARACTER*80 RECORD
      CHARACTER*2 CNUM,NH
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
!     :      /FIXD/NFIX,LFIX(NNNW)
!     :      /ORBA/IORDER(NNNW)
      COMMON/iounit/istdi,istdo,istde
*
*   Print the list of all subshell radial wavefunctions
*
      IEND = 0
      DO 1 I = 1,NW
!         LOC = IORDER(I)
         LOC = I
! For the commenting-out of the IF...ENDIF, see comments in the header
!         IF (.NOT. LFIX(LOC)) THEN
            IF (IEND .GT. 75) THEN
               WRITE (istde,*) RECORD(1:IEND)
               IEND = 0
            ENDIF
            IF (IEND .GT. 0) THEN
               IBEG = IEND+1
               IEND = IBEG
               RECORD(IBEG:IEND) = ' '
            ENDIF
            IBEG = IEND+1
            CALL CONVRT (NP(LOC),CNUM,LENTH)
            IEND = IBEG+LENTH-1
            RECORD(IBEG:IEND) = CNUM(1:LENTH)
            IF (NAK(LOC) .LT. 0) THEN
               LENTH = 1
            ELSE
               LENTH = 2
            ENDIF
            IBEG = IEND+1
            IEND = IBEG+LENTH-1
            RECORD(IBEG:IEND) = NH(LOC)(1:LENTH)
!         ENDIF
    1 CONTINUE
      IF (IEND .GT. 1) WRITE (istde,*) RECORD(1:IEND)
*
      RETURN
      END
