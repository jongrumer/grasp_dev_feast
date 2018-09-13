************************************************************************
*                                                                      *
      SUBROUTINE GETINF
*                                                                      *
*   Interactively determines data governing the generation of MCP co-  *
*   efficients.                                                        *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, GETYN.                                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
*   Modified Xinghong He                  Last revision: 30 Jun 1998   *
*                                                                      *
*   File shared (hard link) by mcpmpi, mcpblk
*
*   Updated to treat ICCUT for block
*
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG,GETYN,LFORDR,YES
      CHARACTER*20 CNUM
*
      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      COMMON/DEFAULT/NDEF
     :      /FOPARM/ICCUT(100)
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
*
*   Determine the physical effects specifications
*
!     Comment out the EAL option
!     IF (NDEF .NE. 0) THEN
!        WRITE  (istde,*) 'Generate MCP coefficients only for'
!    & , ' diagonal matrix elements? '
!        WRITE (istde,*) '(This is appropriate to (E)AL calculation):'
!        DIAG = GETYN ()
!     ELSE
!        DIAG = .FALSE.
!     ENDIF 
      DIAG = .FALSE.
      IF (DIAG) THEN
         LFORDR = .FALSE.
         do i = 1,100
            ICCUT(i) = 0
         end do
      ELSE
         IF (NDEF .NE. 0) THEN
            WRITE (istde,*) 'Treat contributions of some CSFs'
     &,              ' as first-order perturbations?'
            LFORDR = GETYN ()
         ELSE
            LFORDR = .FALSE.
         ENDIF
         IF (LFORDR) THEN
            WRITE (istde,*) 'The contribution of CSFs 1 -- ICCUT will'
     &,              ' be treated variationally;'
            WRITE (istde,*) 'the remainder perturbatively; enter ICCUT:'
            do i = 1,nblock
               write(istde,*) 'Give ICCUT for block',i
    1          READ *, ICCUT(i)
               IF ((ICCUT(i).LE.1).OR.(ICCUT(i).GE.ncfblk(i))) THEN
                  WRITE (istde,*) 'GETINF: ICCUT must be greater than 1'
     &,                 ' and less than ',ncfblk(i)
                  WRITE (istde,*) ' please reenter ICCUT:'
                  GOTO 1
               ENDIF
            end do
         ENDIF
      ENDIF

      RETURN
      END
