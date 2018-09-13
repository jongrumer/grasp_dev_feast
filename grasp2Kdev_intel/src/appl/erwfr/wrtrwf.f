************************************************************************
*                                                                      *
      SUBROUTINE WRTRWF
*                                                                      *
*   Open, write a header and all subshell radial wavefunctions, and    *
*   close the  .rwf  file.                                             *
*                                                                      *
*   Call(s) to: [LIB92]: DALLOC, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*128 FILNAM
      !CHARACTER*11 FORM
      !CHARACTER*3 STATUS
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)

      COMMON/iounit/istdi,istdo,istde
      LOGICAL Imopened

      !FORM = 'UNFORMATTED'
      !STATUS = 'NEW'
      FILNAM = 'rwfn.inp'

      DO newunit = 23, 99             ! 23 is a historical value
         INQUIRE (UNIT=newunit, OPENED=Imopened)
         IF ( .NOT. Imopened) EXIT      ! should be the normal exit point
      ENDDO

      IF (newunit .EQ. 100) THEN
         WRITE (istde,*) 'All unit numbers from 23 to 99 are BUSY!'
         STOP
      ENDIF

      CALL OPENFL (newunit,FILNAM,'UNFORMATTED','NEW',IERR)
      IF (IERR .EQ. 1) THEN
         WRITE (istde,*) 'Error when opening "'
     &                  , filnam(1:LEN_TRIM (filnam)), '"'
         
         STOP
      ENDIF 
*
*   Write the file header
*
      WRITE (newunit) 'G92RWF'
*
*   Write out the radial wavefunctions
*
      DO 2 J = 1,NW
         MFJ = MF(J)
         WRITE (newunit) NP(J),NAK(J),E(J),MFJ
         WRITE (newunit) PZ(J),(PF(I,J),I = 1,MFJ),(QF(I,J),I = 1,MFJ)
         WRITE (newunit) (R(I),I = 1,MFJ)
    2 CONTINUE

      CLOSE (newunit)
*
*   Deallocate the storage for the radial wavefunctions
*
      CALL DALLOC (PNTRPF)
      CALL DALLOC (PNTRQF)
*
      RETURN
      END
