************************************************************************
*
      SUBROUTINE GETMIXBLOCK(NAME,NCI)
*
*    Reads mixing coefficient file from block-structured format 
*
* Note:
*     eav is not compatible with the non-block version if some blocks
*     were not diagonalized
*
*     This is a modified version of cvtmix.f
*     written by Per Jonsson, September 2003
*
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 FILNAM
      CHARACTER*24 NAME
      CHARACTER*11 FORM
      CHARACTER*6 G92MIX
      CHARACTER*3 STATUS
*
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
      POINTER (PIATJP,IATJPO(*))
      POINTER (PIASPA,IASPAR(*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA
      COMMON/iounit/istdi,istdo,istde

*     
*   The  .mix  file is UNFORMATTED; it must exist
*
      K = INDEX(NAME,' ')
      IF (NCI .EQ. 0) THEN
         FILNAM = NAME(1:K-1)//'.cm'
      ELSE
         FILNAM = NAME(1:K-1)//'.m'
      ENDIF
      FORM = 'UNFORMATTED'
      STATUS = 'OLD'
*
      CALL OPENFL (25,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) THEN
         WRITE(istde,*) 'Error when opening',FILNAM
         STOP
      ENDIF
*
*   Check the header of the file; if not as expected, try again
*
      READ (25,IOSTAT = IOS) G92MIX
      IF ((IOS .NE. 0) .OR.
     :    (G92MIX .NE. 'G92MIX')) THEN
         WRITE(istde,*) 'Not a GRASP2K MIXing Coefficients File;'
         CLOSE (25)
         STOP
      ENDIF

      READ (25) nelec, ncftot, nw, nvectot, nvecsiz, nblock
      write (*,*) '   nelec  = ', nelec
      write (*,*) '   ncftot = ', ncftot
      write (*,*) '   nw     = ', nw
      write (*,*) '   nblock = ', nblock
      write (*,*)

************************************************************************
* Allocate memory for old format data
************************************************************************

      CALL ALLOC (PNEVAL,nvectot,8)
      CALL ALLOC (PNEVEC,ncftot*nvectot,8)
      CALL ALLOC (PNIVEC,nvectot,4)
      CALL ALLOC (PIATJP,nvectot,4)
      CALL ALLOC (PIASPA,nvectot,4)

************************************************************************
* Initialize mixing coefficients to zero; others are fine
************************************************************************
      DO i = 1, nvectot*ncftot
         evec(i) = 0.d0
      ENDDO

************************************************************************
* Initialize counters and sum registers
*
*    nvecpat:    total number of eigenstates of the previous blocks
*    ncfpat:     total number of CSF of the previous blocks
*    nvecsizpat: vector size of the previous blocks
*    eavsum:     sum of diagonal elements of the previous blocks where
*                at least one eigenstate is calculated
*    neavsum:    total number CSF contributing to eavsum
************************************************************************

      nvecpat     = 0
      ncfpat      = 0
      nvecsizpat  = 0
      neavsum     = 0
      eavsum      = 0.d0

      write (*,*) '  block     ncf     nev    2j+1  parity'
      DO jb = 1, nblock

         READ (25) nb, ncfblk, nevblk, iatjp, iaspa
			WRITE (*,'(5I8)') nb, ncfblk, nevblk, iatjp, iaspa
         IF (jb .NE. nb) STOP 'jb .NE. nb'

         IF (nevblk .GT. 0) THEN

            READ (25) (ivec(nvecpat+i), i = 1, nevblk)
            DO i = nvecpat + 1, nvecpat + nevblk
               ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
               iatjpo(i) = iatjp
               iaspar(i) = iaspa
            ENDDO

            READ (25) eav, (eval(nvecpat + i ), i = 1, nevblk)

*           ...Construct the true energy by adding up the average
            DO i = 1, nevblk
               eval(nvecpat+i) = eval(nvecpat+i) + eav
            ENDDO
*           ...For overal (all blocks) average energy
            eavsum  = eavsum + eav*ncfblk
            neavsum = neavsum + ncfblk

            READ (25) ((evec( nvecsizpat + ncfpat+i + (j-1)*ncftot ), 
     &                    i = 1, ncfblk), j = 1, nevblk)
         ENDIF

         nvecpat = nvecpat + nevblk
         ncfpat = ncfpat + ncfblk
         nvecsizpat = nvecsizpat + nevblk*ncftot

      ENDDO

*     ...Here eav is the average energy of the blocks where at least
*        one eigenstate is calculated. It is not the averge of the 
*        total Hamiltonian.

      eav = eavsum / neavsum

      IF (ncftot .NE. neavsum) THEN
         PRINT *, 'Not all blocks are diagonalized --- Average E ',
     &            'not correct'
      ENDIF

*     ...Substrct the overal average energy
      DO i = 1, nvectot
         eval(i) = eval(i) - eav
      ENDDO

      CLOSE (25)

      ncf = ncftot
      nvec = nvectot

      RETURN 
      END
