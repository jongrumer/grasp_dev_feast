      PROGRAM cvtmix   
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
************************************************************************
*
* Purpose:
*    Convert a mixing coefficient file from block-structured format to
*    the old format
*
* Routines called:
*    alloc
*
* Note:
*     eav is not compatible with the non-block version if some blocks
*     were not diagonalized
*
* Xinghong He
*
************************************************************************
      POINTER (pntivec, ivec(1))
      POINTER (pntiatjpo, iatjpo(1))
      POINTER (pntiaspar, iaspar(1))
      POINTER (pnteval, eval(1))
      POINTER (pntevec, evec(1))

      CHARACTER*24, infile, outfile, g92mix*6
      DATA nin, nout /25, 26/

************************************************************************
* Format of the block structre data
*
*      READ (nin) g92mix
*      READ (nin) nelec, ncftot, nw, nvectot, nvecsiz, nblock
*      DO jb = 1, nblock
*         READ (nin) nb, ncfblk, nevblk, iatjp, iaspa
*         IF (nevblk .GT. 0) THEN
*            READ (nin) (iccmin(i),i=1,nevblk)
*            READ (nin) EAV,(EVAL(I), I = 1, nevblk)
*            READ (nin) ((EVEC(I+(J-1)*NCFblk), I=1,NCFblk), J=1,nevblk)
*         ENDIF
*      ENDDO
************************************************************************
      PRINT *
      PRINT *, 'Converting mixing coefficient file from block strucure'
     &       , ' to old format'
      PRINT *, 'The input and the output files may have the same name'
      PRINT *, '(overwritten mode, not recommended)'
      PRINT *, 'The input file:'
      READ (*,'(A)') infile
      PRINT *, 'The output file:'
      READ (*,'(A)') outfile
      
      OPEN (nin, FILE = infile, FORM = 'UNFORMATTED', STATUS = 'OLD'
     &      , IOSTAT = IOS)
      IF (IOS .NE. 0) STOP 'Failed to to open input file'

      READ (nin) g92mix
      IF (g92mix .NE. 'G92MIX') 
     &   STOP 'Not a mixing coefficient file'

      READ (nin) nelec, ncftot, nw, nvectot, nvecsiz, nblock
      write (*,*) '   nelec  = ', nelec
      write (*,*) '   ncftot = ', ncftot
      write (*,*) '   nw     = ', nw
      write (*,*) '   nblock = ', nblock
      write (*,*)

************************************************************************
* Allocate memories for old format data
************************************************************************
      CALL alloc (pntivec,   nvectot, 4)
      CALL alloc (pntiatjpo, nvectot, 4)
      CALL alloc (pntiaspar, nvectot, 4)
      CALL alloc (pnteval,   nvectot, 8)
      CALL alloc (pntevec,   nvectot*ncftot, 8)
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
      neavsum = 0
      eavsum = 0.d0

      write (*,*) '  block     ncf     nev    2j+1  parity'
      DO jb = 1, nblock

         READ (nin) nb, ncfblk, nevblk, iatjp, iaspa
			WRITE (*,'(5I8)') nb, ncfblk, nevblk, iatjp, iaspa
         IF (jb .NE. nb) STOP 'jb .NE. nb'

         IF (nevblk .GT. 0) THEN

            READ (nin) (ivec(nvecpat+i), i = 1, nevblk)
            DO i = nvecpat + 1, nvecpat + nevblk
               ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
               iatjpo(i) = iatjp
               iaspar(i) = iaspa
            ENDDO

            READ (nin) eav, (eval(nvecpat + i ), i = 1, nevblk)

*           ...Construct the true energy by adding up the average
            DO i = 1, nevblk
               eval(nvecpat+i) = eval(nvecpat+i) + eav
            ENDDO
*           ...For overal (all blocks) average energy
            eavsum  = eavsum + eav*ncfblk
            neavsum = neavsum + ncfblk

            READ (nin) ((evec( nvecsizpat + ncfpat+i + (j-1)*ncftot ), 
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

      CLOSE (nin)

************************************************************************
* Write-out is straight forward since converting has been done in 
* accordance with the output format below
************************************************************************
      ncf = ncftot
      nvec = nvectot

      OPEN (nout, FILE = outfile, FORM = 'UNFORMATTED',
     &      STATUS = 'UNKNOWN', IOSTAT = IOS)
      IF (IOS .NE. 0) STOP 'Failed to to open output file'

      WRITE (nout) g92mix
      WRITE (nout) nelec, ncf, nw
      WRITE (nout) nvec
      WRITE (nout) (ivec(i), i = 1, nvec)
      WRITE (nout) (iatjpo(i), iaspar(i), i = 1, nvec)
      WRITE (nout) eav, (eval(i), i = 1, nvec)
      WRITE (nout) ((evec(i+(j-1)*ncf), i = 1, ncf), j = 1, nvec)
!            print '(4i)',
!     :  nelec, ncf, nw,nvec
!            print '(3i)',
!     :   (ivec(i), i = 1, nvec)
!            print '(3i)',
!     :   (iatjpo(i), iaspar(i), i = 1, nvec)
!            print '(4f20.15)',
!     :     eav, (eval(i), i = 1, nvec)
!            print '(4f20.15)',
!     :     ((evec(i+(j-1)*ncf), i = 1, ncf), j = 1, nvec)
      CLOSE (nout) 

      CALL dalloc (pntivec)
      CALL dalloc (pntiatjpo)
      CALL dalloc (pntiaspar)
      CALL dalloc (pnteval)
      CALL dalloc (pntevec)

      STOP 'Converting successfully finished'
      END
