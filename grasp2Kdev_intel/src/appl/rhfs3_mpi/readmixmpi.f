************************************************************************
*                                                                      *
      SUBROUTINE READMIXmpi (NAME,INPCI)
*                                                                      *
*   Open and read the mixing coefficient file                          *
*                                                                      *
*   Written by Per Jonsson                                             *
cb adapted for  hfs92MPI by Jacek Bieron
cb 14 April 2008
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

Cww      INTEGER PNTRIQ, PNTIQR
      POINTER(PNTRIQ,RIQDUMMY)

      CHARACTER*24 NAME
c     CHARACTER*128 NAME(2)
      CHARACTER*6 G92MIX

      POINTER(PNEVAL,EVAL(*))
      POINTER(PNEVEC,EVEC(*))
      POINTER(PNIVEC,IVEC(*))
      POINTER(PIATJP,IATJPO(*)),(PIASPA,IASPAR(*))

*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /PRNT/NVEC,PNIVEC,NVECMX
     :      /SYMA/PIATJP,PIASPA

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr

c     print *, ' readmixmpi: myid, NAME = ', myid, NAME

*
*   Read the initial state mixing file    
*
      if (myid .eq. 0) then

      J = INDEX(NAME,' ')
      IF (INPCI.EQ.0) THEN
        OPEN (UNIT = 68,FILE=NAME(1:J-1)//'.cm',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ELSE
        OPEN (UNIT = 68,FILE=NAME(1:J-1)//'.m',FORM='UNFORMATTED',
     :        STATUS='OLD')
      ENDIF
      READ(68,IOSTAT=IOS) G92MIX
      write(*,*) G92MIX
      IF ((IOS.NE.0).OR.(G92MIX.NE.'G92MIX')) THEN
        WRITE(*,*) 'Not a GRASP mixing file'
        STOP
      ENDIF

!     READ(68) N11,N12,N13
!     write(*,*) N11,N12,N13
!     IF ((N11.NE.NELEC).OR.(N12.NE.NCF).OR.(N13.NE.NW)) THEN
!       PRINT *, 'This MIXing Coefficients File is not'
!       PRINT *, 'appropriate for the initial state'
!       STOP
!     ENDIF
*
!     READ(68) NVEC
!     write(*,*) 'nvec',NVEC
      READ (68) nelec, ncftot, nw, nvectot, nvecsiz, nblock
      write (*,*) '   nelec  = ', nelec
      write (*,*) '   ncftot = ', ncftot
      write (*,*) '   nw     = ', nw
      write (*,*) '   nblock = ', nblock
      write (*,*)


      endif  !myid = 0

      CALL MPI_Bcast (nvectot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     CALL MPI_Bcast (NVEC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      CALL ALLOC(PNEVAL,nvectot,8)
      CALL ALLOC(PNEVEC,NCF*nvectot,8)
      CALL ALLOC(PNIVEC,nvectot,4)
      CALL ALLOC(PIATJP,nvectot,4)
      CALL ALLOC(PIASPA,nvectot,4)

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

      if (myid .eq. 0) then

      write (*,*) '  block     ncf     nev    2j+1  parity'

      DO jb = 1, nblock

         READ (68) nb, ncfblk, nevblk, iatjp, iaspa
                        WRITE (*,'(5I8)') nb, ncfblk, nevblk, iatjp, iaspa
         IF (jb .NE. nb) STOP 'jb .NE. nb'

         IF (nevblk .GT. 0) THEN

            READ (68) (ivec(nvecpat+i), i = 1, nevblk)
            DO i = nvecpat + 1, nvecpat + nevblk
               ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
               iatjpo(i) = iatjp
               iaspar(i) = iaspa
            ENDDO

            READ (68) eav, (eval(nvecpat + i ), i = 1, nevblk)

*           ...Construct the true energy by adding up the average
            DO i = 1, nevblk
               eval(nvecpat+i) = eval(nvecpat+i) + eav
            ENDDO
*           ...For overal (all blocks) average energy
            eavsum  = eavsum + eav*ncfblk
            neavsum = neavsum + ncfblk

            READ (68) ((evec( nvecsizpat + ncfpat+i + (j-1)*ncftot ),
     &                    i = 1, ncfblk), j = 1, nevblk)
         ENDIF

         nvecpat = nvecpat + nevblk
         ncfpat = ncfpat + ncfblk
         nvecsizpat = nvecsizpat + nevblk*ncftot

      ENDDO

!     READ(68) (IVEC(I),I=1,NVEC)
!     write(*,*) 'ivec',IVEC(1:NVEC)
!     READ(68) (IATJPO(I),IASPAR(I),I=1,NVEC)
!     write(*,*) 'iatjpo',IATJPO(1:NVEC), IASPAR(1:NVEC)
!     READ(68) EAV,(EVAL(I),I=1,NVEC)
!     write(*,*) 'eav',EAV, EVAL(1:NVEC)
!     READ(68) ((EVEC(I+(J-1)*NCF),I=1,NCF),J=1,NVEC)
!     write(*,*) 'evec',EVEC(1:NCF*NVEC)

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
*
*   Close the initial state mixing  file
*
      CLOSE(68)

      endif  !myid = 0

      CALL MPI_Bcast (IVEC,   nvectot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (IATJPO, nvectot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (IASPAR, nvectot,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (EAV, 1,   MPI_DOUBLE_PRECISION, 0,
     &                          MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (EVAL,   nvectot, MPI_DOUBLE_PRECISION, 0,
     &                          MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (EVEC,   NCF*nvectot, MPI_DOUBLE_PRECISION, 0,
     &                          MPI_COMM_WORLD, ierr)

      NVEC = nvectot
*
      RETURN
      ENd
