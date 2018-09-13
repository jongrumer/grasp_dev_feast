************************************************************************
*                                                                      *
      SUBROUTINE MANEIG (iatjpo, iaspar, NELMNT_a, imethod)
*   eval(), evec(), iatjpo, iaspar
*                                                                      *
*   This module  manages the  operation of the  eigensolvers and the   *
*   storage of the eigenpairs.  There are two principal branches:      *
*                                                                      *
*      (1) Matrix of order 1: the trivial case                         *
*      (2) Matrix of order greater than 1: there are two branches      *
*             (i) Matrices of order less than or equal to IOLPCK:      *
*                 eigenpairs are found using LAPACK SUBROUTINEs        *
*            (ii) Matrices of order greater than IOLPCK: eigenpairs    *
*                 are found using DVDSON; this involves up to three    *
*                 steps:                                               *
*                    (a) The matrix is analysed to determine its       *
*                        block structure (only irreducibe matrices     *
*                        are correctly treated by DVDSON)              *
*                    (b) Eigenpairs are extracted for each block       *
*                    (c) The appropriate eigenpairs are selected and   *
*                        stored                                        *
*                 Different methods of storage and different           *
*                 versions of the matrix-vector multiply are used      *
*                 depending upon the order and density of the matrix   *
*   
*   iatjpo - Output. (2j+1) of the dominant coefficients
*   iaspar - Output, Parity (1 or -1) of the dominant coefficients
*   nelmnt_a - Input, Total number of non-zero elements. Note that
*              the NELMNT in the common block is the number of 
*              non-zero elements of each node. nelmnt_a is the sum of
*              them.
*   imethod - Input. It tries to enforce a particular method by 
*            changing the criteria which are solely used to determine
*            a _suitable_ method.
*
*     value  meaning
*       1    LAPACK
*       2    Davidson, Dense-Memory
*       3    Davidson, Sparse-Memory
*       4    Davidson, Sparse-Disk
*     Other  No enforcement, Program (built-in parameters) decides
*
*   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, posfile,         *
*                        RALLOC.                                       *
*               [RCI92]: DNICMV, SPICMVmpi, SPODMV.               *
*               [DVDSON]: DVDSON.                                      *
*               [AUXBLAS]: DINIT/SINIT.                                *
*               [BLAS]: DCOPY/SCOPY, DSWAP/SSWAP.                      *
*               [LAPACK]: DSPEVX/SSPEVX.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
*   Modified Misha Saparov                                  Feb 1997   *
*            Charlotte F. Fischer                           May 1997   *
*                 Except for the disk version, all matrices have       *
*                 diagonals  shifted by EAV                            *
*  All arrays allocated here are de-allocated except pneval and pnevec
*  which will be de-allocated in matrix.
*   MPI Version  By Xinghong He           Last revision: 29 Jul 1998   *
*
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      EXTERNAL DNICMV,SPICMVmpi,SPODMV

      POINTER (PNEVAL,EVAL(1))
      COMMON/EIGVAL/EAV,PNEVAL

      POINTER (PNEVEC,EVEC(1))
      COMMON/EIGVEC/PNEVEC

!  nposition+1 is the current position of the .res file
!  It is set in matrix and used in maneig, spodmv

      COMMON/fposition/nposition

      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT

      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/ORB2/NCF,NW,PNTRIQ

      POINTER (PNIVEC,IVEC(1))
      COMMON/PRNT/NVEC,PNIVEC,NVECMX

      COMMON/WCHBLK/JBLOCK
      COMMON/WHERE/IMCDF
      COMMON/iounit/istdi,istdo,istde

      LOGICAL HIEND,LDISC,SPARSE
      CHARACTER*8 CNUM

      POINTER (PONTRW,W(1))
      POINTER (PONTRZ,Z(1))
      POINTER (PNWORK,WORK(1))
      POINTER (PIWORK,IWORK(1))
      POINTER (PIFAIL,IFAIL(1))
      POINTER (PNDIAG,DIAG(1))
      POINTER (PJWORK,JWORK(1))

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
      CHARACTER msg*80

!  ABSTOL is the absolute error tolerance for the eigenvalues
!  used by lapack routine

      PARAMETER (ABSTOL = 1.0D-10)
      INTEGER*8 NBRKEV,NDENSE
!-----------------------------------------------------------------------

      IF (myid .EQ. 0) PRINT *, 'Calling maneig...'

!  Default settings for controlling method to use
      IOLPCK = 1000
      NINCOR = 16777216             ! 16*1048576 words = 128 MB
      NBRKEV = (NCF+1)*(NCF+1) / 3
!  To enforce a particular method
      IF (imethod .EQ. 1) THEN       ! Lapack
         IOLPCK = NCF               ! So NCF <= IOLPCK 
      ELSEIF (imethod .EQ. 2) THEN   ! Davidson, Dense-Memory
         IOLPCK = 0
         NBRKEV = 0                 ! So nelmnt_a >= nbrkev
         !NINCOR = (ncf*(ncf+1))/2   ! So nstore < nincor
         NINCOR = (ncf*(ncf+1))   ! So nstore < nincor
      ELSEIF (imethod .EQ. 3) THEN   ! Davidson, Sparse-Memory
         IOLPCK = 0
         NBRKEV = (ncf*(ncf+1))/2 + 1
         !NINCOR = (ncf*(ncf+1))/2
         NINCOR = (ncf*(ncf+1))     !since nstore expression bad.
      ELSEIF (imethod .EQ. 4) THEN   ! Davidson, Sparse-Disk
         IOLPCK = 0
c        NBRKEV = (ncf*(ncf+1))/2 + 1
         NBRKEV = ncf
         NBRKEV = (NBRKEV*(NBRKEV+1))/2 + 1
         NINCOR = 0
      ENDIF
 
!  (nrows+1) is the number of records of the present block's .res file

      nrows = (ncf - myid - 1 + nprocs) / nprocs
      IF (ncf .LT. nprocs) nrows = ncf / (myid+1)
      !CALL posfile (1, imcdf, nrows+1)
      CALL posfile (0, imcdf, nposition)

      IF (NCF .EQ. 1) THEN
!-----------------------------------------------------------------------
!
! (1) - Trivial    ncf = 1
!
!-------------------------------------------------------
         IF (myid .EQ. 0) WRITE (24,*) 'Trivial eigenvalue problem.'

!   Matrix of order 1: the trivial case; we assume that the value
!   of EAV is available

         CALL ALLOC (PNEVAL, 1, 8)
         CALL ALLOC (PNEVEC, 1, 8)
         EVAL(1) = 0.D0
         EVEC(1) = 1.D0

         ! Still read through the .res file
CGG
CGG Gediminas NIST 2005.11.03
CGG         READ (imcdf)
         DO i = 1, nrows + 1
            READ (imcdf)
         ENDDO

      ELSE	!if-2
!-----------------------------------------------------------------------
!
! (2) - Non trivial
!
!-------------------------------------------------------

!   Matrix of order greater than 1; how many elements in a triangle?

         NDENSE = (NCF*(NCF+1))/2

         IF (NCF .LE. IOLPCK) THEN
!-----------------------------------------------------------------------
!
! (2.1) - LAPACK    Dense, Memory, 
!
!-------------------------------------------------------
            IF (myid .EQ. 0) THEN
               PRINT *, 'LAPACK routine DSPEVX'
     :                 //' selected for eigenvalue problem.'
               WRITE (24,*) 'LAPACK routine DSPEVX'
     :                 //' selected for eigenvalue problem.'
            ENDIF

!   Allocate storage for the dense representation of the matrix
!   and initialize emt

            CALL ALLOC (PNTEMT,NDENSE,8)
            CALL DINIT (NDENSE,0.0D 00,EMT,1)

!   Read the matrix into position from the disc file; it's already
!   been properly positioned.

            CALL ALLOC (PNWORK,NCF,8)
            CALL ALLOC (PNIROW,NCF,8)
            READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
            IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &              .OR. nprocsdum .NE. nprocs) 
     &          STOP 'maneig:1'

            DO 2 I = myid + 1, NCF, nprocs
               iofset = (i*(i-1))/2
               READ (IMCDF) NELC,ELSTO,(WORK(IR),IR = 1,NELC),
     :                                 (IROW(IR),IR = 1,NELC)
               ! In the row-mode of the lower triangle, 
               ! diagonal is the last one
               DO 1 IR = 1, NELC - 1
                  EMT(IOFSET+IROW(IR)) = WORK(IR)
    1          CONTINUE
               EMT(IOFSET+IROW(NELC)) = WORK(NELC)-EAV

    2       CONTINUE

            CALL DALLOC (PNWORK)
            CALL DALLOC (PNIROW)

! Let each node have a complete copy of EMT

            CALL gdsummpi (EMT, NDENSE)
!
!   Find the eigenpairs
!
!    ivec() - serial numbers of eigenstates of the current block
!    iccmin() - serial numbers of eigenstates of all blocks.
!    nvecmn - minimum serial number of the eigenstates of the block
!    nvecmx - maximum .............
!    nvex - clear from def: NVECMX-NVECMN+1
!
            NVECMN = NCF
            DO 3 I = 1,NVEC
               NVECMN = MIN (NVECMN,IVEC(I))
    3       CONTINUE
            NVEX = NVECMX-NVECMN+1
            CALL ALLOC (PONTRW,    NVEX,8)
            CALL ALLOC (PONTRZ,NCF*NVEX,8)
            CALL ALLOC (PNWORK,NCF*8   ,8)
            CALL ALLOC (PIWORK,NCF*5   ,4)
            CALL ALLOC (PIFAIL,    NVEX,4)
            CALL DSPEVX ('V','I','U',NCF,EMT,DUMMY,DUMMY,
     :                           NVECMN,NVECMX,ABSTOL,M,W,Z,NCF,
     :                           WORK,IWORK,IFAIL,INFO)
            IF (INFO .NE. 0) THEN
               CALL stopmpi ('maneig: Failure in DSPEVX [LAPACK]', myid)
            ENDIF
            CALL DALLOC (PNWORK)
            CALL DALLOC (PIWORK)
            CALL DALLOC (PIFAIL)
            CALL DALLOC (PNTEMT)

!   Store the eigenpairs in their proper positions EVAL() and EVEC()

            CALL ALLOC (PNEVAL,    NVEC,8)
            CALL ALLOC (PNEVEC,NCF*NVEC,8)
            DO 4 I = 1,NVEC
               LOC = IVEC(I)
               EVAL(I) = W(LOC-NVECMN+1)
               IOFSET = NCF*(I-1)
               LOC = NCF*(LOC-NVECMN)
               CALL DCOPY (NCF,Z(LOC+1),1,EVEC(IOFSET+1),1)
    4       CONTINUE
            CALL DALLOC (PONTRW)
            CALL DALLOC (PONTRZ)

         ELSE
!-----------------------------------------------------------------------
!
! (2.2) - DVDSON --- preparation work
!
!-------------------------------------------------------
            IF (myid .EQ. 0) WRITE (24,*) 'DVDSON routine selected'
     :                 //' for eigenvalue problem;'
!--------------------------------------------------------------
c           IF (myid .EQ. 0) print *,'zou',NELMNT_a, NBRKEV    
            IF (NELMNT_a .LT. NBRKEV) THEN
               SPARSE = .TRUE.
               NSTORE = NELMNT_a + NELMNT_a/2 + (NCF+1)/2
            ELSE
               SPARSE = .FALSE.
               NSTORE = NDENSE
            ENDIF

            CALL ALLOC (PNDIAG,NCF,8)
            CALL dinit (ncf, 0.d0, diag, 1)

            IF (NSTORE .GT. NINCOR) THEN
!-----------------------------------------------------------------------
!
! (2.2.1) - DVDSON --- Disk, load diagonal
!
!-------------------------------------------------------
               IF (myid .EQ. 0) WRITE (24,*) ' matrix stored on disc;'
!
!   Disk storage; necessarily sparse; one column of the matrix in
!   memory

               LDISC = .TRUE.
               SPARSE = .TRUE.
               IMV = 1

!   Load diagonal - Each node will have the same, complete copy
!   after this if block 

               READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
               IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &              .OR. nprocsdum .NE. nprocs) 
     &            STOP 'maneig:2'

               DO I = myid + 1, NCF, nprocs
                  READ (IMCDF) NELC,ELSTO,(dummy,IR=2,NELC), diatmp
     &                                 , (idummy, ir=1,nelc)
                  DIAG(I) = diatmp - EAV
               ENDDO

            ELSE
!-----------------------------------------------------------------------
!
! (2.2.2) - DVDSON --- Memory, load all
!
!-------------------------------------------------------
!
!   Core storage; load matrix into memory

               LDISC = .FALSE.
               IF (SPARSE) THEN
!-----------------------------------------------------------------------
!
! (2.2.2.1) - DVDSON --- Memory, load all, sparse
!
!-------------------------------------------------------
                  IF (myid .EQ. 0) WRITE (24,*) ' matrix stored in'
     :                       //' sparse representation in core;'

                  IMV = 2
                  CALL ALLOC (PNTEMT,NELMNT,8)
                  CALL ALLOC (PNIROW,NELMNT,4)
                  CALL ALLOC (PIENDC,NCF+1,4)
                  IOFSET = 0
                  IENDC(0) = 0
                  READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
                  IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &                  .OR. nprocsdum .NE. nprocs) 
     &               STOP 'maneig:3'
                  DO I = myid + 1, NCF, nprocs
                     READ (IMCDF) NELC,ELSTO,
     :                            (EMT(IR+IOFSET),IR = 1,NELC),
     :                           (IROW(IR+IOFSET),IR = 1,NELC)
                     EMT(NELC+IOFSET) = EMT(NELC+IOFSET) - EAV
                     DIAG(I) = EMT(NELC+IOFSET)
                     IOFSET = IOFSET + NELC
                     IENDC(I) = IOFSET
                  ENDDO
               ELSE
!-----------------------------------------------------------------------
!
! (2.2.2.2) - DVDSON --- Memory, load all, dense
!
!-------------------------------------------------------
                  IF (myid .EQ. 0) WRITE (24,*) ' matrix stored in'
     :                       //' full representation in core;'

                  IMV = 3

! Find NDENSE_L, the number of elements on the node (dense form)

                  NDENSE_L = 0
                  DO i = myid + 1, NCF, nprocs
                     NDENSE_L = NDENSE_L + i
                  ENDDO

                  CALL ALLOC (PNTEMT,NDENSE_L,8)
                  CALL DINIT (NDENSE_L,0.0D 00,EMT,1)
                  CALL ALLOC (PNWORK,NCF,8)
                  CALL ALLOC (PNIROW,NCF,4)

                  READ (IMCDF) ncfdum, iccutdum, myiddum, nprocsdum
                  IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &                  .OR. nprocsdum .NE. nprocs) 
     &               STOP 'maneig:4'

                  IOFSET = 0
                  DO 8 I = myid + 1, NCF, nprocs
                     READ (IMCDF) NELC,ELSTO,
     :                            (WORK(IR),IR = 1,NELC),
     :                            (IROW(IR),IR = 1,NELC)
                     work(nelc) = work(nelc) - eav
                     DIAG(I) = WORK(NELC)
                     DO IR = 1, NELC
                        EMT(IOFSET+IROW(IR)) = WORK(IR)
                     ENDDO
                     iofset = iofset + i
    8             CONTINUE
                  CALL DALLOC (PNWORK)
                  CALL DALLOC (PNIROW)

               ENDIF
!               ...Memory mode - sparse or dense
!-----------------------------------------------------------------------
!  (2.2.2.3e)          *** E n d   o f   D V D S O N   m e m o r y
!-----------------------------------------------------------------------
            ENDIF
!            ...Disk or Memory
!-----------------------------------------------------------------------
!  (2.2.3e)      *** E n d   o f   D V D S O N
!-----------------------------------------------------------------------

! Make diagonals global, no matter it is disk or memory mode

            CALL gdsummpi (DIAG, NCF)

!   Allocate storage for workspace; see the header of DVDSON for
!   the expression below; the value of LIM can be reduced to NVECMX
!   plus a smaller number if storage is severely constrained

            LIM = MIN (NCF,NVECMX+20)
c           LWORK = 2*NCF*LIM+LIM*LIM+(NVECMX+10)*LIM+NVECMX
            LWORK = 2*NCF*LIM+LIM*LIM*2+11*LIM+NVECMX
            CALL ALLOC (PNWORK,LWORK,8)
            work(1:lwork) = 0.0d0 
            LIWORK = 6*LIM+NVECMX
            CALL ALLOC (PIWORK,LIWORK,4)
!*changed by Misha 02/12/97
            CRITE = 1.0D-17
            CRITC = 1.0D-09
            CRITR = 1.0D-09
            ORTHO = max(1D-8,CRITR)
! end of changes

            MAXITR = MAX (NVECMX*100,200)
            CALL ALLOC (PJWORK,LIM,4)

            CALL ALLOC (PNEVAL,    NVECMX,8)
            CALL ALLOC (PNEVEC,NCF*NVECMX,8)

            DMUNGO = 10.D99
            CALL DINIT (NVECMX,DMUNGO,EVAL,1)

!   Compute the eigenpairs in each block

            NVEX = nvecmx
            IF (LDISC) THEN
               MBLOCK = NVEX
            ELSE
               MBLOCK = 1
            ENDIF
            NEND = NCF*NVEX

            ILOW = 1
            IHIGH = NVEX
            NIV = NVEX
!***********************************************************************
!
!   Call Davidson eigensolver
!
            IF (IMV .EQ. 1) THEN
!******************** sparse and matrix on disk **********************
                !print* , ' imv = 1'
               IF (myid .EQ. 0) print *, ' Sparse - Disk, iniestsd'
      CALL posfile (0, imcdf, nposition)  ! was within iniestsd before
	   CALL INIESTSD (1000, ncf, myid, nprocs, NIV, work, IMCDF, EAV)
        if (ncf.gt.1000) then
               CALL GDVD (SPODMV,NCF,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)
         end if
            ELSEIF (IMV .EQ. 2) THEN
                !print*, ' imv = 2'
!******************** sparse and matrix in memory ********************
               IF (myid .EQ. 0) print *, ' Sparse - Memory, iniestmpi'
	       CALL iniestmpi (1000, NCF,NIV,WORK,EMT,IENDC,IROW)
        if(ncf.gt.1000) then
               CALL GDVD (SPICMVmpi,NCF,LIM,DIAG,ILOW,IHIGH,
     :            JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :            WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :            NMV,IERR)
        end if

               CALL dalloc (pntemt)
               CALL dalloc (pnirow)
               CALL dalloc (piendc)

            ELSEIF (IMV .EQ. 3) THEN
              !print*, ' imv = 3'
!*************************** dense and in memory **********************
               IF (myid .EQ. 0) print *, ' Dense - Memory, iniestdm'
	       CALL INIESTDM (1000,NCF,NIV,WORK,EMT)
         if (ncf.gt.1000) then
               CALL GDVD (DNICMV,NCF,LIM,DIAG,ILOW,IHIGH,
     :              JWORK,NIV,MBLOCK,CRITE,CRITC, CRITR,ORTHO,MAXITR,
     :              WORK,LWORK,IWORK,LIWORK,HIEND,NLOOPS,
     :              NMV,IERR)
         end if
               CALL dalloc (pntemt)
            ENDIF
!***********************************************************************
            CALL DALLOC (PNDIAG)
            CALL DALLOC (PIWORK)
            CALL DALLOC (PJWORK)

            IF (myid .EQ. 0) THEN
               WRITE (24,*) ' ', nloops, ' iterations;'
               WRITE (24,*) ' ', nmv,' matrix-vector multiplies.'
            ENDIF

            IF (IERR .NE. 0) THEN
               WRITE (istde,*) 'MANEIG: Returned from DVDSON with'
               WRITE (istde,*) ' IERR = ',IERR,'.'
               CALL stopmpi ('maneig: DVDSON wrong', myid)
            ENDIF

!   Put the eigenpairs in order, overwriting as necessary

            CALL DCOPY (NVEX    , WORK(NEND+1), 1, EVAL, 1)
            CALL DCOPY (NCF*NVEX, WORK(1)     , 1, EVEC, 1)
            CALL DALLOC (PNWORK)

!   Rearrange and reallocate storage for the eigenpairs
!   as necessary

	    IF (NVEC .LT. NVECMX) THEN
               CALL ALLOC (PIWORK,NVECMX,4)
               DO I = 1,NVECMX
                  IWORK(I) = I
               ENDDO
               DO I = 1,NVEC
                  IOFSET = IVEC(I)
                  LOC = IWORK(I)
                  IF (IOFSET .NE. LOC) THEN
                     CALL DSWAP (1,EVAL(IOFSET),1,
     :                             EVAL(I     ),1)
                     IWORK(I) = IWORK(IOFSET)
                     IWORK(IOFSET) = LOC
                     IOFSET = NCF*(IOFSET-1)
                     LOC = NCF*(I-1)
                     CALL DSWAP (NCF,EVEC(IOFSET+1),1,
     :                               EVEC(LOC   +1),1)
                  ENDIF
               ENDDO
               CALL DALLOC (PIWORK)
               CALL RALLOC (PNEVAL,    NVECMX,    NVEC,8)
               CALL RALLOC (PNEVEC,NCF*NVECMX,NCF*NVEC,8)
	    ENDIF

         ENDIF
! (2.3e)              *** E N D   O F    N O N - T R I V I A L   C A S E

      ENDIF
! (3e)               *** E N D   O F    A L L

!--------------------------------------------------------------------
! Only the following quantities are needed after this routine is
! finished:
!   eval(), evec(), iatjpo, iaspar
!--------------------------------------------------------------------
!
!   Clean up eigenvectors; determine their J/P values
!
      DO 23 J = 1, NVEC

!         Find the dominant component of each eigenvector

         IOFSET = (J-1)*NCF

         AMAX = 0.d0
         DO  I = 1, NCF
            WA = ABS (EVEC(I+IOFSET))
            IF (WA .GT. AMAX) THEN
               AMAX = WA
               IA = I
            ENDIF
         ENDDO

!          Find the angular momentum and parity of the dominant component

         iatjpo = ITJPO (IA)
         iaspar = ISPAR (IA)

!          Change sign of eigenvactor if dominant component is negative

         IF (EVEC(IA+IOFSET) .LT. 0.d0) THEN
           DO I = 1, NCF
              EVEC(I+IOFSET) = -EVEC(I+IOFSET)
           ENDDO
         ENDIF
   23 CONTINUE

      RETURN
      END
