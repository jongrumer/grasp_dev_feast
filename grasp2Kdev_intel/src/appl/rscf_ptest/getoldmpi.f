************************************************************************
*                                                                      *
      SUBROUTINE GETOLDmpi (idblk)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*8 idblk(*)
*                                                                      *
*   Interactively determines the data governing OL problem.            *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, GETRSL, GETYN, RALLOC.         *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*   Modified by Xinghong He               Last revision: 10 Jun 1998   *
*                                                                      *
************************************************************************
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

      POINTER (PNTRWT,RWTDUM)
      INTEGER*4 IQADUM
      POINTER (PNTRIQ,IQADUM)
      LOGICAL GETYN,LFIX,NOINVT,ORTHST,YES
      CHARACTER*256 RECORD
*
      DIMENSION indx(NNNW)
*
      POINTER (PCDAMP,CDAMP(*))
      POINTER (PWEIGH,WEIGHT(*))

      POINTER (PCCMIN,ICCMIN(*))
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX
*
!D     CHARACTER*2 NH   ! Test
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF5/PNTRWT,PWEIGH
     :      /DEFAULT/NDEF
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /INVT/NOINVT(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORBA/IORDER(NNNW)
     :      /ORTHCT/ORTHST
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
      COMMON/iounit/istdi,istdo,istde
      LOGICAL lcorre
      COMMON/corre/lcorre(NNNW)
      SAVE /corre/

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (pidxblk, idxblk(*)) ! idx(i= 1,ncmin) is the block where 
      COMMON/blkidx/pidxblk        ! the i_th eigenvalue comes from

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr

!-----------------------------------------------------------------------

      !IF (myid .EQ. 0) WRITE (istde,*) '(E)OL type calculation;'

! lodstate fills
!    nevblk(), ncmaxblk()
!    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)

      CALL alloc (pncmaxblk, nblock, 4)
      CALL alloc (pnevblk, nblock, 4)

      IF (myid .EQ. 0) THEN
         CALL lodstate (nblock, ncfblk(1), idblk, nevblk, ncmaxblk)
      ENDIF

      CALL MPI_Bcast (ncmin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (nevblk(1), nblock, MPI_INTEGER, 0, 
     &                                 MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ncmaxblk(1), nblock, MPI_INTEGER, 0, 
     &                                 MPI_COMM_WORLD, ierr)

      ! actually, ncmin will always be positive
      IF (ncmin .NE. 0) THEN
         IF (myid .NE. 0) THEN   ! lodstate allocated it from node-0
            CALL alloc (pccmin, ncmin, 4)
         ENDIF

         CALL MPI_Bcast (iccmin(1), ncmin, MPI_INTEGER, 0, 
     &                                 MPI_COMM_WORLD, ierr)
      ENDIF
*
*   Allocate the storage for and set the weights
*
      CALL ALLOC (PWEIGH, NCMIN, 8)

      IF (myid .EQ. 0) THEN
         CALL getoldwt ((ndef), (ncmin), weight)
      ENDIF
      CALL MPI_Bcast (weight(1), ncmin, MPI_DOUBLE_PRECISION, 0, 
     &                                 MPI_COMM_WORLD, ierr)
*
*   Eigenvector damping
*
      CALL ALLOC (PCDAMP, NCMIN, 8)
*
      DO I = 1, NCMIN
         CDAMP(I) = 0.D0
      ENDDO
*
*   Print the list of all subshells
*
      IF (myid .EQ. 0) THEN
         WRITE (istde,*) ' Radial functions'
         CALL PRTRSL
      ENDIF
*
*   Determine which orbitals are to be varied, which are fixed.
*   Quantities determined here:
*     nfix, lfix(1:nw), iorder(1:nw), scnsty(1:nw)
*   Instead of broadcasting these quantities, we broadcast 
*   the intermediate result from GETRSL (see below)
*
      DO I = 1, NW
         LFIX(I) = .TRUE.
      ENDDO

            IF (myid .EQ. 0) THEN
      WRITE (istde,*) 'Enter orbitals to be varied (Updating order)'
      CALL GETRSL (indx, NSUBS)
            ENDIF
      CALL MPI_Bcast (nsubs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      IF (nsubs .GT. 0)
     &   CALL MPI_Bcast (indx(1), nsubs, MPI_INTEGER, 0, 
     &              MPI_COMM_WORLD, ierr)
      
      DO I = 1, NSUBS
         LFIX(indx(I)) = .FALSE.
!XHH      give a big value, rather than zero to scnsty()
         scnsty(indx(I)) = 1.D20
      ENDDO
      NFIX = NW - nsubs
      IF (NFIX .EQ. NW) THEN
         IF (myid .EQ. 0)
     &   WRITE (istde,*) 'All subshell radial wavefunctions are fixed;'
     &,                 ' perform CI calculations with RCI92.'
      ENDIF

!   Determine orbital updating order

      NORDER = 0
      DO  I = 1, NW
         IORDER(I) = I
         IF (.NOT. LFIX(I)) THEN
            NORDER = NORDER + 1
            IORDER(I) = indx(NORDER)
         ENDIF
!D        mp = IORDER(I)
!D        WRITE (istde,*) NP(mp),NH(mp),NAK(mp), LFIX(mp)  ! Test
      ENDDO
!----------------------------------------------------------------
!
! Fixed orbitals first
!
!      DO  I = 1, NFIX
!         IF (LFIX(I)) THEN
!            IORDER(I) = I
!         ELSE
!            J = I + 1
!            DO WHILE (.NOT. LFIX(J))
!               J = J + 1
!            ENDDO
!            IORDER(I) = J
!         ENDIF
!D        mp = IORDER(I)
!D        WRITE (istde,*) NP(mp),NH(mp),NAK(mp), LFIX(mp)  ! Test
!      ENDDO
!
! Varying orbitals last
!
!      NORDER = 0
!      DO I = NFIX + 1, NW
!         NORDER = NORDER+1
!         IORDER(I) = indx(NORDER)
!D        mp = IORDER(I)
!D        WRITE (istde,*) NP(mp),NH(mp),NAK(mp), LFIX(mp)  ! Test
!      ENDDO
!D     pause
!****************************************************************
!
!XHH added a array to store the index of the correlation functions
!
      DO i = 1, nw
         lcorre(i) = .TRUE.
      ENDDO

            IF (myid .EQ. 0) THEN
      WRITE (istde,*) 'Which of these are spectroscopic orbitals?'
      CALL GETRSL (indx, NSUBS)
            ENDIF
      CALL MPI_Bcast (nsubs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      IF (NSUBS .GT. 0) THEN
         CALL MPI_Bcast (indx(1), nsubs, MPI_INTEGER, 0, 
     &                                MPI_COMM_WORLD, ierr)
         DO I = 1, NSUBS
            LOC = indx(I)
            IF (.NOT. LFIX(LOC)) THEN
               METHOD(LOC) = 1
               NOINVT(LOC) = .FALSE.
               ODAMP(LOC)  = 0.D0
               lcorre(Loc) = .FALSE.
            ENDIF
        ENDDO
      ENDIF

! Set NSIC. It will be non-zero if all orbitals to be varied are 
! spectroscopic orbitals

      NSIC = 4 + (NW - NFIX) / 4
      DO i = 1, nw
         IF ((.NOT. LFIX(i)) .AND. lcorre(i)) THEN
            NSIC = 0
            EXIT
         ENDIF
      ENDDO
*
      NSCF = 24
      NSOLV = 3
      ORTHST = .TRUE.
*
*   Make the allocation for the auxiliary vector required
*   by SUBROUTINE NEWCO
*
      CALL alloc (pntrwt, ncmin, 8)
*
*   Place the block numbers of the all ncmin eigenstate(wanted)
*   in array idxblk
*
      CALL alloc (pidxblk, ncmin, 4)
      noffset = 0
      DO jblock = 1, nblock
         DO j = 1, nevblk(jblock)
            idxblk (j + noffset) = jblock
         ENDDO
         noffset = noffset + nevblk(jblock)
      ENDDO
      IF (noffset .NE. ncmin) STOP 'getold: ncmin trouble'

      RETURN
      END
