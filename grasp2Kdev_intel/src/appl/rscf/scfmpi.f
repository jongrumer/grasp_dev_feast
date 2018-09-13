************************************************************************
*                                                                      *
      SUBROUTINE scfmpi (EOL, rwffile2)
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL   EOL
      CHARACTER rwffile2*(*)  !RWF output file, passed to orbout
*                                                                      *
*   This  subroutine  performs  the SCF iterations. The procedure is   *
*   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
*   (1986) 290.                                                        *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
*               [RSCF92]: improvmpi, matrixmpi, MAXARR, newcompi,      *
*                         ORBOUT, ORTHSC, setlagMpi.                   *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
*   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
*                                                                      *
************************************************************************
*
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PCCMIN,ICCMINDUM),
     :        (PNEVAL,EVAL(*)),
     :        (PNEVEC,EVEC(*)),
     :        (PNTEMT,EMTDUM),
     :        (PNIECC,IECCDUM), (PNTECV,ECVDUM),
     :        (PNTRIQ,IQADUM), 
     :        (PNTRWT,RWTDUM),
     :        (PNTNDA,NDADUM), (PNTNXA,NXADUM),
     :        (PNTNYA,NYADUM), (PNTRDA,RDADUM),
     :        (PNTRXA,RXADUM), (PNTRYA,RYADUM),
     :        (PIATJP,IATJPO(*)), (PIASPA,IASPAR(*)),
     :        (PNTJQS,JQSDUM), (PNJCUP,JCUPDUM),
     :        (PWEIGH,WEIGHTdum)

      LOGICAL CONVG, LDBPR, LFIX, orthst, lsort
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /LAGR/PNTECV,PNIECC,NEC
     :      /MCPA/KMAXF
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORBA/IORDER(NNNW)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /SYMA/PIATJP,PIASPA
     :      /STAT/PNTJQS,PNJCUP
     :      /ORTHCT/orthst

      COMMON/DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde
      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (pncfpast, ncfpast(*))
      POINTER (pncminpast, ncminpast(*))
      POINTER (pnevecpast, nevecpast(*))
      COMMON/pos/pncfpast,pncminpast,pnevecpast,ncftot,nvecsiz

      POINTER (peavblk, eavblk(*))
      COMMON/peav/peavblk

      POINTER (pidxblk, idxblk(*))   ! idx(i= 1,ncmin) is the block where 
      COMMON/blkidx/pidxblk          ! the i_th eigenvalue comes from

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      ncftot = ncf
      !IF (myid .EQ. 0) PRINT *, '===SCF==='

*=======================================================================
*   Determine Orthonomalization order --- lsort
*=======================================================================

      IF (ndef .EQ. 0) THEN
         lsort = .FALSE.
      ELSE
         IF (myid .EQ. 0) THEN
  123       WRITE (istde,*) 'Orthonomalization order? '
            WRITE (istde,*) '     1 -- Update order'
            WRITE (istde,*) '     2 -- Self consistency connected'
            READ (istdi,*) j
            IF (j .EQ. 1) THEN
               lsort = .FALSE.
            ELSE IF (j .EQ. 2) THEN
               lsort = .TRUE.
            ELSE
               WRITE (istde,*) 'Input is wrong, redo...'
               GOTO 123
            ENDIF
         ENDIF
         CALL MPI_Bcast (lsort, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      ENDIF

*=======================================================================
*   Deallocate storage that will no longer be used
*=======================================================================

      CALL DALLOC (PNTJQS)

*=======================================================================
*   Allocate and fill in auxiliary arrays
*=======================================================================

      CALL alloc (pncfpast, nblock, 4)
      CALL alloc (pncminpast, nblock, 4)
      CALL alloc (pnevecpast, nblock, 4)
      CALL alloc (peavblk, nblock, 8)

      ncfpast(1)   = 0
      ncminpast(1) = 0
      nevecpast(1) = 0
      DO i = 2, nblock
         ncfpast(i)   =   ncfpast(i-1) + ncfblk(i-1)
         ncminpast(i) = ncminpast(i-1) + nevblk(i-1)
         nevecpast(i) = nevecpast(i-1) + nevblk(i-1) * ncfblk(i-1)
      ENDDO

      !*** Size of the eigenvector array for all blocks
      nvecsiz  = nevecpast(nblock) + nevblk(nblock) * ncfblk(nblock)

      IF (EOL) THEN
         CALL alloc (pneval, ncmin, 8)
         CALL alloc (pnevec, nvecsiz, 8)
         CALL alloc (piatjp, ncmin, 4)
         CALL alloc (piaspa, ncmin, 4)
      ENDIF

*=======================================================================
*   
*=======================================================================
      NDDIM = 0
      NXDIM = 0
      NYDIM = 0
*     Setlagmpi should only be called after the call to newcompi which
*     defines the occupation of orbitals.  CFF (June, 2004).
*      CALL setlagmpi (EOL); 
c     pRINT*,'eol=',eol,'?',myid
      !call flush(6)
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr); 
c     pRINT*,'Barrierscf1','?',myid
      !call flush(6)
*   For (E)OL calculations, determine the level energies and
*   mixing coefficients

      IF (EOL) THEN
         CALL matrixmpi
         CALL newcompi(WTAEV)
      ENDIF
      WTAEV0 = 0.0

      DO 4 NIT = 1, NSCF
         IF (myid .EQ. 0) WRITE (*,301) NIT

*   For all pairs constrained through a Lagrange multiplier, compute
*   the Lagrange multiplier

         CALL setlagmpi (EOL)
c     pRINT*,'eol=',eol,'?',myid
      !call flush(6)
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr); 
c     pRINT*,'Barrierscf1','?',myid
      !call flush(6)

*   Improve all orbitals in turn

         DAMPMX = 0.0
         IF (myid .EQ. 0) WRITE (*,302) 
         DO J = 1, NW
            JSEQ = IORDER(J)
            IF (.NOT. LFIX(JSEQ)) THEN
               CALL improvmpi (EOL, JSEQ, lsort,DAMPMX)
            ENDIF
         ENDDO
*
*   For KOUNT = 1 to NSIC: find the least self-consistent orbital;
*   improve it
*
!        write(*,*) 'nsic=',nsic
         DO KOUNT = 1, NSIC
            CALL MAXARR (K)
            IF (K .EQ. 0) THEN
               CONVG = .TRUE.
               GOTO 3
            ELSE
               IF (SCNSTY(K) .LE. ACCY) THEN
                  CONVG = .TRUE.
                  GOTO 3
               ENDIF
            ENDIF
            CALL improvmpi (EOL, K, lsort,DAMPMX)
         ENDDO

         CALL MAXARR (K)

         IF (K .EQ. 0) THEN
            CONVG = .TRUE.
         ELSE
            IF (SCNSTY(K) .LE. ACCY) THEN
               CONVG = .TRUE.
            ELSE
               CONVG = .FALSE.
            ENDIF
         ENDIF

    3    IF (LDBPR(24) .AND. myid .EQ. 0) CALL PRWF (0)

*   Perform Gram-Schmidt process
*   For OL calculation, orthst is true and orbitals are orthonormalized
*   in subroutine improv. For AL calculation, orthst is false.
         IF (.NOT. orthst ) CALL ORTHSC

*   Write the subshell radial wavefunctions to the .rwf file

         IF (myid .EQ. 0) CALL ORBOUT (rwffile2)

         IF (EOL) THEN
            CALL matrixmpi
            CALL newcompi(WTAEV)
         ENDIF
!        print *, "test print: ", WTAEV,WTAEV0,DAMPMX, myid
         IF(ABS(WTAEV-WTAEV0).LT.1.0D-8.and.
     &                 DAMPMX.LT.1.0D-2) CONVG=.true.
         WTAEV0=WTAEV                                                              
         IF (CONVG) THEN
            IF (LDBPR(25) .AND. (.NOT. LDBPR(24)) .AND. myid .EQ. 0)
     &          CALL PRWF (0)
            !IF (EOL) CALL matrixmpi
            GOTO 5
         ENDIF

    4 CONTINUE

      IF (myid .EQ. 0) 
     &   WRITE (istde,*) ' Maximum iterations in SCF Exceeded.'

    5 DO I = 31, 32 + KMAXF
         CLOSE (I)   ! The MCP coefficient files
      ENDDO

      IF (myid .EQ. 0) THEN
         !CLOSE (23)     ! The .rwf file
         CLOSE (25)     ! The .mix file
      ENDIF
*
*   Complete the summary - moved from rscf92 for easier alloc/dalloc
*
      IF (myid .EQ. 0) CALL ENDSUM
*
*   Deallocate storage
*
      call dalloc (pntrwt)		!Either getold or getald

      IF (NEC .GT. 0) THEN
         CALL dalloc (pniecc)
         CALL dalloc (pntecv)
         CALL dalloc (pntriq)
      ENDIF

      IF (NDDIM .GT. 0) THEN
         CALL dalloc (pntrda)
         CALL dalloc (pntnda)
      ENDIF

      IF (NXDIM .GT. 0) THEN
         CALL dalloc (pntrxa)
         CALL dalloc (pntnxa)
      ENDIF

      IF (NYDIM .GT. 0) THEN
         CALL dalloc (pntrya)
         CALL dalloc (pntnya)
      ENDIF

      IF (EOL) THEN
         CALL dalloc (pneval)
         CALL dalloc (pnevec)
         CALL dalloc (piatjp)
         CALL dalloc (piaspa)
         CALL dalloc (pncmaxblk)	! getold.f
         CALL dalloc (peavblk)   ! getold.f
         CALL dalloc (pidxblk)	! Allocated in getold.f
         CALL dalloc (pccmin)		! Allocated in items.f<-getold.f
      ENDIF
*
      CALL dalloc (pncfpast)
      CALL dalloc (pncminpast)
      CALL dalloc (pnevecpast)

  301 FORMAT (/' Iteration number ',1I3
     :        /' --------------------')
  302 FORMAT (41X,'Self-            Damping'
     :        /'Subshell    Energy    Method   P0    '
     :        ,'consistency  Norm-1  factor  JP'
     :        ,' MTP INV NNP'/)

      RETURN
      END
