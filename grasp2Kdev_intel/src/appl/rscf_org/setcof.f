************************************************************************
*                                                                      *
      SUBROUTINE SETCOF (EOL,J)
*                                                                      *
*   This  subroutine  sets  up the coefficients and orbital pointers   *
*   for the direct and exchange potentials for orbital  J .  It also   *
*   sets  up  the  coefficients  and  pointers for the inhomogeneous   *
*   terms arising from off-diagonal I (a,b) integrals.                 *
*                                                                      *
*   Call(s) to: [LIB92]: alloc, dalloc.                                *
*               [RSCF92]: alcsca, dsubrs.                              *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 21 Dec 1992   *
*   Modified by Xinghong He                 Last update: 21 Dec 1997   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)

*      integer*8 nelmnt

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PNTEMT,EMTDUMMY)
      POINTER (PIENDC,IENDCDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL EOL
*** Locals
      POINTER (PCOEFF,COEFF(1))
      POINTER (PICLMN,ICLMN(1))
      POINTER (PINDEX,indx(1))
*
      POINTER (PNIROW,IROW(1))
      POINTER (PNTRDA,DA(1))
      POINTER (PNTRXA,XA(1))
      POINTER (PNTRYA,YA(1))
      POINTER (PNTNDA,NDA(1))
      POINTER (PNTNXA,NXA(1))
      POINTER (PNTNYA,NYA(1))
*
      DIMENSION INDEXS(4)

      CHARACTER*2 NH    ! For printing
      COMMON/ORB10/NH(NNNW)

      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /MCPA/KMAXF
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /SCF1/UCF(NNNW)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
*
      PARAMETER (EPS = 1.0D-10)
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(1))
      POINTER (pncmaxblk, ncmaxblk(1))
      COMMON/hblock2/pnevblk, pncmaxblk

      CHARACTER mcplab*3, idstring*3, msg*128
      CHARACTER*(*), PARAMETER:: myname = 'SETCOF'

      POINTER (pncfpast, ncfpast(1))
      POINTER (pncminpast, ncminpast(1))
      POINTER (pnevecpast, nevecpast(1))
      COMMON/pos/pncfpast,pncminpast,pnevecpast,ncftot,nvecsiz

      COMMON/iounit/istdi,istdo,istde
      COMMON /mpi/ myid, nprocs, ierr

*=======================================================================
*   Initializations
*=======================================================================

      NDIM = 1
      CALL alloc (PCOEFF,NDIM,8)
      CALL alloc (PICLMN,NDIM,4)
      CALL alloc (PINDEX,NDIM,4)

      NDCOF = 0
      NXCOF = 0
      NYCOF = 0

      NAKJ = NAK(J)
      NKJJ = NKJ(J)
      UCFJ = UCF(J)

*=======================================================================
*   Generate YA coefficients that do not require MCP output list.
*   Computation distributed and then collected.
*=======================================================================

      ILABEL = 0
      NWTERM = KEY * (KEY + 1)
      DO 3 IB = 1, NW
         ILABEL = ILABEL + NWTERM
         IF (IB .EQ. J) THEN
            KMAX = NKJJ - 1
         ELSE
            KMAX = 0
         ENDIF
         DO 2 K = 0, KMAX, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            SUMR = 0.D0
            DO nb = 1, nblock
            DO IR = myid + 1, NCFblk(nb), nprocs
                  SUMR = SUMR + dsubrs (EOL, IR, IR, nb) *
     &                     FCO (K, IR + ncfpast(nb), J, IB)
            ENDDO
            ENDDO
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            IF (IB .EQ. J) THEN
               YKAB = 2.D0 * SUMR / UCFJ
            ELSE
               YKAB = SUMR / UCFJ
            ENDIF

           !*** The following IF has to be removed since YKAB 
           !*** is incomplete in multiprocessor case
            !IF (ABS (YKAB) .GT. EPS) THEN
               NYCOF = NYCOF + 1
               IF (NYCOF .GT. NYDIM) THEN
                  IF (NYDIM .GT. 0) THEN
                     CALL alcsca (PNTNYA, PNTRYA, NYDIM, 2)
                  ELSE
                     CALL alcsca (PNTNYA, PNTRYA, NYDIM, 1)
                  ENDIF
               ENDIF
               YA(NYCOF) = YKAB
               NYA(NYCOF) = K + ILABEL
            !ENDIF

    2    CONTINUE
    3 CONTINUE

*=======================================================================
*   Generate XA coefficients that do not require MCP output list
*   Computation distributed and then collected.
*=======================================================================

      ILABEL = KEY * J
      NWTERM = KEY * KEY * (KEY + 1)
      DO 6 IB = 1, NW
         ILABEL = ILABEL + NWTERM
         IF (IB .EQ. J) CYCLE
         NKJIB = NKJ(IB)
         IF (NAKJ * NAK(IB) .GT. 0) THEN
            KMIN = ABS ( (NKJJ - NKJIB) / 2 )
         ELSE
            KMIN = ABS ( (NKJJ - NKJIB) / 2 ) + 1
         ENDIF
         KMAX = (NKJJ + NKJIB) / 2

         DO 5 K = KMIN, KMAX, 2

            !<<< mpi distribute calculation <<<<<<<<<<<<<<<<<<<<<<
            SUMR = 0.D0
            DO nb = 1, nblock
            DO IR = myid + 1, NCFblk(nb), nprocs
                  SUMR = SUMR + dsubrs (EOL, IR, IR, nb) *
     &                     GCO (K, IR + ncfpast(nb), J, IB)
            ENDDO
            ENDDO
            !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            XKAB = SUMR / UCFJ

           !*** The following IF has to be removed since YKAB 
           !*** is incomplete in multiprocessor case
            !IF (ABS (XKAB) .GT. EPS) THEN
               NXCOF = NXCOF + 1
               IF (NXCOF .GT. NXDIM) THEN
                  IF (NXDIM .GT. 0) THEN
                     CALL alcsca (PNTNXA, PNTRXA, NXDIM, 2)
                  ELSE
                     CALL alcsca (PNTNXA, PNTRXA, NXDIM, 1)
                  ENDIF
               ENDIF
               XA(NXCOF) = XKAB
               NXA(NXCOF) = K + ILABEL
            !ENDIF
    5    CONTINUE
    6 CONTINUE

*=======================================================================
*   Subroutine setham called from matrix had gone through the mcp 
*   files once. setcof will do it again. (Go to setmcp to see 
*   contents of these records.
*=======================================================================

      DO nfile = 30, 32 + kmaxf
         REWIND (nfile)
         DO i = 1, 3
            READ (nfile)
            IF (nfile .EQ. 30) READ (nfile)
         ENDDO
      ENDDO

*=======================================================================
*   Generate DA coefficients; these arise from the one-electron
*   integrals
*=======================================================================

      DO jblock = 1, nblock

         !*** Read in IROW from file mcp.30 ***
         READ (30) mcplab, jblockt, ncf
         IF (jblockt .NE. jblock) THEN
            WRITE (istde,*) myname, '1: jblockt .NE. jblock'
            STOP
         ENDIF
         READ (30) NELMNT
         CALL alloc (PNIROW, NELMNT, 4)
         READ (30) (IENDCdum, I = myid + 1, NCF, nprocs),
     &             (IROW(I), I = 1, NELMNT)

         !*** Block info file mcp.31 ***
         READ (31) mcplab, jblockt, ncf, ncoeff
         IF (jblockt .NE. jblock) THEN
            WRITE (istde,*) myname, '2: jblockt .NE. jblock'
            STOP
         ENDIF

         !*** Loop over labels having non-zero coefficients 
         !*** it exits when no more labels for the block

         DO
            READ (31, IOSTAT = IOS) lab, ncontr

            ! 0, 0 marks the end of a block. This is the normal exit

            IF (lab .EQ. 0 .AND. ncontr .EQ. 0) THEN
               CALL dalloc (PNIROW)
               EXIT     ! Actually to next block
            ENDIF

            !*** Decode the labels of I(ab) ***
            IA = MOD (LAB, KEY)
            IB = LAB / KEY
            
            ! At least one orbital should be J in order to have 
            ! non-zero value; otherwise, goto next label.

            IF ((IA .NE. J) .AND. (IB .NE. J)) THEN
               READ (31)   ! No contributions from this integral; skip
               CYCLE       ! to next label
            ENDIF

            IF (NCONTR .GT. NDIM) THEN
               CALL dalloc (PCOEFF)
               CALL dalloc (PICLMN)
               CALL dalloc (PINDEX)
               NDIM = NCONTR
               CALL alloc (PCOEFF, NDIM, 8)
               CALL alloc (PICLMN, NDIM, 4)
               CALL alloc (PINDEX, NDIM, 4)
            ENDIF

            ! Read the column index, the sparse matrix index, and the
            ! coefficient for all contributions from this integral

            READ (31) (ICLMN(I), indx(I), COEFF(I), I = 1, NCONTR)

            ! Add up all the contributions from this integral;
            ! off-diagonal contributions have double the weight

            SUM = 0.D0
            DO I = 1, NCONTR
               IR = IROW(indx(I))
               IC = ICLMN(I)

               CONTR = dsubrs (EOL, IR, IC, jblock) * COEFF(I)

               IF (IR .NE. IC) CONTR = CONTR + CONTR
               SUM = SUM + CONTR
            ENDDO
            SUM = 0.5D0 * SUM / UCFJ

            ! Put coefficients in the list. Since there is always
            ! (almost) some repeatence from different blocks, a check 
            ! and merge is performed. This will significantly reduce
            ! the NDCOF and thus the number of calls to YZK later.

            !*** Find the right counting parameter ***
            IF (IA .EQ. J) THEN
               ithis = IB
            ELSE
               ithis = IA
            ENDIF

            !*** Check it against the previously recorded ***
            IF (jblock .GT. 1) THEN
               DO i = 1, NDCOF
                  IF (NDA(i) .EQ. ithis) THEN   ! found, add the value
                     DA(i) = DA(i) + SUM
                     GOTO 123
                  ENDIF
               ENDDO
            ENDIF

            !*** Not found in the record, add an item ***
            NDCOF = NDCOF + 1
            IF (NDCOF .GT. NDDIM) THEN
               IF (NDDIM .GT. 0) THEN
                  CALL alcsca (PNTNDA, PNTRDA, NDDIM, 2)
               ELSE
                  CALL alcsca (PNTNDA, PNTRDA, NDDIM, 1)
               ENDIF
            ENDIF
            DA(NDCOF) = SUM ! print*, DA(NDCOF), ndcof, myid, 'myid'
            NDA(NDCOF) = ithis

  123       CONTINUE
         ENDDO                ! For labels
      ENDDO                   ! For blocks
     
*=======================================================================
*   Generate YA and XA coefficients; these arise from the two-electron
*   integrals
*=======================================================================

      DO 16 NFILE = 32, 32 + KMAXF

*         ...Re-position file mcp.30

         REWIND (30)
         DO i = 1, 6
            READ (30)
         ENDDO

*=======================================================================
*   Loop over blocks again, this time, for V-coefficients
*=======================================================================

         DO jblock = 1, nblock
*           ...Read in IROW from file mcp.30
            READ (30) mcplab, jblockt, ncf
               IF (jblockt .NE. jblock) THEN
                  WRITE (istde,*) myname, ':3 jblockt .NE. jblock'
                  STOP
               ENDIF
            READ (30) NELMNT
            CALL alloc (PNIROW, NELMNT, 4)
            READ (30) (IENDCdum, I = myid + 1, NCF, nprocs), 
     &                (IROW(I), I = 1, NELMNT)

            READ (nfile) mcplab, jblockt, ncf, ncoeff
               IF (jblockt .NE. jblock) THEN
                  WRITE (istde,*) myname, ':4 jblockt .NE. jblock'
                  STOP
               ENDIF

            K = NFILE - 32    ! multipolarity of the integral

*=======================================================================
*   Attempt to read another block of data
*=======================================================================

  999       READ (NFILE, IOSTAT = IOS) LAB, NCONTR
*
            IF (lab .EQ. 0 .AND. ncontr .EQ. 0) THEN
               CALL dalloc (pnirow)
               CYCLE
            ENDIF
            !***                       k
            !*** Decode the labels of R (abcd)
            INDEXS(4) = MOD (LAB, KEY)
            LAB = LAB / KEY
            INDEXS(2) = MOD (LAB, KEY)
            LAB = LAB / KEY
            INDEXS(3) = MOD (LAB, KEY)
            INDEXS(1) = LAB / KEY

            !*** Determine the number of indices that match
            IRANK = 0
            DO I = 1, 4
               IF (INDEXS(I) .EQ. J) IRANK = IRANK + 1
            ENDDO

            IF (IRANK .EQ. 0) THEN
               READ (nfile)
               GOTO 999
            ENDIF

            !*** At least one subshell index matches; allocate storage
            !*** for reading in the rest of this block
            IF (NCONTR .GT. NDIM) THEN
               CALL dalloc (PCOEFF)
               CALL dalloc (PICLMN)
               CALL dalloc (PINDEX)
               NDIM = NCONTR
               CALL alloc (PCOEFF, NDIM, 8)
               CALL alloc (PICLMN, NDIM, 4)
               CALL alloc (PINDEX, NDIM, 4)
            ENDIF

            !*** Read column index, sparse matrix index, and 
            !*** coefficient for all contributions from this integral
            READ (NFILE) (ICLMN(I), indx(I), COEFF(I), I = 1, NCONTR)

            !*** Add up all the contributions from this integral; 
            !*** off-diagonal contributions have double the weight
            SUM = 0.D0
            DO I = 1, NCONTR
               IR = IROW(indx(I))
               IC = ICLMN(I)
               CONTR = dsubrs (EOL, IR, IC, jblock) * COEFF(I)
               IF (IR .NE. IC) CONTR = CONTR + CONTR
               SUM = SUM + CONTR
            ENDDO
            SUM = 0.5D0 * SUM / UCFJ

            IF (IRANK .EQ. 1) THEN

*=======================================================================
*   One matching index: exchange potential contribution
*=======================================================================

               !*** Similar to DA, find ithis ***
               ithis = -911   ! initialize to an impossible value
                              ! though not necessary
               DO IIND = 1, 4
                  IF (INDEXS(IIND) .EQ. J) THEN    ! at least one
                     IL = IIND + 2
                     IF (IL .GT. 4) IL = IL - 4
                     IORB = INDEXS(IL)
                     IL = IIND + 1
                     IF (IL .GT. 4) IL = IL - 4
                     IYO1 = INDEXS(IL)
                     IL = IIND + 3
                     IF (IL .GT. 4) IL = IL - 4
                     IYO2 = INDEXS(IL)
                     ithis = (((IORB*KEY + IYO2)*KEY) + IYO1)*KEY + K
                     EXIT
                  ENDIF
               ENDDO

               IF (ithis .EQ. -911) STOP 'ithis .EQ. -911'

               !*** Check ithis against the previously recorded ***
               IF (jblock .GT. 1) THEN
                  DO i = 1, NXCOF
                     IF (NXA(i) .EQ. ithis) THEN
                        XA(i) = XA(i) + SUM
                        GOTO 999
                     ENDIF
                  ENDDO
               ENDIF

               !*** Not found in records, add an item ***
               NXCOF = NXCOF + 1
               IF (NXCOF .GT. NXDIM) THEN
                  IF (NXDIM .GT. 0) THEN
                     CALL alcsca (PNTNXA, PNTRXA, NXDIM, 2)
                  ELSE
                     CALL alcsca (PNTNXA, PNTRXA, NXDIM, 1)
                  ENDIF
               ENDIF
               XA(NXCOF) = SUM
               NXA(NXCOF) = ithis

            ELSEIF (IRANK .EQ. 2) THEN

*=======================================================================
*   Two matching indices: either direct or exchange potential
*   contribution
*=======================================================================

               IFOUND = 0
               DO IIND = 1, 4
                  IF (INDEXS(IIND) .EQ. J) THEN
                     IF (IFOUND .EQ. 0) THEN
                        LOC1 = IIND
                        IFOUND = IFOUND + 1
                     ELSEIF (IFOUND .EQ. 1) THEN
                        LOC2 = IIND
                        GOTO 234
                     ENDIF
                  ENDIF
               ENDDO

  234          IF (LOC2 - LOC1 .EQ. 2) THEN
*
*   Direct contribution
*
                  !*** Find ithis ***
                  IL = LOC1 + 3
                  IF (IL .GT. 4) IL = IL - 4
                  IYO2 = INDEXS(IL)
                  IL = LOC1 + 1
                  IF (IL .GT. 4) IL = IL - 4
                  IYO1 = INDEXS(IL)
                  ithis =  (IYO2 * KEY + IYO1) * KEY + K

                  !*** Check it against the previously recorded ***
                  IF (jblock .GT. 1) THEN
                     DO i = 1, NYCOF
                        IF (NYA(i) .EQ. ithis) THEN
                           YA(i) = YA(i) + SUM + SUM
                           GOTO 999
                        ENDIF
                     ENDDO
                  ENDIF

                  !*** Not found, add an item ***
                  NYCOF = NYCOF + 1
                  IF (NYCOF .GT. NYDIM) THEN
                     IF (NYDIM .GT. 0) THEN
                        CALL alcsca (PNTNYA, PNTRYA, NYDIM, 2)
                     ELSE
                        CALL alcsca (PNTNYA, PNTRYA, NYDIM, 1)
                     ENDIF
                  ENDIF
                  YA(NYCOF) = SUM + SUM
                  NYA(NYCOF) = ithis

               ELSE
*
*   Exchange contribution
*
                  !*** Find ithis ***
                  IL = LOC1 + 2
                  IF (IL .GT. 4) IL = IL - 4 
                  IORB = INDEXS(IL)
                  IL = LOC1 + 1
                  IF (IL .GT. 4) IL = IL - 4
                  IYO1 = INDEXS(IL)
                  IL = LOC1 + 3
                  IF (IL .GT. 4) IL = IL - 4
                  IYO2 = INDEXS(IL)
                  ithis = (((IORB*KEY + IYO2)*KEY) + IYO1)*KEY + K

                  !*** Check it against the previously recorded ***
                  IF (jblock .GT. 1) THEN
                     DO i = 1, NXCOF
                        IF (NXA(i) .EQ. ithis) THEN
                           XA(i) = XA(i) + SUM + SUM
                           GOTO 999
                        ENDIF
                     ENDDO
                  ENDIF

                  !*** Not found, add an item ***
                  NXCOF = NXCOF + 1
                  IF (NXCOF .GT. NXDIM) THEN
                     IF (NXDIM .GT. 0) THEN
                        CALL alcsca (PNTNXA, PNTRXA, NXDIM, 2)
                     ELSE
                        CALL alcsca (PNTNXA, PNTRXA, NXDIM, 1)
                     ENDIF
                  ENDIF
                  XA(NXCOF) = SUM + SUM
                  NXA(NXCOF) = ithis

               ENDIF

            ELSEIF (IRANK .EQ. 3) THEN

*=======================================================================
*   Three matching indices: direct and exchange potential contributions
*=======================================================================

               !*** Find ithis AND ithis2
               ithis = -911
               ithis2 = -911
               DO IIND = 1, 4
                  IF (INDEXS(IIND) .NE. J) THEN
                     INDIND = INDEXS(IIND)
                     IYO2 = INDIND
                     IYO1 = J
                     ithis = (IYO2 * KEY + IYO1) * KEY + K
                     IORB = INDIND
                     IYO1 = J
                     IYO2 = J
                     ithis2 = (((IORB*KEY+IYO2)*KEY)+IYO1)*KEY+K
                  ENDIF
               ENDDO

               IF (ithis .EQ. -911 .OR. ithis2 .EQ. -911) STOP 'ithis2'

               !*** Check the previously recorded for YA
               IF (jblock .GT. 1) THEN
                  DO i = 1, NYCOF
                     IF (NYA(i) .EQ. ithis) THEN
                        YA(i) = YA(i) + SUM + SUM
                        GOTO 456
                     ENDIF
                  ENDDO
               ENDIF

               ! Not found, add an item ***
               NYCOF = NYCOF + 1
               IF (NYCOF .GT. NYDIM) THEN
                  IF (NYDIM .GT. 0) THEN
                     CALL alcsca (PNTNYA, PNTRYA, NYDIM, 2)
                  ELSE
                     CALL alcsca (PNTNYA, PNTRYA, NYDIM, 1)
                  ENDIF
               ENDIF
               YA(NYCOF) = SUM + SUM
               NYA(NYCOF) = ithis

  456          CONTINUE

               !*** Check the previously recorded for XA
               IF (jblock .GT. 1) THEN
                  DO i = 1, NXCOF
                     IF (NXA(i) .EQ. ithis2) THEN
                        XA(i) = XA(i) + SUM
                        GOTO 999
                     ENDIF
                  ENDDO
               ENDIF

               ! Not found, add an item ***
               NXCOF = NXCOF + 1
               IF (NXCOF .GT. NXDIM) THEN
                  IF (NXDIM .GT. 0) THEN
                     CALL alcsca (PNTNXA, PNTRXA, NXDIM, 2)
                  ELSE
                     CALL alcsca (PNTNXA, PNTRXA, NXDIM, 1)
                  ENDIF
               ENDIF
               XA(NXCOF) = SUM
               NXA(NXCOF) = ithis2

            ELSEIF (IRANK .EQ. 4) THEN

*=======================================================================
*   Four matching indices: direct potential contribution
*=======================================================================

               !*** Find ithis AND ithis2
               IYO2 = J
               IYO1 = J
               ithis =  (IYO2 * KEY + IYO1) * KEY + K

               !*** Check the previously recorded for YA
               IF (jblock .GT. 1) THEN
                  DO i = 1, NYCOF
                     IF (NYA(i) .EQ. ithis) THEN
                        YA(i) = YA(i) + 4.D0 * SUM
                        GOTO 999
                     ENDIF
                  ENDDO
               ENDIF

               ! Not found, add an item ***
               NYCOF = NYCOF + 1
               IF (NYCOF .GT. NYDIM) THEN
                  IF (NYDIM .GT. 0) THEN
                     CALL alcsca (PNTNYA, PNTRYA, NYDIM, 2)
                  ELSE
                     CALL alcsca (PNTNYA, PNTRYA, NYDIM, 1)
                  ENDIF
               ENDIF
               YA(NYCOF) = 4.D0 * SUM
               NYA(NYCOF) = ithis

            ENDIF

            GOTO 999
         ENDDO          ! loop for V-Coefficients
   16 CONTINUE       ! loop for NFILE (KMAXF)

*=======================================================================
*   Deallocate storage for arrays local to this routine
*=======================================================================
      CALL dalloc (PCOEFF)
      CALL dalloc (PICLMN)
      CALL dalloc (PINDEX)

*=======================================================================
*   Assemble data from different nodes
*=======================================================================
!
!      ! This algorithm is based on the assumption that every node
!      ! has the same index arrays nda(), nxa() and nya()
!
!      WRITE (idstring, '(I3.3)') myid
!
!      WRITE (msg, '(I6.6)') ndcof
!      msg = 'ndcof = ' // msg(1:6) // ' on ' // idstring
!      CALL mpix_printmsg (msg, myid, nprocs)
!
!      WRITE (msg, '(I6.6)') nxcof
!      msg = 'nxcof = ' // msg(1:6) // ' on ' // idstring
!      CALL mpix_printmsg (msg, myid, nprocs)
!
!
!      WRITE (msg, '(I6.6)') nycof
!      msg = 'nycof = ' // msg(1:6) // ' on ' // idstring
!      CALL mpix_printmsg (msg, myid, nprocs)
!
!   Borrow the pointer just deallocated as the mpi buffers
!
!      !i = MAX (ndcof, nxcof - nxcof0, nycof - nycof0)
!      i = MAX (ndcof, nxcof, nycof)
!      IF (i .LE. 0) RETURN
!
!      CALL alloc (PCOEFF, i, 8)
!
!      IF (ndcof .GT. 0) THEN
!         CALL MPI_Allreduce (da(1), coeff(1), ndcof,
!     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!         CALL dcopy (ndcof, coeff(1), 1, da(1), 1)
!      ENDIF
!
!      IF (nxcof .GT. 0) THEN
!         CALL MPI_Allreduce (xa(1), coeff(1), nxcof, 
!     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!         CALL dcopy (nxcof, coeff(1), 1, xa(1), 1)
!      ENDIF
!
!      IF (nycof .GT. 0) THEN
!         CALL MPI_Allreduce (ya(1), coeff(1), nycof, 
!     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!         CALL dcopy (nycof, coeff(1), 1, ya(1), 1)
!      ENDIF
!
!      CALL dalloc (PCOEFF)

      RETURN
      END
