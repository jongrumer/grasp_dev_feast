************************************************************************
*                                                                      *
      SUBROUTINE SETHAM (jblock, myid, nprocs)
      IMPLICIT REAL*8          (A-H, O-Z)
      INTEGER jblock, myid, nprocs
*                                                                      *
*   This  SUBROUTINE  sets up the  Hamiltonian matrix.
*   Works for serial and mpi versions.
*   Matrix is generated in lower-triangle-by-rows way.
*   Evaluation of the average energy is moved out to caller routine.
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CLRX, DALLOC, RINTI, SLATER.           *
*               [RSCF92]: FCO, GCO, IQ.                                *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 21 Dec 1992   *
*   Modified by Xinghong He               Last revision: 22 Jul 1998   *
*                                                                      *
************************************************************************
*

*      integer*8 nelmnt
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)

      !*** needs EMT, IENDC.
      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
      COMMON/HMAT/PNTEMT,PIENDC,PNIROW,NELMNT

      !*** needs ncf, nw
      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/ORB2/NCF,NW,PNTRIQ

      COMMON/MCPA/KMAXF
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)

      POINTER (pncfpast, ncfpast(1))
      POINTER (pncminpast, ncminpast(1))
      POINTER (pnevecpast, nevecpast(1))
      COMMON/pos/pncfpast,pncminpast,pnevecpast,ncftot,nvecsiz

      CHARACTER*3 mcplab
      CHARACTER*(*), PARAMETER:: myname = 'SETHAM'
      LOGICAL SET
      POINTER (PCOEFF,COEFF(1))
      !POINTER (PICLMN,ICLMNdum)
      POINTER (PINDX,INDX(1))

      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
      ncfpat = ncfpast(jblock)

*=======================================================================
*   Accumulate diagonal terms that do not require MCP coefficients
*=======================================================================

*=======================================================================
*   Piece involving I(a,a) integrals
*=======================================================================

      DO IA = 1, NW
         SET = .FALSE.
         DO IR = myid + 1, NCF, nprocs
            itmp = IQ (IA, IR + ncfpat)
            IF (itmp .LE. 0) CYCLE
            !*** Occupation number not zero ...
            IF (.NOT. SET) THEN
               DIAA = RINTI (IA, IA, 0)
               SET = .TRUE.
            ENDIF
            ! IDIAG = IENDC(IR-1)+1
            IDIAG = IENDC(IR)       ! lower-triangle-by-rows mode
            EMT(IDIAG) = EMT(IDIAG) + itmp * DIAA
         ENDDO
      ENDDO

*=======================================================================
*                    0
*   Piece involving F (a,a) integrals
*=======================================================================

      DO IA = 1, NW
         SET = .FALSE.
         DO IR = myid + 1, NCF, nprocs
            COEF = FCO (0, IR + ncfpat, IA, IA)
            IF (COEF .EQ. 0.D0) CYCLE
            !*** Angular coefficient not zero ...
            IF (.NOT. SET) THEN
               F0AA = SLATER (IA, IA, IA, IA, 0)
               SET = .TRUE.
            ENDIF
            ! IDIAG = IENDC(IR-1)+1
            IDIAG = IENDC(IR)
            EMT(IDIAG) = EMT(IDIAG) + COEF * F0AA
         ENDDO
      ENDDO

*=======================================================================
*                    k
*   Piece involving F (a,a) integrals
*=======================================================================

      KM = 0
      K = 0
    6 K = K + 2
      DO IA = 1, NW
         K0 = NKJ(IA) - 1
         IF (K0 .GT. KM) KM = K0
         IF (K .GT. K0) CYCLE
         SET = .FALSE.
         DO IR = myid + 1, NCF, nprocs
            COEF = FCO (K, IR + ncfpat, IA, IA)
            IF (COEF .EQ. 0.D0) CYCLE
            IF (.NOT. SET) THEN
               FKAA = SLATER (IA, IA, IA, IA, K)
               SET = .TRUE.
            ENDIF
            ! IDIAG = IENDC(IR-1)+1
            IDIAG = IENDC(IR)
            EMT(IDIAG) = EMT(IDIAG) + COEF * FKAA
         ENDDO
      ENDDO
      IF (K .LT. KM) GOTO 6

*=======================================================================
*                    0
*   Piece involving F (a,b) integrals
*=======================================================================

      DO 11 IA = 1, NW - 1
      DO 10 IB = IA + 1, NW
         SET = .FALSE.
         DO 9 IR = myid + 1, NCF, nprocs
            COEF = FCO (0, IR + ncfpat, IA, IB)
            IF (COEF .EQ. 0.D0) CYCLE
            IF (.NOT. SET) THEN
               F0AB = SLATER (IA, IB, IA, IB, 0)
               SET = .TRUE.
            ENDIF
            ! IDIAG = IENDC(IR-1)+1
            IDIAG = IENDC(IR)
            EMT(IDIAG) = EMT(IDIAG) + COEF * F0AB
    9    CONTINUE
   10 CONTINUE
   11 CONTINUE

*=======================================================================
*                    k
*   Piece involving G (a,b) integrals
*=======================================================================

      KM = 0
      K = -1
   12 K = K + 1
      DO 15 IA = 1, NW - 1
         NKJIA = NKJ(IA)
         DO 14 IB = IA + 1, NW
            NKJIB = NKJ(IB)
            SET = .FALSE.
            IF (NAK(IA) * NAK(IB) .GT. 0) THEN
               KMIN = ABS ((NKJIA - NKJIB) / 2)
            ELSE
               KMIN = ABS ((NKJIA - NKJIB) / 2) + 1
            ENDIF
            IF (MOD (K - KMIN, 2) .NE. 0) CYCLE

            KMAX = (NKJIA + NKJIB) / 2
            IF (KMAX .GT. KM) KM = KMAX
            IF ((K .LT. KMIN) .OR. (K .GT. KMAX)) CYCLE

            DO IR = myid + 1, NCF, nprocs
               COEF = GCO (K, IR + ncfpat, IA, IB)
               IF (COEF .EQ. 0.D0) CYCLE
               IF (.NOT. SET) THEN
                  GKAB = SLATER (IA, IB, IB, IA, K)
                  SET = .TRUE.
               ENDIF
               ! IDIAG = IENDC(IR-1)+1
               IDIAG = IENDC(IR)
               EMT(IDIAG) = EMT(IDIAG) + COEF * GKAB
            ENDDO
   14    CONTINUE
   15 CONTINUE
      IF (K .LT. KM) GOTO 12

*=======================================================================
*   Local storage for reading mcpXXX files
*=======================================================================

      NDIM = 1
      CALL ALLOC (PCOEFF, NDIM, 8)
      !CALL ALLOC (PICLMN, NDIM, 4)
      CALL ALLOC (PINDX, NDIM, 4)

*=======================================================================
*   Accumulate one-electron terms that require MCP coefficients
*=======================================================================

      READ (31) mcplab, jblockt, ncft, ncoeff

         IF (jblockt .NE. jblock ) THEN
            WRITE (istde,*) myname, ': blk1=', jblockt, ' blk2=', jblock
            STOP
         ENDIF
         IF (ncft .NE. ncf ) THEN
            WRITE (istde,*) myname, ': ncf1 = ', ncft, ' ncf2 = ', ncf
            STOP
         ENDIF

*=======================================================================
*   Loop over non-zero labels which have non-zero elements
*=======================================================================

      READ (31, IOSTAT = IOS) LAB, NCONTR
      IF (IOS .NE. 0) STOP 'IOS .NE. 0 when reading LAB, NCONTR'
      DO WHILE (LAB .NE. 0 .OR. NCONTR .NE. 0)

         !*** decode the label of I(ab)
         IA = MOD (LAB, KEY)
         IB = LAB / KEY

         !*** Compute radial integral I(ab)
         TEGRAL = RINTI (IA, IB, 0)

         ! Read column index, sparse matrix index, and coefficient
         ! for all contributions from this integral.
         IF (NCONTR .GT. NDIM) THEN
            CALL DALLOC (PCOEFF)
            !CALL DALLOC (PICLMN)
            CALL DALLOC (PINDX)
            NDIM = NCONTR
            CALL ALLOC (PCOEFF, NDIM, 8)
            !CALL ALLOC (PICLMN, NDIM, 4)
            CALL ALLOC (PINDX, NDIM, 4)
         ENDIF
         READ (31) (ICLMNdum, INDX(I), COEFF(I), I = 1, NCONTR)

         !*** Store all the contributions from this integral
         DO I = 1, NCONTR
            LOC = INDX(I)
            IF (LOC .GT. NELMNT) THEN
               PRINT *, '  Error in computing 1-e contribution'
               PRINT *, '  LOC = ', LOC, '  NELMNT = ', NELMNT
               STOP
            ENDIF
            EMT(LOC) = EMT(LOC) + TEGRAL * COEFF(I)
         ENDDO

         READ (31, IOSTAT = IOS) LAB, NCONTR
         IF (IOS .NE. 0) STOP 'IOS .NE. 0 when reading LAB, NCONTR'
      ENDDO

*=======================================================================
*   Accumulate two-electron terms that require MCP coefficients
*=======================================================================

      DO 20 NFILE = 32, 32 + KMAXF
         K = NFILE - 32

         READ (nfile) mcplab, jblockt, ncft, ncoeff

         IF (jblockt .NE. jblock ) THEN
            WRITE (istde,*) myname, ': blk3=', jblockt, ' blk4=', jblock
            STOP
         ENDIF
         IF (ncft .NE. ncf ) THEN
            WRITE (istde,*) myname, ': ncf3 = ', ncft, ' ncf4 = ', ncf
            STOP
         ENDIF

*=======================================================================
*   Loop over non-zero labels which have non-zero elements
*=======================================================================

         READ (NFILE, IOSTAT = IOS) LAB, NCONTR
         IF (IOS .NE. 0) STOP 'IOS .NE. 0 when reading LAB, NCONTR 2'
         DO WHILE (LAB .NE. 0 .OR. NCONTR .NE. 0)

            !                         k
            !*** decode the label of R (abcd)
            ID  = MOD (LAB, KEY)
            LAB = LAB / KEY
            IB  = MOD (LAB, KEY)
            LAB = LAB / KEY
            IC  = MOD (LAB, KEY)
            IA  = LAB / KEY

            !*** Compute radial integral
            TEGRAL = SLATER (IA, IB, IC, ID, K)

            ! Read column index, sparse matrix index, and coefficient
            ! for all contributions from this integral.
            IF (NCONTR .GT. NDIM) THEN
               CALL DALLOC (PCOEFF)
               !CALL DALLOC (PICLMN)
               CALL DALLOC (PINDX)
               NDIM = NCONTR
               CALL ALLOC (PCOEFF, NDIM, 8)
               !CALL ALLOC (PICLMN, NDIM, 4)
               CALL ALLOC (PINDX, NDIM, 4)
            ENDIF
            READ (NFILE) (ICLMNdum, INDX(I), COEFF(I), I = 1, NCONTR)

            !*** Store all the contributions from this integral
            DO I = 1, NCONTR
               LOC = INDX(I)
            IF (LOC .GT. NELMNT) THEN
               PRINT *, '  Error in computing 2-e contribution'
               PRINT *, '  LOC = ', LOC, '  NELMNT = ', NELMNT
               STOP
            ENDIF
               EMT(LOC) = EMT(LOC) + TEGRAL * COEFF(I)
            ENDDO

            READ (NFILE, IOSTAT = IOS) LAB, NCONTR
            IF (IOS .NE. 0) STOP 'IOS .NE. 0 when reading LAB, NCONTR 2'
         ENDDO
   20 CONTINUE

*=======================================================================
*   Deallocate local storage
*=======================================================================

      CALL DALLOC (PCOEFF)
      !CALL DALLOC (PICLMN)
      CALL DALLOC (PINDX)

      RETURN
      END
