************************************************************************
*                                                                      *
      SUBROUTINE SETLAGmpi (EOL)
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL EOL
*                                                                      *
*   Sets up the data structure  pertaining to the Lagrange multipli-   *
*   ers  on the first entry;  on subsequent calls it  determines new   *
*   estimates for the multipliers.                                     *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, QUAD, RINTI,                           *
*               [RSCF92]: DACON, SETCOF, XPOT,  YPOT.                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 18 Dec 1992   *
*   MPI version by Xinghong He              Last update: 03 Aug 1998   *
*                                                                      *
************************************************************************
*
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      COMMON/CORE/NCORE

      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*2 NH
      LOGICAL FIRST,FIXLI,FIXLJ,FULLI,FULLJ,LFIX
*
      POINTER (PNTECV,ECV(*))
      POINTER (PNIECC,IECC(*))
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      DIMENSION YPJ(NNNP),YPM(NNNP),
     :          XPJ(NNNP),XPM(NNNP),
     :          XQJ(NNNP),XQM(NNNP)
*
      COMMON/DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /LAGR/PNTECV,PNIECC,NEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF1/UCF(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /ORBA/IORDER(NNNW)
*
      DATA FIRST/.TRUE./

      PARAMETER (P001 = 1.0D-01)
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

      IF (.NOT. FIRST) GOTO 123

*=======================================================================
*   Determine the total number of Lagrange multipliers and store
*   their indeces in IECC(1:NEC). Memories are allocated for IECC
*   and ECV. The only outputs are NEC and IECC(). This is done only
*   for the first time this routine is called.
*
*   Implementation changes:
*   The "DO LIraw ..." is executed twice (if NEC obtained in the 
*   first time > 0) so the memory allocation is done only once.
*
*   This part is not distributed.
*=======================================================================
      EPS = ACCY*0.01D0    ! criterion to see if an orb is occupied
      DO itwice = 1, 2
         NEC = 0
         DO LIraw = 1, NW - 1
            LI = iorder(LIraw)
            LIP1 = MAX (NCORE, LIraw) + 1
            NAKLI = NAK(LI)
            FIXLI = LFIX(LI)
            FULLI = ABS ( UCF(LI)-DBLE (NKJ(LI)+1) ) .LT. EPS
            DO LJraw = LIP1, NW
               LJ = iorder(LJraw)
               FIXLJ = LFIX(LJ)
               FULLJ = ABS ( UCF(LJ)-DBLE (NKJ(LJ)+1) ) .LT. EPS
               IF ( (NAK(LJ) .EQ. NAKLI) .AND.
     &               (.NOT. (FIXLI .AND. FIXLJ)) .AND.
     &               (.NOT. (FULLI .AND. FULLJ)) ) THEN
                  NEC = NEC + 1
                  !*** Encode index at 2nd round ***
                  IF (itwice .EQ. 2) IECC(NEC) = LI + KEY * LJ
               ENDIF
            ENDDO
         ENDDO
         IF (itwice .EQ. 1 .AND. NEC .GT. 0) THEN
            CALL ALLOC (PNTECV, NEC, 8)
            CALL ALLOC (PNIECC, NEC, 4)
         ELSE
            EXIT
         ENDIF
      ENDDO !itwice

*=======================================================================
*   Print information about Lagrange multipliers
*=======================================================================

      IF (myid .EQ. 0) THEN
         IF (NEC .EQ. 0) THEN
            WRITE (*,302)
         ELSE
            WRITE (*,304)
            DO LI = 1, NEC
               !*** Decode index ***
               IECCLI = IECC(LI)
               L1 = IECCLI / KEY
               L2 = IECCLI - KEY * L1
               WRITE (*,305) NP(L2),NH(L2),NP(L1),NH(L1)
            ENDDO
         ENDIF
      ENDIF

      FIRST = .FALSE.

CFF+GG  12/07/05
C     Lagrange multipliers need to be computed also on the first call
C     RETURN

  123 CONTINUE

*=======================================================================
*   Compute Lagrange multipliers for all pairs found above
*=======================================================================

      IF (NEC .EQ. 0) RETURN
      IF (myid .EQ. 0) WRITE (*,306)
      JLAST = 0
      MLAST = 0

      DO 10 LI = 1, NEC
         !*** Decode index ***
         IECCLI = IECC(LI)
         M = IECCLI / KEY
         J = IECCLI - KEY * M
*
         IF (J .NE. JLAST) THEN
            UCFJ = UCF(J)
            CALL SETCOF (EOL, J)
            CALL YPOT (J)
            CALL XPOT (J)
            CALL DACON
            DO I = 1, N
               YPJ(I) = YP(I)
               XPJ(I) = XP(I)
               XQJ(I) = XQ(I)
            ENDDO
            JLAST = J
         ENDIF
*
         IF (M .NE. MLAST) THEN
            UCFM = UCF(M)
            CALL SETCOF (EOL, M)
            CALL YPOT (M)
            CALL XPOT (M)
            CALL DACON
            DO I = 1, N
               YPM(I) = YP(I)
               XPM(I) = XP(I)
               XQM(I) = XQ(I)
            ENDDO
            MLAST = M
         ENDIF
*
         MTP = MAX (MF(J), MF(M))
*
         IF (LFIX(M)) THEN
            TA(1) = 0.D0
            DO I = 2, MTP
               TA(I) = RPOR(I)*( ( PF(I,M)*XQJ(I)
     :                            -QF(I,M)*XPJ(I) )
     :                          *C
     :                          +( PF(I,M)*PF(I,J)
     :                            +QF(I,M)*QF(I,J) ) * YPJ(I) )
            ENDDO
            CALL QUAD (RESULT)
            RIMJ = RINTI (M, J, 1)
            ECV(LI) = (RESULT - RIMJ / nprocs) * UCFJ
         ELSEIF (LFIX(J)) THEN
            TA(1) = 0.D0
            DO I = 2, MTP
               TA(I) = RPOR(I)*( ( PF(I,J)*XQM(I)
     :                            -QF(I,J)*XPM(I) )
     :                          *C
     :                          +( PF(I,J)*PF(I,M)
     :                            +QF(I,J)*QF(I,M) ) * YPM(I) )
            ENDDO
            CALL QUAD (RESULT)
            RIJM = RINTI (J, M, 1)
            ECV(LI) = (RESULT - RIJM / nprocs) * UCFM
         ELSE
            QDIF = ABS ( (UCFJ - UCFM) / MAX (UCFJ, UCFM) )
            IF (QDIF .GT. P001) THEN
               OBQDIF = 1.D0 / UCFJ - 1.D0 / UCFM
               TA(1) = 0.D0
               DO I = 2, MTP
                  TA(I) = RPOR(I)*( ( PF(I,M)*XQJ(I)
     :                               -QF(I,M)*XPJ(I)
     :                               -PF(I,J)*XQM(I)
     :                               +QF(I,J)*XPM(I) )
     :                             *C
     :                             +( YPJ(I)-YPM(I) )
     :                             *( PF(I,M)*PF(I,J)
     :                               +QF(I,M)*QF(I,J) ) )
               ENDDO
               CALL QUAD (RESULT)
               ECV(LI) = RESULT / OBQDIF
            ELSE
               OBQSUM = 1.D0 / UCFJ + 1.D0 / UCFM
               TA(1) = 0.D0
               DO I = 2, MTP
                  TA(I) = RPOR(I)*( ( PF(I,M)*XQJ(I)
     :                               -QF(I,M)*XPJ(I)
     :                               +PF(I,J)*XQM(I)
     :                               -QF(I,J)*XPM(I) )
     :                             *C
     :                             +( YPJ(I)+YPM(I) )
     :                             *( PF(I,M)*PF(I,J)
     :                               +QF(I,M)*QF(I,J) ) )
               ENDDO
               CALL QUAD (RESULT)
               RIMJ = RINTI (M, J, 1)
               ECV(LI) = (RESULT - 2.D0 * RIMJ / nprocs) / OBQSUM
            ENDIF
         ENDIF

*=======================================================================
*  Collect contributions from all nodes.
*  Another alternative is to modify mcpmpi to let every node
*  have the same set of "LAB" and "NCONTR". So an easier aggregation
*  can be done for arrays DA, XA, YA in setcof. And so the computation
*  of YPOT, XPOT, DACON can be distributed.
*=======================================================================
!Rasa start
!     print*,"setlagmpi: myid=",myid," before reduct.ECV(LI)=",ECV(LI)
!     call flush(6)
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr); 
!     if (myid .EQ. 0) read *, aaa 
!Rasa end

         CALL MPI_Allreduce (ECV(LI), tmp, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
         ECV(LI) = tmp
!Rasa start
!     print*,"setlagmpi: myid=",myid," after reduct. ECV(LI)=", ECV(LI)
!     call flush(6)
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr); 
!     if (myid .EQ. 0) read *, aaa 
!Rasa end

         IF (myid .EQ. 0) WRITE (*,307) NP(J),NH(J),NP(M),NH(M),ECV(LI)
         !call flush(6)
         

   10 CONTINUE

  302 FORMAT (/'Lagrange multipliers are not required')
  304 FORMAT (/'Include Lagrange multipliers between:'/)
  305 FORMAT (13X,2(2X,1I2,1A2))
  306 FORMAT (/'Lagrange multipliers:'/)
  307 FORMAT (13X,2(2X,1I2,1A2),2X,1PD16.9)

      RETURN
      END
