************************************************************************
*                                                                      *
      SUBROUTINE orthal
*                                                                      *
*   This routine Schmidt orthogonalises radial wavefunctions in the    *
*   order defined by an index array IORDER                  .          *
*                                                                      *
*                                                                      *
************************************************************************
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPR
      CHARACTER*2 NH

      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))

      DIMENSION OVLAP(NNNW),J(NNNW)

      COMMON/DEBUGR/LDBPR(30)
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /ORBA/IORDER(NNNW)

      LOGICAL changed

      EPS = 0.01D 00*ACCY

      DO 9 Lraw = 2,NW
         L = IORDER(Lraw)

         NAKL = NAK(L)
         KOUNT = 0

!    MTP0 introduced to count the maximum number of points during
!    orthogonalization of the L-th orbital to other orbitals
!    The logical variable changed is initialized to .F.

         MTP0 = MF(L)
         changed = .FALSE.

         DO 8 Kraw = 1, Lraw - 1
            K = IORDER(Kraw)

            IF (NAK(K) .EQ. NAKL) THEN

               changed = .TRUE.

               OVRLAP = RINT (L,K,0)

!   Schmidt orthogonalise

               KOUNT = KOUNT+1
               J(KOUNT) = K
               OVLAP(KOUNT) = OVRLAP

               PZ(L) = PZ(L)-OVRLAP*PZ(K);pRINT *,'orthal1',pz(j)
               MTP = MAX (MF(L),MF(K))
               MTP0 = MAX(MTP0, MF(K))

               DO I = 1,MTP
                  PF(I,L) = PF(I,L)-OVRLAP*PF(I,K)
                  QF(I,L) = QF(I,L)-OVRLAP*QF(I,K)
               ENDDO
            ENDIF
    8    CONTINUE

!   Normalise

!    Use MTP0 to replace MTP and only when the orbital is changed.
!    This is in accordance with the original version which had the
!    normalization etc within the inner K loop.

            IF(changed) THEN
               MTP = MTP0

               MF(L) = MTP
               DNORM = RINT (L,L,0)
               FACTOR = 1.0D 00/SQRT (DNORM)

               PZ(L) = FACTOR*PZ(L);pRINT *,'orthal2',pz(j)
               DO I = 2,MTP
                  PF(I,L) = FACTOR*PF(I,L)
                  QF(I,L) = FACTOR*QF(I,L)
               ENDDO

!   Find new MF(L)

               MTP = MTP+1
    7          MTP = MTP-1
               IF (ABS (PF(MTP,L)) .LT. EPS) THEN
                  PF(MTP,L) = 0.0D 00
                  QF(MTP,L) = 0.0D 00
                  GOTO 7
               ELSE
                  MF(L) = MTP
               ENDIF
            ENDIF

!---------------------------------------------------------------
! Check the orthonormalisation
! KOUNT and OVLAP are destroyed
! should be commented out after seeing the values
!
!      KOUNT = 0
!      DO kraw = 1, Lraw
!         k = iorder(kraw)
!         IF (NAK(K) .EQ. NAKL) THEN
!            KOUNT = KOUNT+1
!            OVLAP(KOUNT) = RINT (L,K,0)
!            J(KOUNT) = K
!         ENDIF
!      ENDDO
!      WRITE (*,301)
!     :      (OVLAP(I),NP(L),NH(L),NP(J(I)),NH(J(I)),I = 1,KOUNT)
!---------------------------------------------------------------
         IF (LDBPR(3) .AND. (KOUNT .GT. 0)) WRITE (99,301)
     :      (OVLAP(I),NP(L),NH(L),NP(J(I)),NH(J(I)),I = 1,KOUNT)

    9 CONTINUE

      RETURN

  301 FORMAT (1P,5(2X,1D10.3,' = <',1I2,1A2,'|',1I2,1A2,'>'))

      END
