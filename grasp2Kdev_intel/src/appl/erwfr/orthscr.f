************************************************************************
*                                                                      *
      SUBROUTINE ORTHSCR
*                                                                      *
*   This routine Schmidt orthogonalises radial wavefunctions.          *
*   It is assumed that the orbitals is in reversed order and thus      *
*   the orthogonalization order is from the end of the list and        *
*   inwards.
*                                                                      *
*   Call(s) to: [LIB92]: RINT.                                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 14 Oct 1992   *
*   Modified by Per Jonsson, at Malmo, December 2013                   *
*  
*                                                                      *
************************************************************************
*
! Normalization of the orbitals moved out of the inner loop
! XHH 1997.02.14
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      DIMENSION OVLAP(NNNW),J(NNNW)
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)

!XHH
      LOGICAL changed
*
*   Set tabulated values of the radial wavefunction to zero
*   if they are less than EPS
*
Cww      EPS = 0.1D 00*ACCY
      EPS = 0.01D 00*ACCY
*
*   Determine the number of interesting overlaps
*
!XHH wasting time, no use at all
!      NOVL = 0
!      DO 2 K = 1,NW-1
!         NAKK = NAK(K)
!         DO 1 L = K+1,NW
!            IF (NAKK .EQ. NAK(L)) NOVL = NOVL+1
!    1    CONTINUE
!    2 CONTINUE
*
!      IF (NOVL .EQ. 0) RETURN
*
*CPer    3 DO 9 L = 2,NW
    3 DO 9 L = NW-1,1,-1
*
         NAKL = NAK(L)
         KOUNT = 0
!XHH MTP0 introduced to count the maximum number of points during
!    orthogonalization of the L-th orbital to other orbitals
!    The logical variable changed is initialized to .F.

         MTP0 = MF(L)
         changed = .FALSE.
*
* CPer    4    DO 8 K = 1,L-1
    4    DO 8 K = NW,L+1,-1
*
            IF (NAK(K) .EQ. NAKL) THEN
!XHH
               changed = .TRUE.
*
*   Compute overlap
*
               OVRLAP = RINT (L,K,0)
*
*   Schmidt orthogonalise
*
               KOUNT = KOUNT+1
               J(KOUNT) = K
               OVLAP(KOUNT) = OVRLAP
*
               PZ(L) = PZ(L)-OVRLAP*PZ(K)
               MTP = MAX (MF(L),MF(K))
               MTP0 = MAX(MTP0, MF(K))

               DO 5 I = 1,MTP
                  PF(I,L) = PF(I,L)-OVRLAP*PF(I,K)
                  QF(I,L) = QF(I,L)-OVRLAP*QF(I,K)
    5          CONTINUE
            ENDIF
    8    CONTINUE
*
*   Normalise
*
!XHH Use MTP0 to replace MTP and only when the orbital is changed.
!    This is in accordance with the original version which had the
!    normalization etc within the inner K loop.

            IF(changed) THEN
               MTP = MTP0

               MF(L) = MTP
               DNORM = RINT (L,L,0)
               FACTOR = 1.0D 00/SQRT (DNORM)
*
               PZ(L) = FACTOR*PZ(L)
               DO 6 I = 2,MTP
                  PF(I,L) = FACTOR*PF(I,L)
                  QF(I,L) = FACTOR*QF(I,L)
    6          CONTINUE
*
*   Find new MF(L)
*
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

!XHH Moved ahead
!            ENDIF
!    8    CONTINUE
*
*   Print overlap information
*
!---------------------------------------------------------------
! debug 97.02.14 done
! Check the orthonormalisation
! KOUNT and OVLAP are destroyed
! should be commented out after seeing the values
!
!      KOUNT = 0
!      DO k = 1, L
!            IF (NAK(K) .EQ. NAKL) THEN
!               KOUNT = KOUNT+1
!               OVLAP(KOUNT) = RINT (L,K,0)
!               J(KOUNT) = K
!            ENDIF
!      ENDDO
!      WRITE (*,301)
!     :      (OVLAP(I),NP(L),NH(L),NP(J(I)),NH(J(I)),I = 1,KOUNT)
!---------------------------------------------------------------
         IF (LDBPR(3) .AND. (KOUNT .GT. 0)) WRITE (99,301)
     :      (OVLAP(I),NP(L),NH(L),NP(J(I)),NH(J(I)),I = 1,KOUNT)
*
    9 CONTINUE
*
      RETURN
*
  301 FORMAT (1P,5(2X,1D10.3,' = <',1I2,1A2,'|',1I2,1A2,'>'))
*
      END
