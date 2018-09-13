************************************************************************
*                                                                      *
      SUBROUTINE ORTHOR (J,INV)
! INV not initialized here but used as
! INV = INV + 1
!
*                                                                      *
*   This routine Schmidt orthogonalizes orbital  J  to  all orbitals   *
*   which  have  better  self-consistency.  Note that fixed orbitals   *
*   have the best self-consistency.                                    *
*                                                                      *
*   Call(s) to: [LIB92]: COUNT, RINT.                                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
*                                                                      *
************************************************************************
*
! A bug regarding IORDER fixed and new schemes introduced. See comments
! Anyway this routine was not used anywhere.
!XHH 1997.02.21
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL CHECK,NOINVT
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INVT/NOINVT(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /OVL/NOVL
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
!XHH Added common /orba/ and /iounit/
      LOGICAL LFIX, changed
      CHARACTER*2 NH
      COMMON/ORB10/NH(NNNW)
      COMMON/ORBA/IORDER(NNNW)
      COMMON/iounit/istdi,istdo,istde
      COMMON/FIXD/NFIX,LFIX(NNNW)

*
Cww      EPS = ACCY*0.1D 00
      EPS = ACCY*0.01D 00
*
      CHECK = .NOT. NOINVT(J)
*
!XHH
! Bug fixed. Here J is actually IORDER(J_raw), thus K should be
! treated the same way.
      changed = .FALSE.
!      DO 4 K = 1,NW
      DO 4 Kraw = 1,NW
         K = IORDER(Kraw)
!         write(istde,*) '***',kraw,k,np(k),nh(k),scnsty(k),'***'
*
!XHH orbitals with higher self-consistency are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) .AND.
!     :        (SCNSTY(K) .LT. SCNSTY(J)) ) THEN
!XHH All orbitals are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) ) THEN
!XHH orbitals before the current and unchanged ones are considered
         IF ( (NAK(K) .EQ. NAK(J)) .AND.
     :        ((K.LT.J) .OR. LFIX(K))   ) THEN

            changed = .TRUE.
*
*   Compute overlap
*
            OVRLAP = RINT (J,K,0)
*
*   Schmidt orthogonalise
*
            PZ(J) = PZ(J)-OVRLAP*PZ(K)
            MTP = MAX (MF(J),MF(K))
            DO 1 I = 1,MTP
               PF(I,J) = PF(I,J)-OVRLAP*PF(I,K)
               QF(I,J) = QF(I,J)-OVRLAP*QF(I,K)
    1       CONTINUE
         ENDIF
    4 CONTINUE
!XHH
      IF(changed) THEN
*
*   Normalise
*
            MF(J) = MTP
            DNORM = RINT (J,J,0)
            FACTOR = 1.0D 00/SQRT (DNORM)
*
*   Determine if inversion is necessary
*
            IF (CHECK) THEN
               CALL COUNT (PF(1,J),MTP,NNCFF,SGN)
               IF (SGN .LT. 0.0D 00) THEN
                  INV = INV+1
                  FACTOR = -FACTOR
               ENDIF
            ENDIF
*
*   Perform normalization and/or inversion
*
            PZ(J) = FACTOR*PZ(J)
            DO 2 I = 2,MTP
               PF(I,J) = FACTOR*PF(I,J)
               QF(I,J) = FACTOR*QF(I,J)
    2       CONTINUE
*
*   Find new MF(J)
*
            MTP = MTP+1
    3       MTP = MTP-1
            IF (ABS (PF(MTP,J)) .LT. EPS) THEN
               PF(MTP,J) = 0.0D 00
               QF(MTP,J) = 0.0D 00
               GOTO 3
            ELSE
               MF(J) = MTP
            ENDIF
*
!XHH two lines moved forward and added an ENDIF
!         ENDIF
!    4 CONTINUE
      ENDIF
*
      RETURN
      END
