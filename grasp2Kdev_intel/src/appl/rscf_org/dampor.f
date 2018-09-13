************************************************************************
*                                                                      *
      SUBROUTINE DAMPOR (J,INV,odampj)
*                                                                      *
*   This subroutine damps the orbital wave function with index J. it   *
*   also  stores  the  previous  determination  of this orbital.       *
*                                                                      *
*   Call(s) to: [LIB92]: COUNT, RINT.                                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PCDAMP
      POINTER (PCDAMP,DAMPDUMMY)
      LOGICAL CHECK,NOINVT
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /INVT/NOINVT(NNNW)
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Initialization
*
Cww      EPS = 0.1D 00*ACCY
      EPS = 0.01D 00*ACCY
      CHECK = .NOT. NOINVT(J)
*
*   Damp orbital J using the damping factor  ABS (ODAMP(J)):  ODAMP(J)
*   is restricted to the open interval (-1,1) by  DATSCF ; the meaning
*   of a  negative  ODAMP(J)  is that  ABS (ODAMP(J))  is the constant
*   damping factor;  a positive  ODAMP(J)  can be reset by  subroutine
*   dampck
*
*   Store previous determination of this orbital in arrays P, Q
*
!XHH odampj goes to the argument
!      ODAMPJ = ABS (ODAMP(J))
*
      IF (ODAMPJ .GT. EPS) THEN
*
         FACTOR = 1.0D 00-ODAMPJ
*
         PZ(J) = FACTOR*P0+ODAMPJ*PZ(J)
*
         MTPO = MF(J)
         MTPN = MTP0
         MTP0 = MTPO
*
         MTP = MAX (MTPN,MTPO)
         DO 1 I = 1,MTP
            POLDI = PF(I,J)
            PF(I,J) = FACTOR*P(I)+ODAMPJ*PF(I,J)
            P(I) = POLDI
            QOLDI = QF(I,J)
            QF(I,J) = FACTOR*Q(I)+ODAMPJ*QF(I,J)
            Q(I) = QOLDI
    1    CONTINUE
*
*   Compute normalization factor
*
         MF(J) = MTP
         DNORM = RINT (J,J,0)
         DNFAC = 1.0D 00/SQRT (DNORM)
*
*   Determine if inversion is necessary
*
         IF (CHECK) THEN
            CALL COUNT (PF(1,J),MTP,NNCFF,SGN)
            IF (SGN .LT. 0.0D 00) THEN
               INV = INV+1
               DNFAC = -DNFAC
            ENDIF
         ENDIF
*
         PZ(J) = PZ(J)*DNFAC
         DO 2 I = 1,MTP
            PF(I,J) = DNFAC*PF(I,J)
            QF(I,J) = DNFAC*QF(I,J)
    2    CONTINUE
*
*   Find new MF(J)
*
         MFJ = MTP+1
    3    MFJ = MFJ-1
         IF (ABS (PF(MFJ,J)) .LT. EPS) THEN
            PF(MFJ,J) = 0.0D 00
            QF(MFJ,J) = 0.0D 00
            GOTO 3
         ELSE
            MF(J) = MFJ
         ENDIF
*
      ELSE
*
         PZ(J) = P0
*
         MTPO = MF(J)
         MTPN = MTP0
         MTP0 = MTPO
*
         MTP = MAX (MTPN,MTPO)
         DO 4 I = 1,MTP
            POLDI = PF(I,J)
            PF(I,J) = P(I)
            P(I) = POLDI
            QOLDI = QF(I,J)
            QF(I,J) = Q(I)
            Q(I) = QOLDI
    4    CONTINUE
*
         MF(J) = MTP
*
      ENDIF
*
      RETURN
      END
