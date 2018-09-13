************************************************************************
*                                                                      *
      LOGICAL FUNCTION TRIANGBREIT2 (IA,IB,IC,ID,L)
*                                                                      *
*     Evaluate the triangle relations for Breit integrals of type 2    *
*     See Grant and McKenzier JPB 13 (1980), 2675-2681 formula 5       *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      include 'parameters.def'
      COMMON/ORB5/NKL(NNNW),NKJ(NNNW)
*
      LA = NKL(IA)
      LB = NKL(IB)
      LC = NKL(IC)
      LD = NKL(ID)
      JA = NKJ(IA)
      JB = NKJ(IB)
      JC = NKJ(IC)
      JD = NKJ(ID)

      TRIANGBREIT2 = .FALSE.

      IF ((MOD(L+LA+LB,2).EQ.0).AND.(ABS(JA-JB).LE.2*L).AND.
     :   (ABS(JA+JB).GE.2*L)) THEN
         IF ((MOD(L+LC+LD,2).EQ.0).AND.(ABS(JC-JD).LE.2*L).AND.
     :      (ABS(JC+JD).GE.2*L)) THEN
            TRIANGBREIT2 = .TRUE.
            RETURN
         END IF
      END IF

      RETURN
      END

