************************************************************************
*                                                                      *
      LOGICAL FUNCTION TRIANGBREIT1 (IA,IB,IC,ID,K)
*                                                                      *
*     Evaluate the triangular relation for the Breit integral of       *
*     See paper by I. Grant and B McKenzie, JPB 13, (1980) 2671-2681   *
*     formula 5                                                        *
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

      TRIANGBREIT1 = .FALSE.

      DO L = K-1,K+1
         IF ((MOD(K+LA+LB,2).EQ.1).AND.(ABS(JA-JB).LE.2*L).AND.
     :      (ABS(JA+JB).GE.2*L)) THEN
            IF ((MOD(K+LC+LD,2).EQ.1).AND.(ABS(JC-JD).LE.2*L).AND.
     :         (ABS(JC+JD).GE.2*L)) THEN
               TRIANGBREIT1 = .TRUE.
               RETURN
            END IF
         END IF
      END DO

      RETURN
      END

