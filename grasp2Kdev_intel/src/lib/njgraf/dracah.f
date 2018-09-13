************************************************************************
*                                                                      *
      SUBROUTINE DRACAH (I,J,K,L,M,N,RAC)
*                                                                      *
*   SUBROUTINE  to calculate Racah coefficients. The arguments I, J,   *
*   K, L, M, N should be twice their actual value. Works for integer   *
*   and  half-integer  values of  angular momenta. The routine makes   *
*   use of the GAM  array, thus  SUBROUTINE FACTT must be called be-   *
*   fore this routine is used.                                         *
*                                                                      *
*   Written by N S Scott                    Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      PARAMETER (MFACT = 500)
*
      COMMON/FACTS/GAM(MFACT)
*
      J1 = I+J+M
      J2 = K+L+M
      J3 = I+K+N
      J4 = J+L+N
      IF (((2*MAX (I,J,M)-J1) .GT. 0) .OR. (MOD (J1,2) .NE. 0)) GOTO 2
      IF (((2*MAX (K,L,M)-J2) .GT. 0) .OR. (MOD (J2,2) .NE. 0)) GOTO 2
      IF (((2*MAX (I,K,N)-J3) .GT. 0) .OR. (MOD (J3,2) .NE. 0)) GOTO 2
      IF (((2*MAX (J,L,N)-J4) .GT. 0) .OR. (MOD (J4,2) .NE. 0)) GOTO 2
      GOTO 1
   2  RAC = 0.0D 00
      RETURN
*
   1  CONTINUE
      J1 = J1/2
      J2 = J2/2
      J3 = J3/2
      J4 = J4/2
      J5 = (I+J+K+L)/2
      J6 = (I+L+M+N)/2
      J7 = (J+K+M+N)/2
      NUMIN = MAX (J1,J2,J3,J4)+1
      NUMAX = MIN (J5,J6,J7)+1
      RAC = 1.0D 00
      ICOUNT = 0
*
      IF (NUMIN .EQ. NUMAX) GOTO 4
      NUMIN = NUMIN+1
*
      DO 3 KK = NUMIN,NUMAX
         KI = NUMAX-ICOUNT
         RAC = 1.0D 00
     :       -(RAC*DBLE(KI*(J5-KI+2)*(J6-KI+2)*(J7-KI+2))/
     :      DBLE((KI-1-J1)*(KI-1-J2)*(KI-1-J3)*(KI-1-J4)))
         ICOUNT = ICOUNT+1
   3  CONTINUE
*
      NUMIN = NUMIN-1
   4  RAC = RAC*((-1.0D 00)**(J5+NUMIN+1))
     : *EXP( (GAM(NUMIN+1)-GAM(NUMIN-J1)
     : -GAM(NUMIN  -J2)-GAM(NUMIN  -J3)-GAM(NUMIN  -J4)-GAM(J5+2-NUMIN)
     : -GAM(J6+2-NUMIN)-GAM(J7+2-NUMIN))+((GAM(J1+1-I)+GAM(J1+1-J)
     : +GAM(J1+1-M)-GAM(J1+2)+GAM(J2+1-K)+GAM(J2+1-L)+GAM(J2+1-M)
     : -GAM(J2+2)+GAM(J3+1-I)+GAM(J3+1-K)+GAM(J3+1-N)-GAM(J3+2)
     : +GAM(J4+1-J)+GAM(J4+1-L)+GAM(J4+1-N)-GAM(J4+2))
     : *0.5D 00  ))
*
      RETURN
      END
