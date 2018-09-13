************************************************************************
*                                                                      *
*      SUBROUTINE RINTDENSVEC (I,J,DINT1VEC,NRNUC)
      SUBROUTINE RINTDENSVEC (I,J,DINT1VEC)
*                                                                      *
*   DINT1VEC CONTAINS                                                  *
*                                                                      *
*                                                                      *
*      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) |                     *
*                       I     J       I     J                          *
*                                                                      *
*                                                                      *
*                                                                      *
*   Written by Per Jonsson  and Jörgen Ekman                           *
*   Last revision: 24 Dec 2013                                         *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
*
      DIMENSION DINT1VEC(NNNW,NNNW,N)
      POINTER (PNTRPF,PF(NNNP,*)),(PNTRQF,QF(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /DEF9/CVAC,PI
!      DO L = 2,NRNUC
      DO L = 2,N
         DINT1VEC(I,J,L) = (PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))/
     :             (4.0D 00*PI*R(L)*R(L))
*         write(*,*) I, J, L, DINT1VEC(I,J,L)
*          write(*,*) L, R(L), I, J, PF(L,I)/R(L), QF(L,I)/R(L)
      END DO
*
      RETURN
      END
