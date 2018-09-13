************************************************************************
*                                                                      *
      SUBROUTINE genintbreit1 (myid, nprocs, NB, j2max)
*                                                                      *
*  Input:
*     myid, nprocs
*  Output:
*   N - Number of integrals
*   j2max - max of 2*J
*
*       Generate the list of Breit Integrals of type 1                 *
*       that could arise from a set of orbitals                        *
*       of orbitals.                                                   *
*                                                                      *
*       This routine is similar to genintrk                            *
*                                                                      *
*     Written by Per Jonsson            October 2014                   *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      PARAMETER (KMAX = 20)

      LOGICAL GEN,TRIANGBREIT1
      POINTER(PNTRIQ,RIQDUMMY(*))
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /KKSTARTBREIT/KSTARTBREIT1(0:KMAX),KSTARTBREIT2(0:KMAX)

      POINTER (PINDT1,INDTP1(*))
      POINTER (PVALT1,VALTP1(*))
      POINTER (PINDT2,INDTP2(*))
      POINTER (PVALT2,VALTP2(*))

      POINTER (PINDT3,INDTP3DUMMY(*))
      POINTER (PVALT3,VALTP3DUMMY(*))
      POINTER (PINDT4,INDTP4DUMMY(*))
      POINTER (PVALT4,VALTP4DUMMY(*))
      POINTER (PINDT5,INDTP5DUMMY(*))
      POINTER (PVALT5,VALTP5DUMMY(*))
      POINTER (PINDT6,INDTP6DUMMY(*))
      POINTER (PVALT6,VALTP6DUMMY(*))


!-----------------------------------------------------------------------

      KEY = NW + 1
      KSTARTBREIT1(0) = 1
*
*   Find 2*JMAX; it should not be greater than PARAMETER KMAX
*
      J2MAX = NKJ(1)
      DO I = 2, NW
         IF (NKJ(I) .GT. J2MAX) J2MAX = NKJ(I)
      ENDDO

      IF (J2MAX .GT. KMAX) THEN
         STOP 'genintbreit1: KMAX too small'
      ENDIF
*
*   Make the breit integrals. Note that there are no symmetry relations that
*   can be utilized. See paper by Grasnt.
*   When GEN is false, sweep through to find dimension
*
      GEN = .FALSE.

  999 NB = 0
      DO K = 0, J2MAX
         DO IA = 1, NW
            DO IB = 1, NW
               DO IC = 1, NW
                  DO ID = 1, NW
                     IF (TRIANGBREIT1(IA,IB,IC,ID,K)) THEN
                        NB = NB + 1
                        IF (GEN .AND.(MOD(NB,nprocs) .EQ. myid)) THEN
                           INDTP1(NB) = ((IA*KEY+IB)*KEY+IC)*KEY+ID
                           VALTP1(NB) = BRINTF (1,IA,IB,IC,ID,K)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         KSTARTBREIT1(K+1) = NB + 1
      ENDDO
*
*     Allocate memory for integral book keeping
*
      IF (.NOT. GEN) THEN
         CALL ALLOC (PINDT1,NB,4)
         CALL ALLOC (PVALT1,NB,8)

! Initialization is necessary in the mpi version

         DO i = 1, NB
            INDTP1(i) = 0
            VALTP1(i) = 0.d0
         ENDDO

         IF (myid .EQ. 0)
     &      PRINT *, 'Computing',NB,' Breit integrals of type 1'

         GEN = .TRUE.
         GOTO 999
      ENDIF

      RETURN
      END


