************************************************************************
*                                                                      *
      SUBROUTINE genintrk (myid, nprocs, N, j2max)
*                                                                      *
*  Input:
*     myid, nprocs
*  Output:
*   N - Number of integrals
*   j2max - max of 2*J
*
*       Generate the list of Rk Integrals that could arise from a set  *
*       of orbitals.                                                   *
*                                                                      *
*     Written by Per Jonsson                                           *
*   MPI version by Xinghong He            Last revision: 22 Jun 1998   *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120) 
      PARAMETER (KMAX = 20)

      LOGICAL GEN,TRIANGRK
      POINTER(PNTRIQ,RIQDUMMY(*))
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)  
     :      /CTEILSRK/PCTEILRK,PCTEVLRK
     :      /KKSTART/KSTART(0:KMAX)   

      POINTER (PCTEVLRK,VALTEIRK(*))
      POINTER (PCTEILRK, INDTEIRK(*))

!-----------------------------------------------------------------------

      KEY = NW + 1
      KSTART(0) = 1
*
*   Find 2*JMAX; it should not be greater than PARAMETER KMAX
*
      J2MAX = NKJ(1)
      DO I = 2, NW
         IF (NKJ(I) .GT. J2MAX) J2MAX = NKJ(I)
      ENDDO

      IF (J2MAX .GT. KMAX) THEN
         STOP 'genintrk: KMAX too small'
      ENDIF
*
*   Make the RK integrals: IA <= IB, IA <= IC, IA <= ID, IB <= ID
*   When GEN is false, sweep through to find dimension
*
      GEN = .FALSE.

  999 N = 0
      DO K = 0, J2MAX
         DO IA = 1, NW
            DO IB = IA, NW
               DO IC = IA, NW
                  IF (TRIANGRK(NKL(IA),K,NKL(IC))) THEN
                     DO ID = IB, NW
                        IF (TRIANGRK(NKL(IB),K,NKL(ID))) THEN
                           N = N + 1
                           IF (GEN .AND. (MOD(N,nprocs) .EQ. myid)) THEN
                              INDTEIRK(N) = ((IA*KEY+IB)*KEY+IC)*KEY+ID
                              VALTEIRK(N) = SLATER (IA,IB,IC,ID,K)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         KSTART(K+1) = N + 1
      ENDDO
*
*     Allocate memory for integral book keeping
*
      IF (.NOT. GEN) THEN
         CALL ALLOC (PCTEILRK,N,4)
         CALL ALLOC (PCTEVLRK,N,8)

! Initialization is necessary in the mpi version

         DO i = 1, N
            INDTEIRK(i) = 0
            VALTEIRK(i) = 0.d0
         ENDDO

         IF (myid .EQ. 0) 
     &      PRINT *, 'Computing ',N,' Rk integrals'

         GEN = .TRUE.
         GOTO 999
      ENDIF

      RETURN
      END
