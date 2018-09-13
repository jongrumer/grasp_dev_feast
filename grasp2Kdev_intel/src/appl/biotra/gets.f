************************************************************************
*                                                                      *
      SUBROUTINE GETS (S,NSHLI,NSHLF)
*                                                                      *
*   This subroutine calculates the initial and final state             *
*   radial overlap matrix                                              *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION S(NSHLI,NSHLF)
      CHARACTER*2 NHII,NHFF

      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /TATB/TA(NNN1),TB(NNN1),MTP

      POINTER (PNTRPFII,PFII(NNNP,*)),(PNTRQFII,QFII(NNNP,*))

      POINTER (PNTRPFFF,PFFF(NNNP,*)),(PNTRQFFF,QFFF(NNNP,*))
*
*  Initial state common
*
      COMMON/ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB10II/NHII(NNNW)
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*
*  Final state common
*
      COMMON/ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)

*
*   Loop over the initial and final state radial functions
*
C      DO I = 1,NSHLI
C        DO J = 2,MFII(I)
C          WRITE(97,*) I,J, PFII(J,I),QFII(J,I),R(J)
C        ENDDO
C      ENDDO

C      DO I = 1,NSHLF
C        DO J = 2,MFFF(I)
C          WRITE(96,*) I,J, PFFF(J,I),QFFF(J,I),R(J)
C        ENDDO
C      ENDDO

      DO I=1,NSHLI
        DO J=1,NSHLF
*
*   Determine the maximum tabulation point for the integrand
*
          MTP=MIN(MFII(I),MFFF(J))
*
*   Tabulate the integrand as required for SUBROUTINE QUAD; the
*   value at the first tabulation point is arbitrary
*
          TA(1)=0.D0
          DO L=2,MTP
            TA(L)=(PFII(L,I)*PFFF(L,J)+
     :      QFII(L,I)*QFFF(L,J))*RP(L)
          ENDDO
*
*   Perform the quadrature
*
          CALL QUAD(RESULT)
          S(I,J)=RESULT
        ENDDO
      ENDDO
*
*   Print out
*
      WRITE(*,*) '********************'
      WRITE(*,*) ' S matrix from GETS'
      WRITE(*,*) '********************'

      CALL WRTMAT(S,NSHLI,NSHLF,NSHLI,NSHLF)
*
      RETURN
      END
