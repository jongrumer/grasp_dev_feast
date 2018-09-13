************************************************************************
*                                                                      *
      SUBROUTINE BRKT
*                                                                      *
*   This subroutine calculates the initial and final state             *
*   radial overlap matrix                                              *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      include 'parameters.def'
cff      PARAMETER (NNNP = 590) 
cff      PARAMETER (NNN1 = NNNP+10)
cff      PARAMETER (NNNW = 120)
*
      DIMENSION BRAKET(NNNW,NNNW)
      CHARACTER*2 NHII,NHFF

      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /TATB/TA(NNN1),TB(NNN1),MTP

      POINTER (PNTRPFII,PFII(NNNP,1)),(PNTRQFII,QFII(NNNP,1))

      POINTER (PNTRPFFF,PFFF(NNNP,1)),(PNTRQFFF,QFFF(NNNP,1))
*
*  Initial state common
*
      COMMON/ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB10II/NHII(NNNW)
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*
*  Final state common
*
      COMMON/ORB2FF/NCFFF,NWFF
     :      /ORB4FF/NPFF(NNNW),NAKFF(NNNW)
     :      /ORB10FF/NHFF(NNNW)
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)

*
      DO I=1,NWII
        DO J=1,NWFF
          IF (NAKII(I) .EQ. NAKFF(J)) THEN
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
     :        QFII(L,I)*QFFF(L,J))*RP(L)
            ENDDO 
*
*   Perform the quadrature
*
            CALL QUAD(RESULT)
            
            BRAKET(I,J)=RESULT

            WRITE(*,9) '<',NPII(I),NHII(I),'|',NPFF(J),NHFF(J),'> =',
     :      BRAKET(I,J)
          ENDIF
        ENDDO
      ENDDO
*
    9 FORMAT(A,I2,A,A,I2,A,A,E20.13)

      RETURN
      END
