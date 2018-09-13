************************************************************************
      SUBROUTINE improv (EOL, J, lsort,DAMPMX)
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL EOL, lsort
      INTEGER J

*   The difference from the mpi version is that it calls serial
*   version subroutines (setlag, cofpot, matrix, newco).
*
*   Improve the orbital J.                                             *
*                                                                      *
*   Call(s) to: [RSCF92]: DACON, DAMPCK, DAMPOR, LAGCON, matrix,
*                         newco, ORTHOR, ROTATE, SETCOF, setlag,
*                         SOLVE, XPOT, YPOT.                           *
*               [LIB92]: ORTHSC, QUAD.                                 *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
*   Modified by Xinghong He                 Last update: 05 Aug 1988   *
*                                                                      *
************************************************************************

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PCDAMP,CDAMPDUM)
      POINTER (PNTRIQ,RIQDUM)
      LOGICAL FAIL,FIRST,ORTHST
      CHARACTER*2 NH

      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
      pointer(PNTRDA,da(NDDIM)),
     :        (PNTNDA,NDADUM), (PNTNXA,NXADUM),
     :        (PNTNYA,NYADUM), (PNTRDA,RDADUM),
     :        (PNTRXA,RXADUM), (PNTRYA,RYADUM) 
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /ORB1/E(NNNW),GAMA(NNNW),PED(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORTHCT/ORTHST
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      PARAMETER (P2    = 2.0D-01,
     :           P005  = 5.0D-03,
     :           P0001 = 1.0D-04)
*
*   C Froese Fischer's IPR and ED1 parameter
*
      DATA IPR /0/
      DATA ED1 /0.D0/
      DATA FIRST /.FALSE./

      LOGICAL lcorre
      COMMON/corre/lcorre(NNNW)
      SAVE /corre/

      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      GAMAJ = GAMA(J)
*
*   C Froese Fischer's parameters IPR, ED1, ED2 are set and
*   used in this routine and in DAMPCK
*
    1 ED2 = E(J)
      ED1 = PED(J)
*
*   Set up the exchange potential and arrays XU, XV as appropriate
*
*   Set coefficients for YPOT, XPOT, DACON
*   Compute direct potential, exchange potential
*   Add in Lagrange-multiplier contribution
*   Add in derivative-terms contribution
*
      npts = N
      CALL cofpot (EOL, J, npts)
*
*   Calculate deferred corrections
*
      CALL DEFCOR (J)
*
*   Solve the Dirac equation
*
      INV = 0
      CALL SOLVE (J, FAIL, INV, JP, NNP)
*
*   Upon failure issue message; take corrective action if possible
*
      IF (FAIL) THEN
         IF (myid .EQ. 0) WRITE (*,300) NP(J),NH(J),METHOD(J)
         IF (METHOD(J) .NE. 2) THEN
            METHOD(J) = 2
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
            CALL ORTHSC
            IF (EOL) THEN
               CALL matrix
               CALL newco(WTAEV)
            ENDIF
            CALL setlag (EOL)
            GOTO 1
         ELSE
            IF (myid .EQ. 0) WRITE (*,301)
            !CALL TIMER (0)
            STOP
         ENDIF
      ENDIF
*
*   Compute norm of radial function
*
      TA(1) = 0.D0
      DO I = 2, MTP0
         TA(I) = (P(I)**2 + Q(I)**2) * RP(I)
      ENDDO
      MTP = MTP0

      CALL QUAD (DNORM)

!   Determine self-consistency [multiplied by SQRT(UCF(J))]

      CALL CONSIS (J)
*
*   Normalize
*
      DNFAC = 1.D0 / SQRT (DNORM)
      P0 = P0 * DNFAC
      DO I = 1, MTP0
         P(I) = P(I) * DNFAC
         Q(I) = Q(I) * DNFAC
      ENDDO
*
*   Check if different method should be used or if improvement
*   count should be reduced
*
      DEL1 = ABS (1.D0 - ED2 / E(J))
      IF (METHOD(J) .EQ. 1) THEN
         DEL2 = MAX (ABS (1.D0 - SQRT (DNORM)),
     :               ABS (DNFAC - 1.D0))
         IF ((DEL1 .LT. P005) .AND. (DEL2 .GT. P2)) THEN
            METHOD(J) = 2
            GOTO 1
         ENDIF
      ELSE
         IF ((DEL1 .LT. P0001) .AND. (NSIC .GT. 1)) NSIC = NSIC-1
      ENDIF
*
*   Damp the orbital --- if not converged
*
       IF (scnsty(J) .GT. ACCY) THEN
          CALL DAMPCK (IPR, J, ED1, ED2)
          odampj = ABS (odamp(j))
       ELSE
          odampj = 0.D0    ! take the whole new orbital
       ENDIF
      CALL DAMPOR (J, INV, odampj)

!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.

      IF (ORTHST) THEN
         !CALL orthor (J, inv)
         nwww = nw
         CALL orthy (nwww, J, lsort)
      ENDIF
*
*   Print details of iteration
*    
      IF (myid .EQ. 0)
     &   WRITE (*,302) NP(J),NH(J),E(J),METHOD(J),PZ(J),SCNSTY(J),
*    &              DNORM-1,ODAMPJ,JP,MF(J),INV,NNP
     &              SQRT(DNORM)-1,ODAMPJ,JP,MF(J),INV,NNP
      DAMPMX = MAX(DAMPMX,ABS(ODAMPJ))

  300 FORMAT (/' Failure; equation for orbital ',1I2,1A2,
     :         ' could not be solved using method ',1I1)
  301 FORMAT (//' ****** Error in SUBROUTINE IMPROV ******'
     :          /' Convergence not obtained'/)
cb
cb 80-col
cb
c 302 FORMAT (1X,1I2,1A2,1P,1D18.9,1x,1I2,D11.3,1D10.2,1D10.2,
  302 FORMAT (1X,1I2,1A2,1P,1D16.7,1x,1I2,D11.3,1D10.2,1D10.2,
c    :        0P,F7.4,1x,1I3,1x,1I3,1x,1I2,2x,1I2)
c    :        0P,F6.3,1x,1I4,1x,1I4,1x,1I2,1x,1I2)
     :        0P,F6.3,1x,1I5,1x,1I5,1x,1I2,1x,1I2)

      RETURN
      END
