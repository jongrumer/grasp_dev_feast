************************************************************************
*                                                                      *
      SUBROUTINE newco(SUM)
*                                                                      *
*   This routine computes the level weights, the generalized occupa-   *
*   tion numbers,  and average energy for  EOL calculations;  this     *
*   information and the eigenvectors are then printed out.             *
*                                                                      *
*   Call(s) to: [LIB92]: IQ.                                           *
*               [RSCF92]: CSFWGT, DSUBRS.                              *
*                                         Last revision: 24 Dec 1992   *
*   Block version by Xinghong He          Last revision: 05 Aug 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)

      POINTER (PIASPA,IASPADUMMY)
      LOGICAL EOL,LDBPG
*
      POINTER (PNTRWT,WT(*))
      POINTER (PWEIGH,WEIGHT(*))
      POINTER (PCCMIN,ICCMIN(*))
      POINTER (PNEVAL,EVAL(*))
      POINTER (PNEVEC,EVEC(*))
      POINTER (PIATJP,IATJPO(*))
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/DEBUGG/LDBPG(5)
     :      /DEF5/PNTRWT,PWEIGH
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /SCF1/UCF(NNNW)
     :      /SYMA/PIATJP,PIASPA
      COMMON/iounit/istdi,istdo,istde

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (pncfpast, ncfpast(*))
      POINTER (pncminpast, ncminpast(*))
      POINTER (pnevecpast, nevecpast(*))
      COMMON/pos/pncfpast,pncminpast,pnevecpast,ncftot,nvecsiz

      POINTER (peavblk, eavblk(*))
      COMMON/peav/peavblk

      POINTER (pidxblk, idxblk(*))
      COMMON/blkidx/pidxblk

      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      !PRINT *, 'NEWCO ...'
! pfeast hack, local nprocss, myidd
      nprocss = 1
      myidd = 0
*
*   Compute weighting factors
*
      SUM = 0.D0
      DO J = 1, NCMIN
         WEITJ = WEIGHT(J)
         IF (WEITJ .EQ. -2.D0) THEN
            WT(J) = 1.D0
         ELSEIF (WEITJ .EQ. -1.D0) THEN
            WT(J) = IATJPO(J)
         ELSE
            WT(J) = WEIGHT(J)
         ENDIF
         SUM = SUM + WT(J)
      ENDDO

      DO J = 1, NCMIN
         WT(J) = WT(J) / SUM
      ENDDO
*
*   Compute generalised occupation numbers
*   <----- Distributed ----->
*
      EOL = .TRUE.

      DO J = 1, NW
         SUM = 0.D0
         DO nb = 1, nblock
            DO I = myidd + 1, NCFblk(nb), nprocss
               SUM = SUM + DSUBRS (EOL, I, I, nb) * 
     &                     IQ (J, I + ncfpast(nb))
            ENDDO
         ENDDO
         UCF(J) = SUM
      ENDDO
*
*   Write out level energies and weights
*
      if(myid==0) WRITE (*,300)
      SUM = 0.D0
      noff = 0
      DO Jall = 1, NCMIN
         nb = idxblk(jall)		! Block number of this state
         EE = EAVblk(nb) + EVAL(Jall)
         if(myid == 0) then
         WRITE (*,301) iccmin(jall), EE, WT(Jall)
         IF (LDBPG(5)) THEN
            WRITE (99,*) jall, nb, NCFblk(nb), nevecpast(nb)
            WRITE (99,302)
            WRITE (99,303) (EVEC(I+noff), I = 1, NCFblk(nb))
            noff = noff + NCFblk(nb)
         ENDIF
         endif
         SUM = SUM + WT(Jall) * EE
      ENDDO

      CALL CSFWGT (.TRUE.)
*
*   Write out average energy
*
      IF (NCMIN .GT. 1 .and. myid == 0) WRITE (*,304) SUM
*
*   Write out generalized occupation numbers
*
      if(myid==0) then
      WRITE (*,305)
      WRITE (*,303) (UCF(I),I = 1,NW)
      endif

  300 FORMAT (/'Optimise on the following level(s):'/)
  301 FORMAT ('Level ',1I2,4X,'Energy = ',1P,1D19.12,
     :                      4X,'Weight = ',  1D12.5)
  302 FORMAT (/'Configuration mixing coefficients:')
  303 FORMAT (1X,1P,6D12.4)
  304 FORMAT (/'Weighted average energy of these levels = ',1PD18.10)
  305 FORMAT (/'Generalised occupation numbers:'/)

      RETURN
      END
