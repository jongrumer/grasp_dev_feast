************************************************************************
*                                                                      *
      SUBROUTINE genmat2 (irestart, nelmnt_a, elsto)
      IMPLICIT REAL*8          (A-H, O-Z)
*
*   Get eav and do writings to the summary file .csum
*   The mpi version (genmat2mpi) also gets nelmnt_a and elsto
*
* Xinghong He 98-06-15
*
************************************************************************

      integer*8 nelmnt_a, nelmnttmp

      COMMON/setham_to_genmat2/CUTOFFtmp,
     &  NCOEItmp, NCOECtmp, NCTEItmp, NCTECtmp, NTPItmp(6), NMCBPtmp, 
     &  NCOREtmp, NVPItmp, NKEItmp, NVINTItmp, NELMNTtmp, NCFtmp

      POINTER (PNEVAL,EVALDUMMY)
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /EIGVAL/EAV,PNEVAL

      DIMENSION NTPI_a(6)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Transfer parameters around and print out on node-0.
!  .Temporary variables (ending with "_a") are used to store the 
!   sums over nodes;
!  .Outputs to stream 24 are done on node-0;
! The following parameters are accumulated in setham which will not 
! contain the correct values in restart mode. And thus is skipped.
!-----------------------------------------------------------------------
      IF (irestart .NE. 0) THEN ! non-restart mode
         WRITE (24,301) CUTOFFtmp
         WRITE (24,302) NCOEItmp
         WRITE (24,303) NCOECtmp
         WRITE (24,304) NCTEItmp
         WRITE (24,305) NCTECtmp
         IF (LTRANS) THEN
            WRITE (24,306) NTPItmp
            WRITE (24,307) NMCBPtmp
            WRITE (24,308) NCOREtmp
         ENDIF
         IF (LVP) WRITE (24,309) NVPItmp
         IF (LNMS) WRITE (24,310) NKEItmp
         IF (LSMS) WRITE (24,311) NVINTItmp
      ELSE
			WRITE (24,*) 'Restart mode --- no report on radial integrals'
      ENDIF !(irestart .NE. 0) ! non-restart mode


!-----------------------------------------------------------------------
! ELSTO, EAV are not only for print-out, but also used later.
! density of the Hamiltonian matrix is only for print-out.
! At this place Hamiltonian matrix elements (EMT on each node) do
!   _not_ contain ELSTO. ELSTO will be added to the total energy
!  later with EAV.
!-----------------------------------------------------------------------

      NELMNT_a = NELMNTtmp
      
      DENSTY = DBLE (NELMNTtmp) / DBLE ((NCFtmp*(NCFtmp+1))/2)

      EAV = EAV/ DBLE (ncftmp) + ELSTO

      WRITE (24,312) NELMNTtmp
      WRITE (24,313) DENSTY
      WRITE (24, *)
      WRITE (24, 300) eav
      WRITE (24, *) eav,elsto,ncftmp

  300 FORMAT ('Average energy = ',1PD19.12,' Hartrees.')
  301 FORMAT ('CUTOFF set to ',1PD17.10)
  302 FORMAT ('Dirac-Coulomb one-e radial integrals:',1I8)
  303 FORMAT ('One-e angular integrals that exceed CUTOFF: ',1I8)
  304 FORMAT ('Coulomb two-e radial integrals: ',1I8)
  305 FORMAT ('Two-e angular integrals that exceed CUTOFF: ',1I8)
  306 FORMAT ('Transverse two-e radial integrals: '/6I8)
  307 FORMAT ('MCBP coefficients that exceed CUTOFF: ',1I8)
  308 FORMAT ('Core coefficients that exceed CUTOFF: ',1I8)
  309 FORMAT ('Vacuum polarisation integrals: ',1I8)
  310 FORMAT ('Kinetic energy integrals: ',1I8)
  311 FORMAT ('Vinti integrals: ',1I8)
  312 FORMAT ('Elements that exceed CUTOFF in the lower',
     :        ' triangle of the H matrix: ',1I8)
  313 FORMAT ('Density of the H(amiltonian) matrix: ',1PD22.15)

      RETURN
      END
