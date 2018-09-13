************************************************************************
*                                                                      *
      SUBROUTINE genmat2 (irestart, nelmnt_a, elsto)
      IMPLICIT REAL*8          (A-H, O-Z)
*
*   Get eav and do writings to the summary file .csum
*   The mpi version (here) also gets nelmnt_a and elsto
*
* Xinghong He 98-06-15
*
************************************************************************
      COMMON/setham_to_genmat2/CUTOFFtmp,
     &  NCOEItmp, NCOECtmp, NCTEItmp, NCTECtmp, NTPItmp(6), NMCBPtmp, 
     &  NCOREtmp, NVPItmp, NKEItmp, NVINTItmp, NELMNTtmp, NCFtmp

      POINTER (PNEVAL,EVALDUMMY)
      LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
      COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
     :      /EIGVAL/EAV,PNEVAL

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr

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
      CALL MPI_Reduce (NCOEItmp, NCOEI_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 
      CALL MPI_Reduce (NCOECtmp, NCOEC_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 
      CALL MPI_Reduce (NCTEItmp, NCTEI_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 
      CALL MPI_Reduce (NCTECtmp, NCTEC_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 

      IF (LTRANS) THEN
         CALL MPI_Reduce (NTPItmp, NTPI_a, 6, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 
         CALL MPI_Reduce (NMCBPtmp, NMCBP_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 
         CALL MPI_Reduce (NCOREtmp, NCORE_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 
      ENDIF

      IF (LVP)
     &   CALL MPI_Reduce (NVPItmp, NVPI_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 

      IF (LNMS)
     &   CALL MPI_Reduce (NKEItmp, NKEI_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &                                 MPI_COMM_WORLD, ierr) 

      IF (LSMS)
     &   CALL MPI_Reduce (NVINTItmp, NVINTI_a, 1, MPI_INTEGER, MPI_SUM,
     &                                 0, MPI_COMM_WORLD, ierr) 

      IF (myid .EQ. 0) THEN
         WRITE (24,301) CUTOFFtmp
         WRITE (24,302) NCOEI_a
         WRITE (24,303) NCOEC_a
         WRITE (24,304) NCTEI_a
         WRITE (24,305) NCTEC_a
         IF (LTRANS) THEN
            WRITE (24,306) NTPI_a
            WRITE (24,307) NMCBP_a
            WRITE (24,308) NCORE_a
         ENDIF
         IF (LVP) WRITE (24,309) NVPI_a
         IF (LNMS) WRITE (24,310) NKEI_a
         IF (LSMS) WRITE (24,311) NVINTI_a
      ENDIF
            ELSE
      IF (myid .EQ. 0) THEN
         WRITE (24,*) 'Restart mode --- no report on radial integrals'
      ENDIF
            ENDIF !(irestart .NE. 0) ! non-restart mode


!-----------------------------------------------------------------------
! ELSTO, EAV are not only for print-out, but also used later.
! density of the Hamiltonian matrix is only for print-out.
! Make ELSTO available to all nodes. In this mpi version, 
!     possibly non-zero ELSTO resides only on node-0.
! Make EAV a global average (including the ELSTO).
! At this place Hamiltonian matrix elements (EMT on each node) do
!   _not_ contain ELSTO. ELSTO will be added to the total energy
!  later with EAV.
!-----------------------------------------------------------------------

      CALL MPI_Allreduce (NELMNTtmp, NELMNT_a, 1, MPI_INTEGER, 
     &                         MPI_SUM, MPI_COMM_WORLD, ierr) 
c     DENSTY = DBLE (NELMNT_a) / DBLE ((NCFtmp*(NCFtmp+1))/2)
      DENSTY = DBLE (NELMNT_a) / (DBLE(NCFtmp)*DBLE(NCFtmp+1)/2.)

      CALL MPI_Bcast (ELSTO, 1, MPI_DOUBLE_PRECISION, 0, 
     &                      MPI_COMM_WORLD, ierr)
      CALL MPI_Allreduce (EAV, EAV_a, 1, MPI_DOUBLE_PRECISION, 
     &                        MPI_SUM, MPI_COMM_WORLD, ierr) 
      EAV = EAV_a / DBLE (ncftmp) + ELSTO

      IF (myid .EQ. 0) THEN
         WRITE (24,312) NELMNT_a
         WRITE (24,313) DENSTY
         WRITE (24, *)
         WRITE (24, 300) eav
      ENDIF

  300 FORMAT ('Average energy = ',1PD19.12,' Hartrees.')
  301 FORMAT ('CUTOFF set to ',1PD17.10)
  302 FORMAT ('Dirac-Coulomb one-e radial integrals:',1I8)
  303 FORMAT ('One-e angular integrals that exceed CUTOFF: ',1I8)
  304 FORMAT ('Coulomb two-e radial integrals: ',1I8)
  305 FORMAT ('Two-e angular integrals that exceed CUTOFF: ',1I11)
  306 FORMAT ('Transverse two-e radial integrals: '/6I8)
  307 FORMAT ('MCBP coefficients that exceed CUTOFF: ',1I8)
  308 FORMAT ('Core coefficients that exceed CUTOFF: ',1I8)
  309 FORMAT ('Vacuum polarisation integrals: ',1I8)
  310 FORMAT ('Kinetic energy integrals: ',1I8)
  311 FORMAT ('Vinti integrals: ',1I8)
  312 FORMAT ('Elements that exceed CUTOFF in the lower',
     :        ' triangle of the H matrix: ',1I11)
  313 FORMAT ('Density of the H(amiltonian) matrix: ',1PD22.15)

      RETURN
      END
