************************************************************************
*                                                                      *
      SUBROUTINE maneigmpi (lprint, jblock, 
     &                   ncfpat, ncminpat, nevecpat, ncftot)
*                                                                      *
*   This module  manages the  operation of the  eigensolvers and the   *
*   storage of the eigenpairs.  There are two principal branches:      *
*                                                                      *
*      (1) Matrix of order 1: the trivial case                         *
*      (2) Matrix of order greater than 1: eigenpairs are found        *
*          using DVDSON; this involves up to three steps:              *
*                    (a) The matrix is analysed to determine its       *
*                        block structure (only irreducibe matrices     *
*                        are correctly treated by DVDSON)              *
*                    (b) Eigenpairs are extracted for each block       *
*                    (c) The appropriate eigenpairs are selected and   *
*                        stored                                        *
*                                                                      *
*   We  assume that  the sparse representation  of the matrix  is in   *
*   core.                                                              *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, ISPAR, ITJPO, RALLOC.          *
*               [RSCF92]: POSNFL, SPICMVmpi.                      *
*               [DVDSON]: DVDSON                                       *
*               [AUXBLAS]: DINIT/SINIT                                 *
*               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL, DSWAP/SSWAP          *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
*   Modified by Xinghong He               Last revision: 17 Aug 1998   *
*                                                                      *
************************************************************************

      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      POINTER (PNTRIQ,RIQDUM)
      LOGICAL HIEND,LPRINT

      EXTERNAL spicmvmpi

      POINTER (PNWORK,WORK(1))
      POINTER (PIWORK,IWORK(1))
      POINTER (PNDIAG,DIAG(1))
      POINTER (PJWORK,JWORK(1))

      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PNEVAL,EVAL(1))
      POINTER (PNEVEC,EVEC(1))
      POINTER (PNTEMT,EMT(1))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(1))
      POINTER (PIATJP,IATJPO(1))
      POINTER (PIASPA,IASPAR(1))

      COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /EIGVAL/EAV,PNEVAL
     :      /EIGVEC/PNEVEC
     :      /HMAT/PNTEMT,PIENDC,PNIROW,NELMNT
     :      /ORB2/NCF,NW,PNTRIQ
     :      /SYMA/PIATJP,PIASPA
     :      /WCHBLK/JBLOCKK
     :      /WHERE/IMCDF,NREC

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(1))
      POINTER (pncmaxblk, ncmaxblk(1))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (ptmp, atmp(1))

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      !IF (myid .EQ. 0) PRINT *, 'maneigmpi ...'

*      ...SPICMVmpi needs this COMMON /WCHBLK/JBLOCKK
      jblockk = jblock
*
*=======================================================================
*   Trivial case
*=======================================================================
      IF (ncf .EQ. 1) THEN
         eval(ncminpat+1) = 0.D0
         evec(nevecpat+1) = 1.D0
         GOTO 123    ! Don't like big ELSE
      ENDIF
*
*=======================================================================
*   Non-trivial case - Use Davidson eigensolver
*=======================================================================
*
*   Allocate storage for workspace; see the header of DVDSON for
*   the expression below; the value of LIM can be reduced to NVECT
*   plus a smaller number if storage is severely constrained
*
      nvect = ncmaxblk(jblock)
      lim   = MIN (ncf, nvect + 200)
!     lim   = MIN (ncf, nvect + 20)
      lwork = 2*ncf*lim + lim*lim*2 + 11*lim + nvect

      CALL alloc (pnwork, lwork, 8)

      !...At most 14 ? restriction removed xhh 98-05-19
      !nvex = MIN (nvect,ncfblk(jblock),14)
      nvex = MIN (nvect, ncfblk(jblock))
      niv = nvex
      maxitr = MIN (nvect*100, ncf)
      n1000 =1000
*
*   Initial estimates for eigenvectors
*
	   CALL iniestmpi (n1000, ncf, niv, work, emt, iendc, irow)

* iniest looks for eigenvectors of n1000*n1000 matrix so there 
* is no need to call dvdson if block size <= n1000

	   IF (ncf .GT. n1000 ) THEN
         IF (myid .EQ.0) WRITE (*,*) 'Calling dvdson!!!', maxitr,nvect
 
         ! Call Davidson eigensolver

         mblock = 1
         ilow = 1
         ihigh = nvex
         liwork = 6*lim + nvect
         crite = 1.0D-17
         critc = 1.0D-08
         critr = 1.0D-08
         ortho = MAX (1D-8, critr)
*
*   Store the diagonals in a separate array and make it global
*
         CALL alloc (pndiag, ncf, 8)
         CALL alloc (ptmp, ncf, 8)
         DO i = 1, ncf
            atmp(i) = 0.D0
            diag(i) = 0.D0    ! this one may not be necessary
         ENDDO

         DO ic = myid + 1, ncf, nprocs
            atmp(ic) = emt(iendc(ic))
         ENDDO
         CALL MPI_Allreduce (atmp, diag, ncf, MPI_DOUBLE_PRECISION,
     &                         MPI_SUM, MPI_COMM_WORLD, ierr)
         CALL dalloc (ptmp)

         CALL alloc (piwork, liwork, 4)
         CALL alloc (pjwork, lim, 4)
         if (ncf.gt.1000) then
         CALL gdvd (spicmvmpi,ncf,lim,diag,ilow,ihigh,
     :            jwork,niv,mblock,crite,critc, critr,ortho,maxitr,
     :            work,lwork,iwork,liwork,hiend,nloops,
     :            nmv,ierr)
         end if
         CALL dalloc (pndiag)
         CALL dalloc (piwork)
         CALL dalloc (pjwork)

         IF (myid .EQ. 0) THEN
            WRITE (*,301) nloops, nmv
            IF (ierr .NE. 0) THEN
               WRITE (*,302) ierr
            ENDIF
         ENDIF
      ENDIF
*
*   Pick up the eigen pairs and store in EVAL and EVEC
*
      nend = ncf * nvex
      DO j = 1, nevblk(jblock)
         eval(ncminpat+j) = work( nend + iccmin(j+ncminpat) )
         CALL dcopy (ncf, work( ncf*(iccmin(j+ncminpat)-1) + 1 ), 1,
     &                    evec( nevecpat + ncf*(j-1) + 1)    ,    1)
      ENDDO
*
*   Deallocate storage
*
      CALL dalloc (pnwork)

  123 CONTINUE
*
*   Clean up eigenvectors; determine their J/P values
*
      DO 17 Jstate = 1, nevblk(jblock)
*
*   Find the dominant component of each eigenvector
*
         iofset = nevecpat + ncf * (Jstate - 1)

         amax = 0.D0
         DO i = 1, ncf
            wa = ABS (evec(i+iofset))
            IF (wa .GT. amax) THEN
               amax = wa
               ia = i
            ENDIF
         ENDDO
*
*   Find the angular momentum and parity of the dominant component
*
         iatjpo(Jstate+ncminpat) = itjpo (ia + ncfpat)
         iaspar(Jstate+ncminpat) = ispar (ia + ncfpat)
*
*   Redefine eigenvectors so that the dominant component
*   is positive
*
         IF (evec(ia+iofset) .LT. 0.D0) THEN
            dnfac = -1.D0
            CALL dscal (ncf, dnfac, evec(iofset+1), 1)
         ENDIF
!===============================================================

   17 CONTINUE

  301 FORMAT ('DVDSON: ',1I3,' loops; ',
     :                   1I3,' matrix-vector multiplies.')
  302 FORMAT (' Returned from DVDSON with IERR = ',1I4)
  303 FORMAT (/' ***** WARNING *****'
     :       //' The angular momentum and parity of level ',1I2,
     :         ' have changed:'
     :        /' Last iteration: (2J+1) = ',1I2,', parity = ',1I2,';'
     :        /' this iteration: (2J+1) = ',1I2,', parity = ',1I2,'.')

      RETURN
      END
