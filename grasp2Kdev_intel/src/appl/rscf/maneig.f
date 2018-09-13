************************************************************************
*                                                                      *
      SUBROUTINE maneig (dvdfirst, lprint, jblock, 
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
*               [RSCF92]: POSNFL, SPICMV2.
*               [DVDSON]: DVDSON                                       *
*               [AUXBLAS]: DINIT/SINIT                                 *
*               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL, DSWAP/SSWAP          *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Sep 1993   *
*   Modified by Xinghong He               Last revision: 17 Aug 1998   *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8          (A-H, O-Z)

*      integer*8 nelmnt


      POINTER (PNTRIQ,RIQDUM)
CFF	... add variable first
      LOGICAL HIEND,LPRINT, FIRST, dvdfirst

      EXTERNAL spicmv2

c     POINTER (PNWORK,WORK(29999))
      POINTER (PNWORK,WORK(*))
      POINTER (PIWORK,IWORK(*))
c     POINTER (PNDIAG,DIAG(19999))
      POINTER (PNDIAG,DIAG(*))
      POINTER (PJWORK,JWORK(*))

      POINTER (PCCMIN,ICCMIN(*))
c     POINTER (PNEVAL,EVAL(19999))
      POINTER (PNEVAL,EVAL(*))
c     POINTER (PNEVEC,EVEC(1999))
      POINTER (PNEVEC,EVEC(*))
c     POINTER (PNTEMT,EMT(1999))
      POINTER (PNTEMT,EMT(*))
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))
      POINTER (PIATJP,IATJPO(*))
      POINTER (PIASPA,IASPAR(*))

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

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      !PRINT *, 'maneig ...'

*      ...spicmv2 needs this COMMON /WCHBLK/JBLOCKK
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
!     lim   = MIN (ncf, nvect + 20)
      lim   = MIN (ncf, 2*nvect + 40)
      lwork = 2*ncf*lim + lim*lim*2 + 11*lim + nvect

      CALL alloc (pnwork, lwork, 8)

      !...At most 14 ? restriction removed xhh 98-05-19
      !nvex = MIN (nvect,ncfblk(jblock),14)
      nvex = MIN (nvect, ncfblk(jblock))
      niv = nvex
!     maxitr = MIN (nvect*100, ncf)
      maxitr = MIN (nvect*200, ncf)
!     n1000 = 2000
      n1000 = 1000
*
*   Initial estimates for eigenvectors
*
CFF
      if (dvdfirst .or. (ncf .LE. n1000) ) then
	   CALL iniest2 (n1000, ncf, niv, work, emt, iendc, irow)
      else
CFF 	   .. use current estimates
      nend = ncf * nvex
         DO j = 1, nevblk(jblock)
            work( nend + iccmin(j+ncminpat) ) = eval(ncminpat+j)
            CALL dcopy (ncf, evec( nevecpat + ncf*(j-1) + 1)    ,    1,
     &                    work( ncf*(iccmin(j+ncminpat)-1) + 1 ), 1)
         ENDDO
      ENDIF


* iniest2 looks for eigenvectors of n1000*n1000 matrix so there 
* is no need to call dvdson if block size <= n1000

	   IF (ncf .GT. n1000 ) THEN
         WRITE (*,*) 'Calling dvdson!!!', maxitr
 
         ! Call Davidson eigensolver

         mblock = 1
         ilow = 1
         ihigh = nvex
         liwork = 6*lim + nvect
         crite = 1.0D-17
!         critc = 1.0D-08
!         critr = 1.0D-08
!         ortho = MAX (1D-8, critr)
         critc = 1.0D-09
         critr = 1.0D-09
         ortho = MAX (1D-9, critr)
*
*   Store the diagonals in a separate array and make it global
*
         CALL alloc (pndiag, ncf, 8)

         DO ic = myid + 1, ncf, nprocs
            diag(ic) = emt(iendc(ic))
         ENDDO

         CALL alloc (piwork, liwork, 4)
         CALL alloc (pjwork, lim, 4)
         CALL gdvd (spicmv2,ncf,lim,diag,ilow,ihigh,
     :            jwork,niv,mblock,crite,critc, critr,ortho,maxitr,
     :            work,lwork,iwork,liwork,hiend,nloops,
     :            nmv,ierr)

         CALL dalloc (pndiag)
         CALL dalloc (piwork)
         CALL dalloc (pjwork)

         WRITE (*,301) nloops, nmv
         IF (ierr .NE. 0) THEN
            WRITE (*,302) ierr
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
c     print *, ncminpat,(eval(ncminpat+j),j=1,nevblk(jblock)),
c    1 'zou,from maneig'
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
      write(507,*) evec(1:ncf*nevblk(jblock))
      write(507,*) '+++++++++++++++++++++++'

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
