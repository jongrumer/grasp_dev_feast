************************************************************************
*                                                                      *
      SUBROUTINE matrix(dvdfirst)
      IMPLICIT REAL*8          (A-H, O-Z)
*                                                                      *
*   Calls routines to form the Hamiltonian matrix and to diagonalise   *
*   it. The total angular momenta and parity of each  ASF  is found;   *
*   the eigenvectors are normalised so  that the sign of the largest   *
*   element is positive.                                               *
*                                                                      *
*   Call(s) to: [LIB92] ALLOC, DALLOC.                                 *
*               [RSCF92]: SETHAM, HMOUT.                               *
*               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL                       *
*                                                                      *
*                                         Last revision: 24 Dec 1992   *
* Block version by Xinghong He                           07 Aug 1998   *
*                                                                      *
************************************************************************
*
      include 'parameters.def'

*      integer*8 nelmnt
*       first call to dvdson or not
	logical  dvdfirst


CGG      PARAMETER (NNNW = 120)

      !*** Orbital damping and eigenvector damping
      POINTER (pcdamp,cdamp(*))
      COMMON/damp/odamp(NNNW),pcdamp

      LOGICAL ldbpg
      COMMON/debug/ldbpg(5)

      !*** eigenstate indeces
      POINTER (pccmin,iccmin(*))
      COMMON/def7/pccmin,ncmin,ncmax

      !*** eigenvalues and eigenvectors
      POINTER (pneval,eval(*))
      POINTER (pnevec,evec(*))
      COMMON/eigval/eav,pneval
     :      /eigvec/pnevec

      !*** hamiltoniam matrix: elements and indeces
      POINTER (pntemt,emt(*))
      POINTER (piendc,iendc(0:*))
      POINTER (pnirow,irow(*))
      COMMON/hmat/pntemt,piendc,pnirow,nelmnt

      POINTER (pntriq,riqdum)
      COMMON/orb2/ncf,nw,pntriq

      !*** J and parity
      POINTER (piatjp,iatjpo(*))
      POINTER (piaspa,iaspar(*))
      COMMON/syma/piatjp,piaspa

      !*** Block info
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

      COMMON/mcpa/kmaxf

      LOGICAL FIRST
      DATA FIRST /.TRUE./
      SAVE FIRST

      CHARACTER*3 mcplab
      POINTER (pncmvl,cmvl(*))

      COMMON /mpi/ myid, nprocs, ierr
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------

      !IF (myid .EQ. 0) PRINT *, 'matrix ...'
      IF (myid .EQ. 0) PRINT *

*   Allocate memory for CMVL once (the maximum size)
*   Save previous estimate of eigenvectors

      IF (.NOT. FIRST) THEN
         CALL alloc (pncmvl, nvecsiz, 8)
         CALL dcopy(nvecsiz, evec, 1, cmvl, 1)
      ENDIF

*=======================================================================
*   Position the files - MCP files (unit NFILE) for reading 
*   and mixing coefficients file (unit 25) for writing
*=======================================================================

      DO nfile = 30, 32 + kmaxf
         REWIND (nfile)
         IF (nfile .EQ. 30) THEN
            READ (nfile)
            READ (nfile)
            READ (nfile)
         ENDIF
         READ (nfile)
         READ (nfile)
         READ (nfile)
      ENDDO

      ! To put in ncmin and nvecsiz. Values read here are the same
      ! as those from elsewhere (shch as common blocks)
      IF (myid .EQ. 0) THEN
         REWIND (25)
         READ (25)   ! 'G92MIX'
         READ (25) nelec, ncftot, nw, ntmp, ntmp, nblock
         BACKSPACE (25)
         WRITE (25) nelec, ncftot, nw, ncmin, nvecsiz, nblock
      ENDIF

*=======================================================================
*   Do the job block by block
*=======================================================================

               !------------------------------------------------
       			DO 456 jblock = 1, nblock ! block do-loop
               !------------------------------------------------

*=======================================================================
*   Read indeces of non-zero elements from mcp.30 file. Note the
*   format has been changed to lower-triangle-by-rows.
* Length of iendc can be reduced 
*=======================================================================

      READ (30) mcplab, jblockt, ncf
      IF (jblockt .NE. jblock .OR. ncf .NE. ncfblk(jblock))
     &   STOP 'matrx: jblockt .NE. jblock .OR. ncf1 .NE. ncf2'
      READ (30) nelmnt

      CALL alloc (pnirow, nelmnt, 4)
      CALL alloc (pntemt, nelmnt, 8)
      CALL alloc (piendc, ncf + 1, 4)

      DO i = 0, ncf     ! may not be necessary if iendc is ALWAYS used
         iendc(i) = 0   ! the way it is assigned here.
      ENDDO

      !...EMT will be accumulated in setham
      DO i = 1, nelmnt
         emt(i) = 0.D0
      ENDDO

      READ (30) (iendc(i), i = myid + 1, ncf, nprocs), 
     &           (irow(i), i = 1, nelmnt)

      ncfpat = ncfpast(jblock)
      ncminpat = ncminpast(jblock)
      nevecpat = nevecpast(jblock)

*=======================================================================
*   Skip current block if no eigenlaue is required
*=======================================================================

      IF (nevblk(jblock) .EQ. 0) THEN
         DO nfile = 31, 32 + kmaxf
            READ (nfile) mcplab, jblockt, ncft, ncoeff
            IF (jblockt .NE. jblock) STOP 'matrx: jblockt .NE. jblock'
            IF (ncft .NE. ncf) STOP 'matrx: ncft .NE. ncf'

            READ (nfile) lab, ncontr
            DO WHILE (lab .NE. 0 .OR. ncontr .NE. 0)
               READ (nfile) (itmp, itmp, tmp, i = 1, ncontr)
               READ (nfile) lab, ncontr
            ENDDO
         ENDDO

         CALL dalloc (piendc)
         CALL dalloc (pnirow)

         CYCLE

      ENDIF

*=======================================================================
*   Generate the Hamiltonian matrix - average energy is removed here
*=======================================================================

      CALL setham (jblock, myid, nprocs)
*
*   Determine average energy
*
      eav = 0.D0
      DO ir = myid + 1, ncf, nprocs
         eav = eav + emt(iendc(ir))
      ENDDO

      eav = eav / ncf
      eavblk(jblock) = eav

      ! Print Hamiltonian matrix and average energy
      ! hmout is not general
      !call hmout (0, 1, ncf,eav)

      IF (myid .EQ. 0) THEN
         WRITE (*,302) eav
      ENDIF

      ! Subtract the average energy from the diagonal elements
      ! to reduce the condition number of the matrix
      DO i = myid + 1, ncf, nprocs
         idiag = iendc(i)     ! new mode: each row ends in diagonal
         emt(idiag) = emt(idiag) - eav
      ENDDO

*=======================================================================
*   Compute and store eigenpairs
*=======================================================================

      CALL maneig (dvdfirst, ldbpg(3), 
     &          jblock, ncfpat, ncminpat, nevecpat, ncftot)

*=======================================================================
*   Damp and Schmidt orthogonalise eigenvectors for OL calculations
*=======================================================================

      IF (.NOT. FIRST) THEN

         DO 345 J = 1, nevblk(jblock)

            iofset = (j-1) * ncf + nevecpat
            jother = j

            ! cdamp has the original non-block feature
            cdampj = cdamp(j+ncminpat)
            IF (cdampj .EQ. 0.D0) CYCLE                  ! So SURE ???

            omcdaj = 1.D0 - cdampj

           !...Damp eigenvector and determine the new dominant component
  123       amax = 0.D0
            DO i = 1, ncf
               evecij = omcdaj * evec(i+iofset) +
     :                  cdampj * cmvl(i+iofset)
               evec(i+iofset) = evecij
               wa = ABS (evecij)
               IF (wa .GT. amax) THEN
                  amax = wa
                  ia = i
               ENDIF
            ENDDO

            !...compute the normalization factor
            sum = 0.D0
            DO i = 1, ncf
               sum = sum + evec(i+iofset)**2
            ENDDO
            dnfac = 1.D0 / SQRT (sum)

            !...Renormalize and invert as necessary
            IF (evec(ia+iofset) .LT. 0.D0) dnfac = - dnfac
            CALL dscal(ncf, dnfac, evec(iofset+1), 1)

            !...Schmidt orthogonalise
  234       jother = jother - 1
            IF (jother .GE. 1) THEN
               jofset = (jother - 1) * ncf + nevecpat
               ovrlap = ddot (ncf, evec(iofset+1), 1,
     :                             evec(jofset+1), 1)
               IF (ovrlap .NE. 0.D0) THEN                ! So SURE ???
                  omcdaj = 1.D0
                  cdampj = - ovrlap
                  CALL dcopy(ncf, evec(jofset+1), 1,
     :                             cmvl(iofset+1), 1)
                  GOTO 123
               ELSE
                  GOTO 234
               ENDIF
            ENDIF
  345    CONTINUE
      ENDIF

*   Write out the eigenpair information: ASF symmetries, eigenvalues, 
*   and eigenvectors to GRASP92 mixing coefficients File

      IF (nevblk(jblock) .EQ. 0) THEN
         iattmp = 999
         iastmp = 999
      ELSE
         iattmp = iatjpo(ncminpat+1)
         iastmp = iaspar(ncminpat+1)
      ENDIF

      IF (myid .EQ. 0) THEN
         WRITE (25) jblock, ncf, nevblk(jblock), iattmp, iastmp
         WRITE (25) (iccmin(i+ncminpat), i = 1, nevblk(jblock) )
         WRITE (25) eav, (eval(i+ncminpat), i = 1, nevblk(jblock))
         WRITE (25) ((evec(i+(j-1)*ncf+nevecpat),
     &             i = 1, ncf), j = 1, nevblk(jblock))
      ENDIF

      CALL dalloc (pntemt)
      CALL dalloc (piendc)
      CALL dalloc (pnirow)

!----------------------
  456 			CONTINUE
!----------------------
*
*   Deallocate the temporary storage
*
      IF (.NOT. FIRST) THEN
         CALL dalloc (pncmvl)
      ELSE
         FIRST = .FALSE.
      ENDIF

  302 FORMAT (' Average energy = ',1PD18.10,' Hartrees')

      RETURN
      END
