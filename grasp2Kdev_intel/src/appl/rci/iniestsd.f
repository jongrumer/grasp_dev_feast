************************************************************************
*
      SUBROUTINE INIESTSD (nmax, ncf, myid, nprocs, 
     &    NIV, BASIS, IMCDF, EAV)
*
*    Routine for providing initial estimates from upper left corner
*	  of the matrix. It is exact (not estimates) if ncf <= nmax which
*    is currently set to be 400 in the calling routine.
*
*    Matrix is sparse and on the disk
*
*   Block version by Xinghong He            Last revision: 14 Dec 1998
*
************************************************************************
      IMPLICIT REAL*8          (A-H, O-Z)
      DIMENSION basis(*)
      ! Locals
      POINTER (iqap,ap(*)),(iqeig,eigval(*)),(iqvec,vec(*))
      POINTER (iqwork,work(*)),(iqiwork,iwork(*)),(iqif,IFAIL(*))
      POINTER (phmx, hmx(*)), (pirow, irow(*))
!-----------------------------------------------------------------------
      NS = min (nmax, ncf)

      CALL alloc (iqap, (NS*(NS+1))/2, 8)
      CALL dinit ((NS*(NS+1))/2, 0.d0, ap, 1)

***** separate upper left block of size NS*NS

      CALL alloc (phmx, ncf, 8)
      CALL alloc (pirow, ncf, 4)
      READ (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
      IF (myid .EQ. 0) PRINT *, 'iniestsd ...........'
      IF (ncf .NE. ncfdum .OR.  myid .NE. myiddum
     &              .OR. nprocsdum .NE. nprocs) 
     &   STOP 'iniestsd: ncf read wrong'

*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.

      DO j = myid + 1, ns, nprocs
         joff = (j*(j-1))/2
         READ (IMCDF) NELC,ELSTO,(HMX(IR),IR=1,NELC),
     :                          (IROW(IR),IR=1,NELC) 
         HMX(NELC) = HMX(NELC) - EAV ! Shift the diagonal
         DO ir = 1, nelc
            ap(irow(ir) + joff) = hmx(ir)
         ENDDO
      ENDDO

! Let each node have a complete copy of ap

!      CALL gdsummpi (ap, (NS*(NS+1))/2)

! To be in step with other cases, go through the whole block.
!
! This is not necessary since currently the file pointer is moved
! to the absolute position of the .res files which is always counted
! from the begining of the .res files of each node. Besides, the
! following segment seems not working properly for the last block.
! Xinghong He 98-12-14

      !mylast = j - nprocs
      !DO j = mylast, ncf, nprocs
      !   READ (imcdf)
      !ENDDO

      CALL dalloc (phmx)
      CALL dalloc (pirow)

      CALL alloc (iqeig,NS,8)
      CALL alloc (iqvec,NS*NIV,8)
      CALL alloc (iqwork,8*NS,8)
      CALL alloc (iqiwork,5*NS,4)
      CALL alloc (iqif,NS,4)

      CALL DSPEVX ('Vectors also','In a range','Upper triangular',
     :          NS,AP,-1.,-1.,1,NIV,0.d0,
     :          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      IERR = -ABS (INFO)

*******************************************************************

*       ..Build the Basis. 

      CALL DINIT (ncf*NIV, 0.D0, BASIS, 1)
*       ...scatter the vectors
      DO J = 1, NIV
         CALL dcopy (ns, vec(ns*(j-1)+1),1, basis(ncf*(j-1)+1), 1)
      ENDDO
      CALL dcopy (NIV, EIGVAL,1,BASIS(NIV*ncf+1),1)

      CALL dalloc (iqap)
      CALL dalloc (iqeig)
      CALL dalloc (iqvec)
      CALL dalloc (iqwork)
      CALL dalloc (iqiwork)
      CALL dalloc (iqif)

      RETURN
      END
