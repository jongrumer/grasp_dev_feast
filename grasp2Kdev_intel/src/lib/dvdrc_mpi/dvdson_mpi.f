		SUBROUTINE DVDSON(IRC,IREV,N,LIM,NOC,
     :			 ILOW,IHIGH,ISELEC,NIV,MBLOCK,
     :			 CRITE,CRITC,CRITR,MAXITER, 
     :			 WORK,IWRSZ,IWORK,IIWSZ,
     :			 HIEND,NLOOPS,IERR) 
*=======================================================================
*
*	Author: Andreas Stathopoulos, Charlotte F. Fischer
*	
*	Computer Science Department 
*	Vanderbilt University
*	Nashville, TN 37212
*	andreas@vuse.vanderbilt.edu
*	cff@vuse.vanderbilt.edu
*
*	June 1995
*
*	Copyright (c) by Andreas Stathopoulos and Charlotte F. Fischer
*
*	References: A Stathopoulos and C F Fischer,
*		    Comp. Phys. Commun. 79 (1994) 268-290.
*
*                   A Stathopoulos, Y. Saad and C F Fischer,
*                   J. Computational and applied Mathematics (in print).
*
*	DVDSON is a Fortran77 program that finds a few selected 
*	eigenvalues and their eigenvectors at either end of spectrum of 
*	a large, symmetric (and usually sparse) matrix, denoted as A. 
*	The matrix A is only referenced externally through the user 
*	supplied routine which implements a block matrix-vector 
*	operation(see below) and the user supplied preconditioning
*       routine (also perfomed externally). Either the range of the 
*       eigenvalues wanted or an array of the indices of selected ones 
*       can be specified. 
*	DVDSON is a front-end routine for setting up arrays, and initial
*	guess. It also performs detailed error checking.
*	DVDRVR is the driver routine that implements a version of the 
*	Davidson algorithm. The characteristics of this version are:
*	 o  Use of REVERSE COMMUNICATION to perform the preconditioning
*	    and the matrix vector multiply.
*	 o  All arrays used by the program are stored in MEMORY.
*   	 o  BLOCK method (many vectors may be targeted per iteration.)
* 	 o  Eigenvectors are targeted in an optimum way without
*	    the need to compute all unconverged residuals,
*	 o  Uses two Modified Gram-Schmidt passes for orthogonality loss.
*	 o  Finds HIGHEST eigenpairs by using the negative of the A.
*	 o  Finds SELECTED eigenpairs specified by the user.
*	 o  If not enough initial vectors are available (NIV<NUME)
*           it builds the first NUME-NIV Lanczos from those given.
*	 o  The user can provide STOPPING CRITERIA for eigenvalues, 
*	    and residuals and coefficients and he  can CONTROL block size.
*	 o  On exit INFORMATION is given about the convergence status 
*	    of eigenpairs.
*
*	The program consists of the following routines:
*	DVDSON, INITDVD, DVDRVR, ADDS, TSTSEL, MULTBC, OVFLOW, 
*       NEWVEC, MGS_NRM

* 	It also calls some basic BLAS routines:
*	DCOPY, DSCAL, DDOT, DAXPY, IDAMAX, DGEMV, DINIT

*	For solving the small eigenproblem, the routine DSPEVX from 
*	LAPACK is used. DSPEVX is obtainable from NETLIB, together 
*	with a series of subroutines that it calls.
*
*	All the routines have IMPLICIT REAL*8          (A-H,O-Z)
*
*-----------------------------------------------------------------------
        IMPLICIT REAL*8          (A-H,O-Z)
	DIMENSION WORK(IWRSZ),IWORK(IIWSZ)
	DIMENSION ISELEC(LIM),IREV(*)
	LOGICAL HIEND

        save
*       (this is for the reverse communication)

*-----------------------------------------------------------------------
*  (Important to the following is the concept of NUME, the distance of
*   the index of the eigenpair wanted which is farthest from the 
*   extremes,i.e.,
*      if  lowest  eigepairs i1<i2<...<ik are wanted, NUME=ik
*      if highest eigenpairs i1<i2<...<ik are wanted, NUME=N-i1+1
*   where i1,...,ik are the indices of the wanted eigenpairs.
*   Obviously, NUME.GE.(No. of EiGenpairs wanted). )

*   on entry
*   -------
*   IRC         REVERSE COMMUNICATION integer.
*               IRC=0   --->>  Exit
*               IRC=1   --->>  Preconditioning
*               IRC=2   --->>  Matrix vector for initial setups.
*               IRC=3   --->>  Matrix vector in DVDRVR
*	   NOTE Usually in reverse communication the routine called
*	   ----	gives out two arrays w1 and w2. The first is used
*		as input for calculation Aw1->w2, or solve Mw2=w1.
*		Then these are passed back to the routine.
*		TO avoid extra work space, AND to be able to 
*		perform block operations (i.e., where w1 is m blocks)
*		we use the same WORK array, where BASIS and AB are 
*		stored in. An integer array with 6 integers is passed out, 
*		except IRC: IREV. 
*		IREV(1)=NB Simply NB operations must be performed
*		           from either the preconditioner or the matvec.
*			   Also the number of targeted eigenpairs.
*		IREV(2)=iw1
*		IREV(3)=iw2
*		        A*WORK(iw1:iw1+NB*N) --> WORK(iw2:iw1+NB*N)
*		        Solve NB rhsides: 
*		        M*WORK(iw2:iw1+NB*N) = WORK(iw1:iw1+NB*N)
*		        Obviously: iw1=ibasis+("indexB"-1)*N
*			           iw2=iab   +("indexD"-1)*N
*		for the matvec, while for the preconditioner it is 
*		in reverse order.
*		IREV(4)=ieigval  (To denote where the eigenvalues are)
*		IREV(5)=iscra3   (where the eigenvalue index starts)
*		        Simply the targeted eigenvalues are:
*		    	do i=1,NB
*			   indx = IWORK(iscra3+i-1)
*			   L_i  = WORK(ieigval+indx-1)
*                          solve (M-L_i) w2 = w1
*			enddo
*               The following is additional information that advanced
*		users may need.
*		IREV(6)=iscra1   (where the residuals of the corresponding
*			      eigenpairs are. According to the above index)
*		IREV(7)=iscra1+LIM  this will hold the \delta eps shifts
*		from robust preconditioning
*		
*   N           Order of the matrix.
*   LIM 	The upper limit on the dimension of the expanding basis.
*		NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed
*		only for LIM=NUME=N. The choice of LIM depends on the 
*		available workspace (see below). If the space is 
*		available it is preferable to have a large LIM, but not
*		larger than NUME$+$40.
*   NOC 	Number of Orthogonalization Constraints (vectors).
*	        These NOC vectors are given on input in the first
*		N*NOC elements of WORK.
*   ILOW	The index of the lowest eigepair to be computed. If 
*		(ILOW.LE.0).or.(ILOW.GT.N), the selected eigenpairs 
*		to be computed should be contained in array ISELEC.
*		(Modified on exit).
*   IHIGH 	The index of the highest eigenpair to be computed. 
*		Considered ONLY when ILOW is in the range 
*		(0.LT.ILOW.LE.N). (Modified on exit).
*   ISELEC	Array of size LIM holding the user specified indices
*		for the eigenpairs to be computed. Considered only when
*		(ILOW.LE.0).or.(ILOW.GT.N). The indices are read from 
*		the first position until a non positive integer is met.
*		   Example: if N=500, ILOW=0, and ISELEC(1)=495, 
*		   ISELEC(2)=497, ISELEC(3)=-1, the program will find 
*		   2 of the highest eigenpairs, pairs 495 and 497.
*		Any order of indices is acceptable (Modified on exit).
*   NIV		Number of Initial Vector estimates provided by the user.
*		If LIM > NIV > 0 then the NIV columns of size N of WORK 
*		coming after the NOC constraint	vectors,  should contain 
*               the estimates (see below). 
*               If NIV = 0, the vector (1,1,1,...,1)T is chosen.
*   MBLOCK	Number of vectors to be targeted in each iteration. 
*		1.LE.MBLOCK.LE.(No. EiGenpairs wanted) should hold.
*		Large block size reduces the number of iterations
*		(matrix acceses) but increases the matrix-vector
*		multiplies. It should be used when the matrix accese
*		is expensive (disc, recomputed or distributed).
*   CRITE	Convergence threshold for eigenvalues. 
*		If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted
*		eigenvalues, convergence is signaled.
*   CRITC       Convergence threshold for the coefficients of the last
*               added basis vector(s). If all of those corresponding
*               to unconverged eigenpairs are less than CRITC convergence
*               is signaled.
*   CRITR	Convergence threshold for residual vector norms. If 
*		all the residual norms ||Ax_i-l_ix_i|| of the targeted
*		x_i are less than CRITR convergence is signaled.
*		If ANY of the criteria are satisfied the algorithm stops
*   MAXITER	Upper bound on the number of iterations of the 
*		algorithm. When MAXITER is exceeded the algorithm stops.
*		A typical MAXITER can be MAX(200,NUME*40), but it can 
*		be increased as needed.
*   WORK	Real array of size IWRSZ. Used for both input and output
*		If NIV is in (1.LE.(NIV).LE.(LIM)), on input, WORK
*		must have the NIV initial estimates. These NIV N-element
*		vectors start from WORK(1) and continue	one after the 
*		other. They must form an orthonormal basis.
*   IWRSZ	The size of the real workspace. It must be at least as 
*		large as:
*
*		Not this any more:
*		     (2*LIM+NOC)*N + LIM*LIM + (NUME+10)*LIM + NUME 
*		New workspace:
*     --> this       (2*LIM+NOC)*N + LIM*LIM*2 +      11*LIM + NUME 
*
*   IWORK 	Integer work array of size IIWSZ. Used as scrath array
*		for indices and for use in the LAPACK routines.
*   IIWSZ	The size of the integer workspace. It must be at least 
*		as large as:
*			             6*LIM + NUME
*		
*		If LIM or NUME needs to be increased, the space should
*               also be increased accordingly. For given IWRSZ and 
*		IIWSZ one can calculate how big a problem one can 
*		solve (LIM,NUME).
*
*   on exit
*   -------
*   WORK(NOC*N+1) 
*		The first NUME*N locations contain the approximations to
*		the NUME extreme eigenvectors. If the lowest eigenpairs
*		are required, (HIEND=false), eigenvectors appear in 
*		ascending order, otherwise (HIEND=false), they appear in
*		descending order. If only some are requested, the order
*		is the above one for all the NUME extreme eigenvectors,
*		but convergence has been reached only for the selected 
*		ones. The rest are the current approximations to the 
*		non-selected eigenvectors.
*   WORK(NOC*N+NUME*N+1)
*		The next NUME locations contain the approximations to 
*		the NUME extreme eigenvalues, corresponding to the above
*		NUME eigenvectors. The same ordering and convergence 
*		status applies here as well.
*   WORK(NOC*N+NUME*N+NUME+1)
*		The next NUME locations contain the corresponding values
*		of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of
*		the last step of the algorithm.
*   WORK(NOC*N+NUME*N+NUME+NUME+1)
*		The next NUME locations contain the corresponding 
*		residual norms of the NUME above eigenvectors, of the 
*		last step.
*   NLOOPS      Number of iterations (accesses to matrix)
*
*   HIEND	Logical. If .true. on exit the highest eigenpairs are
*		found in descending order. Otherwise, the lowest 
*		eigenpairs are arranged in ascending order.
*   IERR	An integer denoting the completions status:
*		IERR = 0 	denotes normal completion.
*		IERR = -k 	denotes error in DSPEVX (k eigenpairs 
*			 	not converged)
*		0<IERR<=2048	denotes some inconsistency as follows:
*	 If (INT( MOD(IERR,  2)/1  ) N < LIM
*	 If (INT( MOD(IERR,  4)/2  ) LIM < 1
*	 If (INT( MOD(IERR,  8)/4  ) ISELEC(1)<1, and no range specified
*	 If (INT( MOD(IERR, 16)/8  ) IHIGH > N (in range or ISELEC)
*	 If (INT( MOD(IERR, 32)/16 ) IHIGH < ILOW (Invalid range)
*	 If (INT( MOD(IERR, 64)/32 ) NEIG >= LIM (Too many wanted)
*	 If (INT( MOD(IERR,128)/64 ) Probable duplication in ISELEC
*	 If (INT( MOD(IERR,256)/128) NUME >= LIM (max eigen very far)
*	 If (INT( MOD(IERR,512)/256) MBLOCK is out of bounds
*	 If (INT( MOD(IERR,1024)/512) IWRSZ or IIWSZ is not enough
*	 If (INT( MOD(IERR,2048)/1024) Orthogonalization Failed
*	 If (INT( MOD(IERR,4096)/2048) NLOOPS > MAXITER
*	 If (INT( MOD(IERR,8192)/4096) Not proper number of NIV
*	
*		The program will also print an informative message to 
*		the standard output when NIV is not proper 
*-----------------------------------------------------------------------
c	SUBROUTINE DVDSON(IRC,IREV,N,LIM,NOC,
c     :			 ILOW,IHIGH,ISELEC,NIV,MBLOCK,
c     :			 CRITE,CRITC,CRITR,MAXITER, 
c     :			 WORK,IWRSZ,IWORK,IIWSZ,
c     :			 HIEND,NLOOPS,IERR)      
*-----------------------------------------------------------------------
*       Reverse Communication
*-----------------------------------------------------------------------
        If ((IRC.eq.1).or.(IRC.eq.3)) then
*          ..(=1)preconditioning. Go into DVDRVR
*          ..(=3)or Matrix vector in DVDRVR
           goto 200
        else if (IRC.eq.2) then
*          ..matrix vector for setting up.
           goto 100
        endif
*-----------------------------------------------------------------------
*
* Checking user input errors, and setting up the problem to solve.
*
	IERR=0
	IF (LIM.GT.N) IERR=IERR+1
	IF (LIM.LE.0) IERR=IERR+2

	HIEND=.false.

	IF ((ILOW.LE.0).OR.(ILOW.GT.N)) THEN
*          ..Look for user choice of eigenpairs in ISELEC
	   IF (ISELEC(1).LE.0) THEN
*             ..Nothing is given in ISELEC
	      IERR=IERR+4
	   ELSE 
*             ..Find number of eigenpairs wanted, and their 
*	      ..min/max indices
	      NEIG=1
	      ILOW=ISELEC(1)
	      IHIGH=ISELEC(1)
	      DO 10 I=2,LIM
		 IF (ISELEC(I).LE.0) GOTO 20
		 ILOW=MIN(ILOW,ISELEC(I))
		 IHIGH=MAX(IHIGH,ISELEC(I))
		 NEIG=NEIG+1
 10	      CONTINUE
* 	      ..Check if a very large index is asked for
 20	      IF (IHIGH.GT.N) IERR=IERR+8
	   ENDIF
     	ELSE
*          ..Look for a range between ILOW and IHIGH
*          ..Invalid range. IHIGH>N
	   IF (IHIGH.GT.N) IERR=IERR+8
	   NEIG=IHIGH-ILOW+1
*          ..Invalid range. IHIGH<ILOW
	   IF (NEIG.LE.0) IERR=IERR+16
	   IF (NEIG.GT.LIM) THEN
* 	      ..Not enough Basis space. Increase LIM or decrease NEIG
              IERR=IERR+32
	   ELSE 
*  	      ..Fill in the ISELEC with the required indices
	      DO 40 I=1,NEIG
 40		 ISELEC(I)=ILOW+I-1  
	   ENDIF
	ENDIF

	IF (IERR.NE.0) RETURN

	NUME=IHIGH
*       ..Identify if few of the highest eigenpairs are wanted.
        IF ((ILOW+IHIGH-1).GT.N) THEN
           HIEND=.true.
           NUME=N-ILOW+1
*          ..Change the problem to a minimum eipenpairs one
*          ..by picking the corresponding eigenpairs on the
*          ..opposite side of the spectrum.
           DO 50 I=1,NEIG
 50           ISELEC(I)=N-ISELEC(I)+1
        ENDIF
*       ..duplications in ISELEC
        IF (NEIG.GT.NUME) IERR=IERR+64
*       ..Not enough Basis space. Increase LIM or decrease NUME
        IF ((NUME.GT.LIM).OR.((NUME.EQ.LIM).AND.(NUME.NE.N)))
     :     IERR=IERR+128
*	..Size of Block out of bounds
	IF ( (MBLOCK.LT.1).OR.(MBLOCK.GT.NEIG) ) IERR=IERR+256

*       ..Check for enough workspace for Dvdson
	IF ((IWRSZ.LT.(LIM*(2*N+LIM+(NUME+10))+NUME+NOC*N)).OR.
     :      (IIWSZ.LT.(6*LIM+NUME))) IERR=IERR+512

	IF (NIV.GT.LIM) THEN
* 	   ..Check number of initial estimates NIV is lower than LIM.
	   PRINT*,'WARNING: Too many initial estimates.?'
           PRINT*,'The routine needs at most:',LIM
	   IERR=IERR+4096
	ELSEIF (NIV.LT.NUME) THEN
* 	   ..check if enough initial estimates. 
     	   PRINT*,'WARNING: Not enough initial estimates'
           PRINT*,NUME-NIV,' Lanczos vectors will be added'
	ENDIF

	IF (IERR.NE.0) RETURN
*  
* Assigning space for the real work arrays
*
	iOrtho    =1
	iBasis    =iOrtho  + N*NOC
	ieigval   =iBasis  + N*LIM
	iAB	  =ieigval + LIM
	iS	  =iAB     + N*LIM
  	itempS	  =iS      + LIM*(LIM+1)/2
	iSvec	  =itempS  + LIM*(LIM+1)/2
CCC	iscra1    =iSvec   + LIM*(NUME+1)
	iscra1    =iSvec   + LIM*LIM
	ioldval   =iscra1  + 8*LIM
*
* Assigning space for the integer work arrays
*
	iscra2    =1
	iscra3    =iscra2  +5*LIM
	iIcv      =iscra3  +LIM
*
* Initialize tha basis, the AB, the S. 
*
	
 100	call initdvd(IRC,IREV,N,NOC,NIV,NUME+1,LIM,HIEND,WORK(iscra1),
     :		WORK(iOrtho),WORK(iBasis),WORK(iAB),WORK(iS))
*       ----------------------------------------------------------------
*       ..Reverse Communication for possible matrix vector
        if (IRC.eq.2) then
           IREV(2) = iBasis + (IREV(2)-1)*N
           IREV(3) = iAB    + (IREV(3)-1)*N
           return
        endif
*       ----------------------------------------------------------------
*
* Call main driver routine.
*
	NLOOPS=1

 200	CALL DVDRVR(IRC,IREV,N,HIEND,LIM,MBLOCK,NOC,
     :		   NUME,NIV,NEIG,ISELEC,
     :	    	   CRITE,CRITC,CRITR,MAXITER,
     :		   WORK(ieigval),WORK(iBasis),WORK(iOrtho),WORK(iAB),
     :		   WORK(iS),WORK(itempS),WORK(iSvec),
     :		   WORK(iscra1),IWORK(iscra2),IWORK(iscra3),
     :		   IWORK(iIcv),WORK(ioldval),
     :		   NLOOPS,IERR)
*       ----------------------------------------------------------------
*       some Reverse Communication
        if (IRC.eq.1) then
*	   ..Preconditioning 
	   IREV(2) = iAB    + (IREV(2)-1)*N
           IREV(3) = iBasis + (IREV(3)-1)*N
	   IREV(4) = ieigval
	   IREV(5) = iscra3
	   IREV(6) = iscra1
	   IREV(7) = iscra1 + LIM
           return
	else if (IRC.eq.3) then
*	   ..Matrix-vector
	   IREV(2) = iBasis + (IREV(2)-1)*N
           IREV(3) = iAB    + (IREV(3)-1)*N
	   return
	endif
*       ----------------------------------------------------------------

	IF (HIEND) CALL DSCAL(NUME,-1.D0,WORK(ieigval),1)
*
* -Copy the eigenvalues after the eigenvectors
* -Next, copy the difference of eigenvalues between the last two steps
* -Next, copy the residuals for the first NUME estimates
*
	CALL DCOPY(NUME,WORK(ieigval),1,WORK(iBasis+N*NUME),1)
	CALL DCOPY(NUME,WORK(ioldval),1,WORK(iBasis+(N+1)*NUME),1)
	CALL DCOPY(NUME,WORK(iscra1),1,WORK(iBasis+(N+2)*NUME),1)
*
* Set IRC=0 for normal exit with no reverse communication
*
        IRC=0
 	RETURN
	END
*=======================================================================
	SUBROUTINE ADDS(N,LIM,KPASS,NNCV,BASIS,AB,S)
*=======================================================================
*	Called by: DVDSON
*
*	Calculates the new column in the S matrix. S has a 
*	new row and column, but being symmetric only the new column is 
*	stored. S(i,kpass+1)=B(i)^T D(kpass+1) for all i.
*
*	subroutines called:
*	DDOT, DSCAL
*	
*-----------------------------------------------------------------------
	IMPLICIT REAL*8          (A-H,O-Z)
	DIMENSION BASIS(N*LIM),AB(N*LIM)
	DIMENSION S(LIM*(LIM+1)/2)
*-----------------------------------------------------------------------
*   on entry
*   -------
*   N		The order of the matrix A
*   kpass	The current dimension of the expanding sub-basis
*   NNCV	Number of new basis vectors.
*   Basis	the basis vectors, including the new NNCV ones.
*   AB		The matrix D=AB. (with new NNCV columns)
*   on exit
*   -------
*   S           The small matrix with NNCV new columns at the last part
*-----------------------------------------------------------------------
*
* The new S is calculated by adding the new last columns
* S(new)=B^T D(new).
*
        IDSTART=KPASS*N+1
	ISSTART=KPASS*(KPASS+1)/2
	DO 20 IV=1,NNCV
	   IBSTART=1
	   DO 10 IBV=1,KPASS+IV
	       SS=DDOT(N,BASIS(IBSTART),1,AB(IDSTART),1)
	       S(ISSTART + IBV)=SS
	       IBSTART=IBSTART+N
 10        CONTINUE
	   ISSTART=ISSTART+KPASS+IV
	   IDSTART=IDSTART+N
 20	CONTINUE

	RETURN
	END
*=======================================================================
	SUBROUTINE DVDRVR(IRC,IREV,N,HIEND,LIM,MBLOCK,NOC,
     :			  NUME,NIV,NEIG,ISELEC,
     :			  CRITE,CRITC,CRITR,MAXITER,
     :                    EIGVAL,BASIS,OrthoBasis,AB,S,TEMPS,SVEC,
     :			  SCRA1,ISCRA2,INCV,ICV,OLDVAL,
     :			  NLOOPS,IERR)
*=======================================================================
*	called by DVDSON
*
*	Driver routine implementing Davidson's main loop. On entry it 
*	is given the Basis, the work matrix D=AB and the small symmetric
*	matrix to be solved, S=B^TAB (as found by ADDS). In each step 
*	the small problem is solved by calling DSPEVX. 
*	TSTSEL tests for eigenvalue convergence and selects the next 
*	pairs to be considered for targeting (as a block). 
*      	NEWVEC computes the new vectors (block) to be added in the 
*	expanding basis, and tests for residual convergence.
*	After Preconditioning ((M-l_i)d=res), the basis is orthogonalized.
*       Then the matrix-vector multiplication finds the new vectors of D 
*       (Dnew=ABnew), and the new small problem S, is calculated. 
*       The algorithm is repeated. 
*	In case of a large expanding basis (KPASS=LIM) the Basis, AB, 
*	SVEC and S are collapsed.
*	At the end the current eigenvector estimates are computed as 
*	well as the residuals and eigenvalue differences.
*
*	Subroutines called:
*	DSPEVX, MULTBC, TSTSEL,	OVFLOW, NEWVEC, ADDS,
*	DCOPY, DDOT, DAXPY
*-----------------------------------------------------------------------
        IMPLICIT REAL*8          (A-H,O-Z)
	DIMENSION IREV(*)
	DIMENSION S(LIM*(LIM+1)/2),TEMPS(LIM*(LIM+1)/2)
CCC	DIMENSION SVEC(LIM*(NUME+1)),EIGVAL(LIM)
	DIMENSION SVEC(LIM*LIM),EIGVAL(LIM)
	DIMENSION ISELEC(NEIG)
	DIMENSION BASIS(N*LIM),OrthoBasis(N*LIM+NOC*N),AB(N*LIM)
	DIMENSION SCRA1(8*LIM),ISCRA2(5*LIM),INCV(LIM)
	DIMENSION ICV(NUME+1),OLDVAL(NUME+1)
	LOGICAL FIRST,DONE,HIEND,TSTSEL
	INTEGER rest

        save
        !include 'mpif.h'
        common/mpi/myid, nprocs, irrr
*       (this is for the reverse communication)

*-----------------------------------------------------------------------
*
*   on entry
*   -------
*
*   N           The order of the matrix A
*   HIEND	Logical. True only if the highest eigenpairs are needed.
*   LIM 	The limit on the size of the expanding Basis
*   MBLOCK	Number of vectors to be targeted in each iteration.
*   NOC 	Number of orthogonalization constraints (vectors)
*   NUME        The largest index of the eigenvalues wanted.
*   NIV		Starting dimension of expanding basis.
*   NEIG 	Number of eigenvalues wanted.
*   ISELEC	Array containg the indices of those NEIG eigenpairs.
*   CRITE       Convergence thresholds for eigenvalues, coefficients
*   CRITC,CRITR	and residuals.
*   BASIS 	Array with the basis vectors.
*   AB       	Array with the vectors D=AB
*   S		Array keeping the symmetric matrix of the small problem.
*   TEMPS	scratch array
*   SVEC	Array for holding the eigenvectors of S
*   SCRA1 	Srcatch array used by DSPEVX. It also holds residuals
*		and eigenvalue shifts for preconditioning.
*   ISCRA2	Integer Srcatch array used by DSPEVX.
*   INCV	Srcatch array used in DSPEVX. Also used in TSTSEL and 
*		NEWVEC where it holds the Indices of uNConVerged pairs
*   ICV		It contains "1" to the locations of ConVerged eigenpairs
*   OLDVAL	Array keeping the previous' step eigenvalue estimates.
*
*   on exit
*   -------
*
*   EIGVAL	Array containing the NUME lowest eigenvalues of the 
*		the matrix A (or -A if the highest are sought).
*   Basis 	On exit Basis stores the NUME corresponding eigenvectors
*   OLDVAL 	On exit it stores the final differences of eigenvalues.
*   SCRA1	On exit it stores the NUME corresponding residuals.
*   NLOOPS	Number of loops taken by the algorithm
*
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*	Reverse Communication
        if (IRC.eq.1) then
*          ..Came from preconditioning
           goto 100
        else if (IRC.eq.3) then
*          ..came from matrix vector multiply
           goto 200
        endif
*-----------------------------------------------------------------------

 	DO 5 I = 1,NUME
	   EIGVAL(I) = 1.D30
  5        ICV(I) = 0
	FIRST = .true.
	KPASS = NIV
	NNCV  = KPASS
C
C Decide HERE how many to restart with, not more than LIM-3.
C
 	rest = max(lim/2,2*nume)
	rest = min(rest, lim-3)

10 	CONTINUE 
*	(iterations for kpass=min(NIV,NUME),LIM)
*
* Diagonalize the matrix S. Find only the IFIND smallest eigenpairs
*
C
C Decide on how many to solve for:
C
	   IFIND = MIN(kpass,rest)
              tol = slamch('S')
	   CALL DCOPY(IFIND,EIGVAL,1,OLDVAL,1)
	   CALL DCOPY((KPASS*(KPASS+1))/2,S,1,TEMPS,1)
	   CALL DSPEVX('Vectors also','In a range','Upper triangular',
     :		KPASS,TEMPS,-1.,-1.,1,IFIND,tol,
     :		NFOUND,EIGVAL,SVEC,KPASS,SCRA1,ISCRA2,INCV,INFO)
	   IERR = -ABS(INFO)
	   IF (IERR.NE.0) GOTO 60
*
* TeST for convergence on the absolute difference of eigenvalues between
* successive steps. Also SELect the unconverged eigenpairs and sort them
* by the largest magnitude in the last added NNCV rows of Svec.
*
	   DONE = TSTSEL(KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV,
     :		       CRITE,CRITC,SCRA1,ISCRA2,OLDVAL,NNCV,INCV)
           IF ((DONE).OR.(KPASS.GE.N)) GOTO 30
*
* Maximum size for expanding basis. Truncate basis, D, and S, Svec
* Consider the basis vectors found in TSTSEL for the newvec. KPASS=NUME
*

*          Change suggested by Anreas, March 18, 1996
           IF (KPASS.EQ.LIM) THEN
!      PRINT*,'collapsing the basis: lim,nume,rest',lim,nume,rest
!        PRINT*,'myid = ', myid, '  nprocs = ', nprocs
!        PRINT*,'collapsing the basis: lim,nume,rest',lim,nume,rest,myid
!23456789012345678901234567890123456789012345678901234567890123456789012
              CALL OVFLOW(N,rest,KPASS,SCRA1,BASIS,AB,S,SVEC,EIGVAL)
*             CALL OVFLOW(N,NUME,KPASS,SCRA1,BASIS,AB,S,SVEC,EIGVAL)
	   END IF
*
* Compute and add the new residuals. NNCV is set to the number of new
* vectors that have not converged. If none, DONE=true, exit.
*
 	   CALL NEWVEC(N,NUME,LIM,MBLOCK,KPASS,CRITR,NNCV,INCV,
     :	      SVEC,EIGVAL,OLDVAL,AB,BASIS,ICV,SCRA1,SCRA1(LIM+1),DONE)

	   IF (DONE) GOTO 30
*-----------------------------------------------------------------------
* Preconditioning the NNCV (residuals-deps x) stored in AB(kpass+1).
*          ..Robust Preconditioning (Eigenvalue shift). 
*          ..Look for lowest, so Li-|eps_i|
           do i=1,nncv
	      EIGVAL(incv(i)) = EIGVAL(incv(i)) - scra1(LIM+incv(i))
	   enddo
*
*          ..Change sign of eigenvalues if HIEND.
           IF (HIEND) CALL DSCAL(NUME,-1.D0,EIGVAL,1)

* Use of Reverse Communication.
*
	   IREV(1) = NNCV
           IREV(2) = kpass + 1
           IREV(3) = kpass + 1
	   IRC = 1
	   if (IRC.eq.1) return
*          ..Continue from preconditioning
 100       continue
*
*          ..Change back sign of eigenvalues if HIEND.
           IF (HIEND) CALL DSCAL(NUME,-1.D0,EIGVAL,1)
*
*          ..Shift the eigenvalues back to what they were.
           do i=1,nncv
	      EIGVAL(incv(i)) = EIGVAL(incv(i)) + scra1(LIM+incv(i))
	   enddo
*-----------------------------------------------------------------------
*
* Orthonormalization of the previous vectors to the Basis and to any 
* orthogonalization constraints. The not-yet-filled
* spaces of AB (from NEWSATRT onwards) are used for scratch.
*
	   NEWSTART = KPASS*N+1

	   call mgs_nrm(N,noc+kpass,nncv,scra1(LIM+1),OrthoBasis)
*-----------------------------------------------------------------------
* Use of Reverse Communication. Add new columns in D through matrix vector 
* multiplication CALL OP(N,NNCV,BASIS(NEWSTART),AB(NEWSTART))
*
           IREV(1) = NNCV
           IREV(2) = kpass + 1
           IREV(3) = kpass + 1
           IRC = 3
           return

*          ..Continue from matrix-vector multiply
 200 	   continue
*-----------------------------------------------------------------------
*
* 	   ..If highest pairs are sought, use negative of the matrix
*
	   IF (HIEND) CALL DSCAL(N*NNCV,-1.D0,AB(NEWSTART),1)
*
* Add new column in S, from the NNCV new vectors.
*
	   CALL ADDS(N,LIM,KPASS,NNCV,BASIS,AB,S)

           KPASS  = KPASS+NNCV
	   NLOOPS = NLOOPS+1

	IF (NLOOPS.LE.MAXITER) GOTO 10
	IERR   = IERR+2048
	NLOOPS = NLOOPS-1
	KPASS  = KPASS-NNCV
 30	CONTINUE
*
* Calculate final results. EIGVAL contains the eigenvalues, BASIS the
* eigenvectors, OLDVAL the eigenvalue differences, and SCRA1 residuals.
*
	DO 40 I = 1,NUME
 40	   OLDVAL(I) = ABS(OLDVAL(I)-EIGVAL(I))
	
	CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,BASIS)
	CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,AB)
*
* i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi
* temporarily stored in AB(NUME*N+1)
*
	DO 50 I = 1,NUME
	   CALL DCOPY(N,AB((I-1)*N+1),1,AB(NUME*N+1),1)
	   CALL DAXPY(N,-EIGVAL(I),BASIS((I-1)*N+1),1,AB(NUME*N+1),1)
	   SCRA1(I) = DDOT(N,AB(NUME*N+1),1,AB(NUME*N+1),1)
	   SCRA1(I) = SQRT(SCRA1(I))
 50	CONTINUE
*
* Set IRC=0 for normal exit with no reverse communication
*
       if (myid.eq.0) print*, 'DVDSON::NLOOPS =', NLOOPS
 60     IRC = 0
        RETURN
	END
*=======================================================================
        subroutine initdvd(IRC,IREV,N,NOC,NIV,NUME,LIM,HIEND,SCRA1,
     :			OrthoBasis,Basis,AB,S)
*=======================================================================
*	Initializes the basis and the auxiliary arrays AB and S.
*	If not enough initial estimates exist the basis will be 
*	supplemented with Lanczos vectors of the current NIV vectors.
*	eg. (b1,b2) --> (b1,Ab1,AAb1, b2,Ab2,AAb2) for a NUME of 6.

*  	OrthoBasis is the Basis with the NOC orthogonalization 
*	constraint vectors in the begining. This equialence holds:
*	   	  equivalence(Basis(1),OrthoBasis(NOC*N+1))
*-----------------------------------------------------------------------
	implicit REAL*8          (a-h,o-z)
	dimension IREV(*)
	dimension OrthoBasis(N*(NOC+LIM)),Basis(N*LIM),AB(N*LIM)
	dimension S(*),SCRA1(*)
	logical HIEND

	save
*-----------------------------------------------------------------------
*       ..Reverse Communication
	if (IRC.eq.2) goto 100
*-----------------------------------------------------------------------
*
* If no initial estimates pick one 
*
        if (NIV.eq.0) then
           call dinit(N,1.D0/SQRT(DBLE(N)),Basis,1)
	   niv = 1
        endif
*
* Compute AB. Also fill basis with orthonormalized ABs until enough NIVs.
*
        ist  = 1
        iend = niv
	new  = niv
 10 	continue
*       ----------------------------------------------------------------
*       ..Reverse Communication. Matrix Vector of A and Basis.
*
        IRC     = 2
        IREV(1) = new
        IREV(2) = ist
        IREV(3) = ist
        return
 100    continue
*       ----------------------------------------------------------------

	if (iend.lt.NUME) then
	   new = min( nume-iend, iend-ist+1)
           call dcopy(N*new,AB((ist-1)*N+1),1,Basis(1+iend*N),1)
*            ..orthonormalize OrthoBasis i.e.,
*            B-1..B-noc,B1,...,Biend,( Biend+1,...Biend+new )
	   call mgs_nrm(N,noc+iend,new,scra1,OrthoBasis)
	   ist  = iend + 1
           iend = iend+new
	   goto 10
	endif
	niv = iend
	!xhh print*, 'niv=',niv, 'nume=', nume
*
* Scale if HIEND for highest eigepairs
*
        IF (HIEND) CALL DSCAL(N*NIV,-1.D0,AB,1)
*
* Also find the small matrix S = B^TAB.
*
        KPASS = 0
        CALL ADDS(N,LIM,KPASS,NIV,Basis,AB,S)

	IRC = 0
        return
	end

*=======================================================================
        SUBROUTINE MULTBC(N,K,M,C,TEMP,B)
*=======================================================================
*	called by: DVDRVR
*
*	Multiplies B(N,K)*C(K,M) and stores it in B(N,M)
*	Used for collapsing the expanding basis to current estimates,
*	when basis becomes too large, or for returning the results back

*       Subroutines called
*       DINIT, DGEMV, DCOPY
*-----------------------------------------------------------------------
        IMPLICIT REAL*8          (A-H,O-Z)
	DIMENSION B(N*K),C(K*M),TEMP(M)
*-----------------------------------------------------------------------
        DO 10 IROW=1,N
	   CALL DGEMV('Transp',K,M, 1.D0, C,K,B(IROW),N, 0.D0 ,TEMP,1)
           CALL DCOPY(M,TEMP,1,B(IROW),N)
  10    CONTINUE

        RETURN
        END
*=======================================================================
	SUBROUTINE NEWVEC(N,NUME,LIM,MBLOCK,KPASS,CRITR,NNCV,INCV,
     :		  SVEC,EIGVAL,OLDVAL,AB,BASIS,ICV,SCRA1,EPSIL,DONE)
*=======================================================================
*
*	Called by: DVDRVR
*
*	It calculates the new expansion vectors of the basis. 
*	For each one of the vectors in INCV starting with the largest
*	megnitude one, calculate its residual Ri= DCi-liBCi and check 
*	the ||Ri|| for convergence. If it is converged do not add it 
*	but look for the immediate largest coefficient and its vector.
*	The above procedure continues until MBLOCK vectors have been
*	added to the basis, or the upper limit has been encountered.
*	Thus only  the required MBLOCK residuals are computed. 
*	For robust Preconditioning, each residual is then modified as 
*		         \delta_eps x - res
*	where \delta_eps is an estimate of the eigenvalue error.
*	A different estimate is used for the eigenvalue shift.
*
*	Subroutines called:
*	DNRM2, DGEMV
*
*-----------------------------------------------------------------------
	IMPLICIT REAL*8          (A-H,O-Z)
	DIMENSION INCV(NUME)
	DIMENSION ICV(NUME)
	DIMENSION BASIS(N*LIM),AB(N*LIM)
	DIMENSION SVEC(LIM*NUME)
	DIMENSION EIGVAL(LIM),OLDVAL(NUME),SCRA1(LIM),EPSIL(LIM)
	LOGICAL DONE
*-----------------------------------------------------------------------
*   on entry
*   -------- 
*   N           The order of the matrix A
*   NUME        The largest index of the eigenvalues wanted.
*   LIM         The limit on the size of the expanding Basis
*   MBLOCK	Maximum number of vectora to enter the basis
*   KPASS	the current dimension of the expanding basis
*   CRITR	Convergence threshold for residuals
*   NNCV	Number of Non ConVerged pairs (MBLOCK will be targeted)
*   INCV        Index to the corresponding SVEC columns of these pairs.
*   SVEC,EIGVAL	Arrays holding the eigenvectors and eigenvalues of S
*   AB		Array with the vectors D=AB
*   BASIS	the expanding basis having kpass vectors
*   ICV		Index of converged eigenpairs (ICV(i)=1 <=>i converged)

*   on exit
*   -------
*   NNCV	The number of vectors finally added to the basis.
*   BASIS 	The eigenector estimated reside at the end.
*   AB          The right hand sides for preconditioning.
*   EPSIL       The eigenvalue shifts for preconditioning,
*   ICV		Index of converged eigenpairs (updated)
*   DONE 	logical, if covergance has been reached.
*		the Basis should be collapsed to current approximations.
*-----------------------------------------------------------------------
        !include 'mpif.h' 
        common/mpi/myid, nprocs, irrr
	DONE 	= .FALSE.
	NEWSTART= KPASS*N+1
	NADDED  = 0
	ICVC    = 0
	LIMADD  = MIN( LIM, MBLOCK+KPASS )
	ICUR    = NEWSTART
*
* Compute RESIDUALS for the MBLOCK of the NNCV not converged vectors.
*
	DO 10 I=1,NNCV
	   INDX = INCV(I)
*
*          ..Compute b = BASIS*Svec_indx
*	   ..Compute d = AB*Svec_indx
*   	   ..Daxpy d'= d - eigval b  gives the residual
*
	   CALL DGEMV('N',N,KPASS,1.D0,BASIS,N,SVEC((INDX-1)*KPASS+1),1,
     :		       0.d0,BASIS(ICUR),1)
	   CALL DGEMV('N',N,KPASS,1.D0,AB,N,SVEC((INDX-1)*KPASS+1),1,
     :                 0.d0,AB(ICUR),1)
	   CALL DAXPY(N,-EIGVAL(INDX),BASIS(ICUR),1,AB(ICUR),1)
*
*	   ..Compute the norm of the residual
*	   ..and check for convergence
*
	   SQRES = DDOT(N,AB(ICUR),1,AB(ICUR),1)
	   SCRA1(INDX) = SQRT( SQRES )

          if (myid.eq.0) 
     :        print '(A11,f22.16,i2,A10,f18.16)', ' EIGVAL(i) ', 
     :        EIGVAL(INDX),indx,' Res.Norm ',SCRA1(INDX)

	   IF (SCRA1(INDX).LT.CRITR) THEN
*	      ..Converged,do not add. Go for next non converged one
	      ICVC = ICVC+1
	      ICV( INDX ) = 1
	      IF (ICVC.LT.NNCV) GOTO 10
*	      ..All have converged.
		!print*, 'converged by critr'
		!print*, 'converged by critr',myid
	      DONE = .TRUE.
	      RETURN
	   ELSE
*	      ..Not converged. Consider it for preconditioning
*          ---------------  ROBUST MODIFICATION ---------------------
*	      ..Daxpy d'= d'- \delta_eps b  gives the rhs for precond.
*	      ..It is stored in AB.
*             ..Deps are also stored in EPSIL
*
              if (indx.eq.1) then 
                 gaplow = 1.0D+99
		 gap = abs( eigval(2) - eigval(1) )
*                 print*, 'h1,h2:',eigval(1),eigval(2)
	      else
		 gaplow = abs(eigval(indx)-eigval(indx-1))
		 gapup  = abs(eigval(indx+1)-eigval(indx))
                 gap = min( gaplow,gapup )
	      endif
	      dl = abs(OLDVAL(indx)-eigval(indx))
	      res = scra1(indx)

              if ( gap.gt.res) then
*                 EPSIL(indx) = min(dl,res,gaplow)
                 rhseps = sqrt(dl*res)
              else 
*                EPSIL(indx) = min( res, gaplow)
	         rhseps = min(dl,sqrt(dl))
	      endif

	      EPSIL(indx) = 0.d0
              CALL DAXPY(N,-rhseps,BASIS(ICUR),1,AB(ICUR),1)

*          -------------- END OF ROBUST MODIFICATIONS -----------------
*
	      NADDED = NADDED+1
	      INCV(NADDED) = INDX
	      IF ((NADDED+KPASS).EQ.LIMADD) GOTO 20
*	      ..More to be added in the block
	      ICUR = ICUR+N
	   ENDIF
 10	CONTINUE

 20	NNCV=NADDED

        RETURN
	END
*=======================================================================
	SUBROUTINE OVFLOW(N,NUME,KPASS,SCRA1,BASIS,AB,S,SVEC,EIGVAL)
*=======================================================================
*	Called by: DVDRVR
*	Called when the upper limit (LIM) has been reached for the basis
*	expansion (by KPASS). The basis is truncated B'=BC and 
*	similarly the array D'=DC. The new S is computed as 
*	S'(i,j)=l(i)delta(i,j) where l(i) eigenvalues, and delta of 
*	Kronecker, i,j=1,NUME. The new eigenvectors of the small matrix 
*	are the unit vectors.
*
*	Subroutines called:
*	DCOPY, DINIT, MULTBC
*-----------------------------------------------------------------------
	IMPLICIT REAL*8          (A-H,O-Z)
        DIMENSION SVEC(KPASS*NUME),S((KPASS*(KPASS+1))/2)
        DIMENSION EIGVAL(KPASS),BASIS(N*KPASS),AB(N*KPASS)
*-----------------------------------------------------------------------
*   on entry
*   -------
*   NUME	The largest index of eigenvalues wanted.
*   SVEC 	the kpass eigenvectors of the smaller system solved
*   EIGVAL 	the eigenvalues of this small system
*   on exit
*   -------
*   Basis
*   AB
*   S 		The new small matrix to be solved.
*-----------------------------------------------------------------------
* Truncate  the basis and the AB array.
*
	CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,BASIS)
        CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,AB)
*
* calculation of the new upper S=diag(l1,...,l_NUME) and 
* its matrix Svec of eigenvectors (e1,...,e_NUME)
*
	CALL DINIT((NUME*(NUME+1))/2,0.d0,S,1)
	CALL DINIT(NUME*NUME,0.d0,SVEC,1)
	IND=0
	ICUR=0
	DO 10 I=1,NUME
	   S(IND+I)=EIGVAL(I)
	   SVEC(ICUR+I)=1
	   ICUR=ICUR+NUME
   10	   IND=IND+I

	KPASS = NUME

	RETURN
	END
*=======================================================================
	LOGICAL FUNCTION TSTSEL(KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV,
     :			     CRITE,CRITC,ROWLAST,IND,OLDVAL,NNCV,INCV)
*=======================================================================
*
*	Called by: DVDRVR

*	It first checks if the wanted eigenvalues have reached 
* 	convergence and updates OLDVAL. Second, for each wanted and non
*	converged eigenvector, it finds the largest absolute coefficient
*       of the NNCV last added vectors (from SVEC) and if not coverged, 
*	places it in ROWLAST. IND has the corresponding indices. 
*	Third, it sorts ROWLAST in decreasing order and places the 
*	corresponding indices in the array INCV. The index array INCV
*       and the number of unconverged pairs NNCV, are passed to DVDRVR.
*		Later in NEWVEC only the first MBLOCK of NNCV pairs
*		will be targeted, since if ROWLAST(i)>ROWLAST(j) 
*		then approximately RESIDUAL(i)>RESIDUAL(j)
*
*	Subroutines called
*	IDAMAX
*-----------------------------------------------------------------------
	IMPLICIT REAL*8          (A-H,O-Z)
	LOGICAL DONE
	DIMENSION SVEC(KPASS*NUME),EIGVAL(NUME)
	DIMENSION ICV(NUME)
        DIMENSION ROWLAST(NEIG),IND(NEIG),OLDVAL(NUME)
	DIMENSION INCV(NEIG),ISELEC(NEIG)
*-----------------------------------------------------------------------
*
*   on entry
*   -------
*   KPASS       current dimension of the expanding Basis
*   NUME        Largest index of the wanted eigenvalues.
*   NEIG        number of wanted eigenvalues of original matrix
*   ISELEC	index array of the wanted eigenvalues.
*   SVEC        the eigenvectors of the small system
*   EIGVAL	The NUME lowest eigenvalues of the small problem
*   ICV		Index of converged eigenpairs.ICV(i)=1 iff eigenpair i
*		has converged, and ICV(i)=0 if eigenpair i has not.
*   CRITE,CRITC	Convergence thresholds for eigenvalues and coefficients
*   ROWLAST   	scratch array, keeping the largest absolute coefficient
*		of the NNCV last rows of Svec.
*   IND		scratch array, temporary keeping the indices of Rowlast
*   OLDVAL	The previous iteration's eigenvalues.
*   
*   on exit
*   -------
*   NNCV 	 Number of non converged eigenvectors (to be targeted)
*   INCV	 Index to these columns in decreasing order of magnitude
*   TSTSEL   	 true if convergence has been reached 
*
*-----------------------------------------------------------------------
        !include 'mpif.h'
        common/mpi/myid, nprocs, irrr

	DONE=.False.
*
* Test all wanted eigenvalues for convergence under CRITE
*
	NNCE=0
	DO 10 I=1,NEIG
	   IVAL=ISELEC(I)
 10    	   IF (ABS(OLDVAL(IVAL)-EIGVAL(IVAL)).GE.CRITE) NNCE=NNCE+1
	IF (NNCE.EQ.0) THEN
	   if(myid.eq.0)  print*, 'converged by crite'
	   TSTSEL=.TRUE.
	   RETURN
	ENDIF
*
* Find the maximum element of the last NNCV coefficients of unconverged
* eigenvectors. For those unconverged coefficients, put their indices
* to IND and find their number NNCV
*
	ICNT=0
	DO 30 I=1,NEIG
	   IF (ICV(ISELEC(I)).EQ.0) THEN
*             ..Find coefficient and test for convergence
              ICUR=KPASS*ISELEC(I)
              TMAX=ABS( SVEC(ICUR) )
              DO 20 L=1,NNCV-1
 20              TMAX=MAX( TMAX, ABS(SVEC(ICUR-L)) )
              IF (TMAX.LT.CRITC) THEN
*                ..this  coefficient converged
                 ICV(ISELEC(I))=1
              ELSE
*                ..Not converged. Add it to the list.
                 ICNT=ICNT+1
                 IND(ICNT)=ISELEC(I)
                 ROWLAST(ICNT)=TMAX
              ENDIF
           ENDIF
 30     CONTINUE

	NNCV=ICNT
	IF (NNCV.EQ.0) DONE=.TRUE.
	!IF (NNCV.EQ.0) print*, 'converged by critc', kpass
	IF (NNCV.EQ.0) print*, 'converged by critc', kpass,myid
*
* Sort the ROWLAST elements interchanging their indices as well
*
	DO 40 I=1,NNCV
           INDX=IDAMAX(NNCV-I+1,ROWLAST(I),1)
	   INCV(I)=IND(INDX+I-1)
	
	   TEMP=ROWLAST(INDX+I-1)
	   ROWLAST(INDX+I-1)=ROWLAST(I)
	   ROWLAST(I)=TEMP
	   ITEMP=IND(INDX+I-1)
	   IND(INDX+I-1)=IND(I)
	   IND(I)=ITEMP
 40	CONTINUE

	TSTSEL=DONE
	RETURN
	END
*=======================================================================
	subroutine mgs_nrm(N,kp,new,scra,B)
*=======================================================================
*       Orthogonalizis the new vectors in B (from kp+1...kp+new)
*	to the previous kp B vectors and to themselves
*	using Modified Gram Schmidt. Then normalizes them.
*       The procedure is repeated twice.
*-----------------------------------------------------------------------
	implicit REAL*8          (a-h,o-z)
	dimension B((kp+new)*N),scra(new)
      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr

*
* MODIFIED GRAM-SCHMIDT  (twice)
*
      if (myid == 0) then
	do i=1,2

	newstart = kp*N + 1
*
*  First record contribution from the kp vectors to the new ones.	
*
	kcur = 1
	do k=1,kp
	   jcur = newstart
	   call dgemv('T',N,new,1.d0,B(jcur),N,
     :                 B(kcur),1, 0.d0,scra,1)
	   do j=1,new
!	      call daxpy(N,-scra(j),B(kcur),1,B(jcur),1)
              do mm = 0, N-1
                  B(jcur+mm) = B(jcur+mm) - scra(j)*B(kcur+mm)
              end do
	      jcur = jcur + N
	   enddo
	   kcur = kcur + N
	enddo
*
*  Then orthogonalize the new ones among themselves.
*
	do k=1,new
           jcur = kcur + N
*  The current vector should be normalized
*
	   dnm = ddot(N,B(kcur),1,B(kcur),1)
	   dnm = sqrt(dnm)
	   call dscal(n,1/dnm,B(kcur),1)
*
           call dgemv('T',N,new-k,1.d0,B(jcur),N,
     :                 B(kcur),1, 0.d0,scra,1)
	   do j=k+1,new
!              call daxpy(N,-scra(j-k),B(kcur),1,B(jcur),1)
              do mm = 0, N-1
                  B(jcur+mm) = B(jcur+mm) - scra(j-k)*B(kcur+mm)
              end do
              jcur = jcur + N
	   enddo
           kcur = kcur + N
	enddo
	enddo
        end if 

        call MPI_BCAST(B,(kp+New)*n,MPI_DOUBLE_PRECISION,0,
     :                 MPI_COMM_WORRLD,ierr)
	return
        end


