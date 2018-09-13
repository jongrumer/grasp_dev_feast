************************************************************************
*                                                                      *
      SUBROUTINE GETSCD (EOL, idblk, isofile, rwffile)
      IMPLICIT REAL*8          (A-H,O-Z)
      LOGICAL   EOL
      CHARACTER idblk(*)*8, isofile*(*), rwffile*(*)
*                                                                      *
*   Interactively determines the data governing the SCF problem.       *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, GETRSL, GETYN, NUCPOT, RADGRD,        *
*                        SETISO, SETQIC, SETRWF.                       *
*               [RSCF92]: GETALD, GETOLD
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 27 Dec 1992   *
*
* Xinghong He 98-08-06
*
************************************************************************

*
* Grid parameters are tradeoff between cost and precision
* 
* If you want a 'system dependent default grid', you should consider
* both Z and (Z-N) (ionicity? ionization?) dependency.
* The radial extent of neutral atom changes by a factor of 10 from He to Cs,
* while He-like Cs ion is much smaller than neutral He atom.
*
* The grid values which we recommend here are for neutral atoms.
* Default RNT (first grid point) should depend on Z, but the dependence
* on ionization stage (Z-N) is equally important.
* Step size 'H' should also depend on (Z-N), and also on Z.
* He-like Cs is small, but it still has to squeeze all those nodes within
* its radial extent, which means small RNT, small H, large NNNP (many points).
*
* RNT, H, and NNNP (number of tabulation points) are not independent.
* A simplified procedure for developing an appropriate grid
* is to start with very large NNNP, small RNT, and small H.
* Then slowly increase RNT and H until one of the two things happen:
* desired accuracy becomes unattainable or the code becomes unstable.
* This procedure should then be repeated with correlation orbitals,
* which makes the excercise somewhat tricky, but it can be done with
* resonably limited CSF expansion.
* Correlation orbitals require finer grid than spectroscopic orbitals.
* On the other hand, the 'desired accuracy' is usually much lower for
* high-n correlation orbitals, as compared to spectroscopic orbitals.
* When both RNT and H are set (one way or another), the number of
* tabulation points NNNP is determined by: RNT, H, and the radial extent
* of the system in question. It can be seen on output.
* NNNP and NNNW should be trimmed before large-scale runs.
*
* Jacek Bieron 2011 07 13
*
************************************************************************

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      CHARACTER*(*), PARAMETER:: MYNAME = 'GETSCD' 

      COMMON/iounit/istdi,istdo,istde

      INTEGER*4 IQADUM
      POINTER (PNTRIQ,IQADUM)
      LOGICAL DIAG,GETYN,LFIX,LFORDR,NOINVT,ORTHST,YES
      CHARACTER*80 RECORD
      CHARACTER*20 CNUM
      CHARACTER*2 NH
*
      DIMENSION indx(NNNW)
*
      POINTER (PCCMIN,ICCMIN(1))
      POINTER (PCDAMP,CDAMP(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/COUN/THRESH
     :      /DAMP/ODAMP(NNNW),PCDAMP
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEF7/PCCMIN,NCMIN,NCMAX
     :      /DEF9/CVAC,PI
     :      /DEFAULT/NDEF
     :      /FIXD/NFIX,LFIX(NNNW)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INVT/NOINVT(NNNW)
     :      /MCPB/DIAG,LFORDR
      COMMON/NODE/NNODEP(NNNW)
     :      /NPAR/PARM(2),NPARM
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB1/E(NNNW),GAMA(NNNW),PED(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORBA/IORDER(NNNW)
     :      /ORTHCT/ORTHST
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)

      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
*
*   Open, check, load data from, and close the  .iso  file
*
      CALL setiso (isofile) ! NOW ALSO GETS THE RADIAL PARAMETERS, Jon Grumer (Lund, 2013)
*
*   Set default speed of light and grid parameters
*
      C = CVAC
!      IF (NPARM .EQ. 0) THEN
!         RNT = EXP (-65.D0/16.D0) / Z
!         H = 0.5D0**4
!         N = MIN (220, NNNP)
!      ELSE
!         default comes here  ! DEFAULT IS SET ABOVE IN SETISO (from isodata), JG (Lund, 2013)
!         RNT = 2.D-6
!         H = 5.D-2
!         N = NNNP
!      ENDIF
!      HP = 0.D0
*
*     Print grid parameters to screen
*     
C      PRINT*,
C      PRINT*, 'DEFAULT GRID PARAMETERS'
C      PRINT*, '----------------------------------------------------'
C      PRINT*, 'RNT =', RNT 
C      PRINT*, 'H   =', H 
C      PRINT*, 'HP  =', HP
C      PRINT*, 'N   =', N
C      PRINT*, '----------------------------------------------------'
C      PRINT*,
cb
cb  default ACCY
*   ACCY is an estimate of the accuracy of the numerical procedures
*
      ACCY = H**6
cb
      IF (NDEF .NE. 0) THEN

         WRITE (istde,'(A)',ADVANCE='NO') 'Change the default speed of'
     &         ,' light or radial grid parameters?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
c
cbieron - shorter questions
cb
cb           WRITE (istde,*) 'The physical speed of light in'
cb    &          , ' atomic units is',CVAC,';', ' revise this value?'
            WRITE (istde,*) 'Speed of light = ',CVAC,';', ' revise ?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE (istde,*) 'Enter the revised value:'
               READ *, C
            ENDIF
*
*   Determine the parameters controlling the radial grid
*
            WRITE (istde,*) 'The default radial grid parameters for '
     &,                    'this case are:'
            WRITE (istde,*) ' RNT = ',RNT
            WRITE (istde,*) ' H   = ',H
            WRITE (istde,*) ' HP  = ',HP
            WRITE (istde,*) ' N   = ',N
            WRITE (istde,*) ' revise these values?'
            YES = GETYN ()
            IF (YES) THEN
               WRITE (istde,*) 'Enter RNT:'
               READ *, RNT
               WRITE (istde,*) 'Enter H:'
               READ *, H
               WRITE (istde,*) 'Enter HP:'
               READ *, HP
               WRITE (istde,*) 'Enter N:'
               READ *, N
cb Revised grid
               WRITE (istde,*) 'Revised RNT = ', RNT
               WRITE (istde,*) 'Revised H   = ', H
               WRITE (istde,*) 'Revised HP  = ', HP
               WRITE (istde,*) 'Revised N   = ', N
            ENDIF
         ENDIF
*
cb
cb read ACCY on input
cb
         WRITE (istde,*) 'Revise the default ACCY = ', ACCY
         YES = GETYN ()
            IF (YES) THEN
               WRITE (istde,*) 'Enter ACCY:'
               READ *, ACCY
cb Revised ACCY
               WRITE (istde,*) 'Revised ACCY = ', ACCY
         ENDIF
      ENDIF
*
*   Set up the coefficients for the numerical procedures
*
      CALL SETQIC
*
*   Generate the radial grid and all associated arrays
*
      CALL RADGRD
*
*   Generate $- r \times V_ (r)$
*
      CALL NUCPOT
*
*   Load the subshell radial wavefunction estimates
*
      CALL setrwfa (rwffile)
*
*   Set some defaults
*
      THRESH = 0.05D0

      DO I = 1, NW
!         IORDER(I) = I    ! Completely determined in GETOLD
         METHOD(I) = 3
         NOINVT(I) = .TRUE.
!CFF 2009    ODAMP(I) = 1.D0
         ODAMP(I) = 0.D0
         SCNSTY(I) = 0.D0
         PED(I)    = 0.D0
         IF (NAK(I) .LT. 0) THEN
            NNODEP(I) = NP(I) + NAK(I)
         ELSE
            NNODEP(I) = NP(I) - NAK(I) - 1
         ENDIF
      ENDDO

      IF (DIAG) THEN
         EOL = .FALSE.
         CALL GETALD          ! (E)AL type calculation,
                              ! H(DC) will not be diagonalised
      ELSEIF (LFORDR) THEN
         EOL = .TRUE.
         CALL GETOLD (idblk)  ! EOL type calculation
      ELSE
C         WRITE (istde,'(A)',ADVANCE='NO')
C     &            'EOL type calculation?  (y/n) '
         EOL = .true.
C         EOL = GETYN ()
         IF (EOL) THEN
C            IF (NDEF.EQ.0) THEN
C              WRITE(734,'(A)') 'y            ! EOL calculation'
C            END IF
            CALL GETOLD (idblk)
         ELSE
C            IF (NDEF.EQ.0) THEN
C              WRITE(734,'(A)') 'n            ! EOL calculation'
C            END IF
            CALL GETALD    ! (E)AL type calculation,
                              ! H(DC) will not be diagonalised
         ENDIF
      ENDIF

      WRITE (istde,*) 'Enter the maximum number of SCF cycles:'
      READ (*,*) NSCF

      IF (NDEF.EQ.0) THEN
         WRITE(734,*) NSCF,'! Number of SCF cycles'
      END IF

*
*   Allow the user to modify other defaults
*
      IF (NDEF .NE. 0) THEN
         WRITE (istde,'(A)',ADVANCE='NO') 
     &          'Modify other defaults?  (y/n) '
         YES = GETYN ()
      ELSE
         YES = .FALSE.
      ENDIF

      IF (.NOT. YES) RETURN
!=======================================================================
! From here to end, "other defaults" are handled. For simplicity
! We'll let node-0 do the job and then broadcast results to all
! nodes.
!=======================================================================
                  !-------------------------------------------
                  IF (myid .EQ. 0) THEN   ! This is a _big_ IF
                  !-------------------------------------------
*
*   THRESH
*
      WRITE (istde,*) 'An oscillation in the large-component of the '
     &,              'radial wavefunction is diregarded'
      WRITE (istde,*) 'for the purposes of node counting if its '
     &,              'amplitude is less than 1/20 the'
      WRITE (istde,*) 'maximum amplitude.   Revise this?'
      YES = GETYN ()
      IF (YES) THEN
    3    WRITE (istde,*) 'Enter the new threshold value:'
         READ (*,*) THRESH
         IF (THRESH .LE. 0.D0) THEN
            WRITE (istde,*) MYNAME, ': This must exceed 0;'
            GOTO 3
         ENDIF
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*
*   IORDER - has been determined at the place where the varying 
*            orbitals are specified
*
!      WRITE (istde,*) 'The subshells will be improved in the order'
!      CALL PRTRSL
!      WRITE (istde,*) 'Revise this order?'
!      YES = GETYN ()
      YES = .FALSE.
!      IF (YES) THEN
!         WRITE (istde,*) 'Revised order:'
!    4    CALL GETRSL (indx,NSUBS)
!         NORDER = 0
!         DO 6 I = 1,NSUBS
!            NPLOC = NP(indx(I))
!            NAKLOC = NAK(indx(I))
!            DO 5 J = 1,NW
!               IF ((NP(J) .EQ. NPLOC) .AND.
!     :             (NAK(J) .EQ. NAKLOC) .AND.
!     :             (.NOT. LFIX(J))) THEN
!                  NORDER = NORDER+1
!               ENDIF
!    5       CONTINUE
!    6    CONTINUE
!         IF (NORDER .NE. NW-NFIX) THEN
!            WRITE (istde,*) MYNAME, ': All subshells that are to be '
!     &,                    'improved must appear in the list.'
!            GOTO 4
!         ELSE
!            NORDER = 0
!            DO 7 I = 1,NW
!               IF (.NOT. LFIX(I)) THEN
!                  NORDER = NORDER+1
!                  IORDER(I) = indx(NORDER)
!               ENDIF
!    7       CONTINUE
!         ENDIF
!      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*
*   METHOD
*

! Piece only for printing...

      DO I = 1, 4

         DO J = 1, NW
            IF ((METHOD(J) .EQ. I) .AND. (.NOT. LFIX(J))) THEN
               WRITE (istde,*) 'Method ',I,' is used for '
     &,                       'integrating the radial differential '
     &,                       'equation for subshells'
               GOTO 9
            ENDIF
         ENDDO

         CYCLE

    9    IEND = 0
         DO J = 1, NW
            IF ((METHOD(J) .EQ. I) .AND. (.NOT. LFIX(J))) THEN
               IBEG = IEND + 1
               IEND = IBEG
               RECORD(IBEG:IEND) = ' '
               CALL CONVRT(NP(J), CNUM, LENTH)
               IBEG = IEND + 1
               IEND = IBEG + LENTH - 1
               RECORD(IBEG:IEND) = CNUM(1:LENTH)
               IBEG = IEND + 1
               IF (NAK(J) .LT. 0) THEN
                  IEND = IBEG
                  RECORD(IBEG:IEND) = NH(J)(1:1)
               ELSE
                  IEND = IBEG + 1
                  RECORD(IBEG:IEND) = NH(J)(1:2)
               ENDIF
            ENDIF
            IF (IEND .GT. 76) THEN
               WRITE (istde,*) RECORD(1:IEND)
               IEND = 0
            ENDIF
         ENDDO
         IF (IEND .GT. 0 .AND. myid .EQ. 0)
     &      WRITE (istde,*) RECORD(1:IEND)
      ENDDO

! Reads user inputs and fills array indx(1:nsubs) where nsubs itself
! is an output from getrsl.
! METHOD(1:4) is the only output. indx() and nsubs are discarded

      WRITE (istde,*) 'Select a different integration method for '
     &,              'any subshell radial wavefunction?'
      YES = GETYN ()
      IF (YES) THEN
         DO I = 1, 4
            WRITE (istde,*) 'Method ',I,':'
            CALL GETRSL (indx, NSUBS)
            DO J = 1, NSUBS
               LOC = indx(J)
               IF (.NOT. LFIX(LOC)) METHOD(LOC) = I
            ENDDO
         ENDDO
      ENDIF
*
*   NOINVT
*
      WRITE (istde,*) 'The first oscillation of the large component'

      DO I = 1, NW
         IF (NOINVT(I) .AND. (.NOT. LFIX(I))) THEN
            WRITE (istde,*) 'of the following radial wavefunctions '
     &,                    'will be required to be positive'
            GOTO 15
         ENDIF
      ENDDO

      WRITE (istde,*) 'of all radial wavefunctions will be required '
     &,              'to be positive.   Revise this?'
      YES = GETYN ()
      GOTO 17
   15 IEND = 0
      DO I = 1, NW
         IF (NOINVT(I) .AND. (.NOT. LFIX(I))) THEN
            IBEG = IEND + 1
            IEND = IBEG
            RECORD(IBEG:IEND) = ' '
            CALL CONVRT(NP(I),CNUM,LENTH)
            IBEG = IEND + 1
            IEND = IBEG + LENTH - 1
            RECORD(IBEG:IEND) = CNUM(1:LENTH)
            IBEG = IEND + 1
            IF (NAK(I) .LT. 0) THEN
               IEND = IBEG
               RECORD(IBEG:IEND) = NH(I)(1:1)
            ELSE
               IEND = IBEG + 1
               RECORD(IBEG:IEND) = NH(I)(1:2)
            ENDIF
         ENDIF
         IF (IEND .GT. 76) THEN
            WRITE (istde,*) RECORD(1:IEND)
            IEND = 0
         ENDIF
      ENDDO
      IF (IEND .GT. 0) WRITE (istde,*) RECORD(1:IEND)
      WRITE (istde,*) 'Revise this?'
      YES = GETYN ()
   17 IF (YES) THEN
         WRITE (istde,*) 'Suppressing enforcement of positive first '
     &,                 'oscillation:'
         CALL GETRSL (indx, NSUBS)
         DO I = 1, NSUBS
            LOC = indx(I)
            IF (.NOT. LFIX(LOC)) NOINVT(LOC) = .TRUE.
         ENDDO
      ENDIF
*
*   ODAMP
*
      DO I = 1, NW
         IF ((ODAMP(I) .NE. 0.D0) .AND.
     :       (.NOT. LFIX(I))) THEN
            WRITE (istde,*) 'Subshell accelerating parameters have '
     &,                    'been set.   Revise these?'
            YES = GETYN ()
            GOTO 20
         ENDIF
      ENDDO
      WRITE (istde,*) 'Set accelerating parameters for subshell '
     &,              'radial wavefunctions?'
      YES = GETYN ()
   20 IF (YES) THEN
         WRITE (istde,*) 'Different accelerating parameters for '
     &,                 'different subshell radial wavefunction?'
         YES = GETYN ()
         IF (YES) THEN
   21       WRITE (istde,*) 'Enter an accelerating parameter'
            WRITE (istde,*) ' (0< ODAMP < 1 allows ODAMP to be '
     &,                    'reduced as convergence is approached;'
            WRITE (istde,*) ' -1 < ODAMP < 0 implies |ODAMP| is '
     &,                    'held constant):'
            READ (*,*) ODAMPU
            IF ((ABS (ODAMPU) .EQ. 0.D0) .OR.
     :          (ABS (ODAMPU) .GE. 1.D0)) THEN
               WRITE (istde,*) MYNAME, ': Value out of range ...'
               GOTO 21
            ELSE
               CALL GETRSL (indx, NSUBS)
               DO I = 1, NSUBS
                  LOC = indx(I)
                  IF (.NOT. LFIX(LOC)) ODAMP(LOC) = ODAMPU
               ENDDO
            ENDIF
         ELSE
   23       WRITE (istde,*) 'Enter the accelerating parameter'
            WRITE (istde,*) ' (0< ODAMP < 1 allows ODAMP to be '
     &,                    'reduced as convergence is approached;'
            WRITE (istde,*) ' -1 < ODAMP < 0 implies |ODAMP| is '
     &,                    'held constant):'
            READ (*,*) ODAMPU
            IF ((ABS (ODAMPU) .EQ. 0.D0) .OR.
     :          (ABS (ODAMPU) .GE. 1.D0)) THEN
               WRITE (istde,*) MYNAME, ': Value out of range ...'
               GOTO 23
            ELSE
               DO I = 1, NW
                  IF (.NOT. LFIX(I)) ODAMP(I) = ODAMPU
               ENDDO
            ENDIF
         ENDIF
      ENDIF
*
*   CDAMP
*
      WRITE (istde,*) 'Set accelerating parameters for the '
     &,              'eigenvectors?'
      YES = GETYN ()
      IF (YES) THEN
         WRITE (istde,*) 'Different accelerating parameters '
     &,                 'for each eigenvector?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE (istde,*) 'Enter an accelerating parameter for'
            CALL CONVRT (NCMIN, RECORD, LENTH)
            WRITE (istde,*) ' each of the '//RECORD(1:LENTH)//
     &                     ' levels :'
            READ (*,*) (CDAMP(I), I = 1, NCMIN)
         ELSE
            WRITE (istde,*) 'Enter the accelerating parameter:'
            READ (*,*) CDAMPU
            DO I = 1, NCMIN
               CDAMP(I) = CDAMPU
            ENDDO
         ENDIF
      ENDIF
*
*   NSIC
*
      WRITE (istde,*) 'Following the improvement of each of the '
     &,              'subshell radial wavefunctions in turn, '
      WRITE (istde,*) 'the ',NSIC,' least self-consistent'
     &,              ' functions will be improved at the'
      WRITE (istde,*) 'end of the first SCF cycle. Revise this '
     &,              'setting?'
      YES = GETYN ()
      IF (YES) THEN
         WRITE (istde,*) 'Enter the number of additional '
     &,                 'improvements:'
         READ (*,*) NSIC
      ENDIF
*
*   NSOLV
*
      WRITE (istde,*) 'The maximum number of cycles in attempting '
     &,              'to solve each radial equation is '
      WRITE (istde,*) NSOLV,' times the principal quantum'
     &,             ' number of the radial'
      WRITE (istde,*) 'wave-function to be estimated.   '
     &,              'Revise this setting?'
      YES = GETYN ()
      IF (YES) THEN
         WRITE (istde,*) 'Enter the factor that multiplies the '
     &,                 'principal quantum number:'
         READ (*,*) NSOLV
      ENDIF
*
*   Orthogonalisation
*
      IF (ORTHST) THEN
         WRITE (istde,*) 'Subshell radial wavefunctions will be '
     &,                 'Schmidt orthogonalised immediately'
         WRITE (istde,*) 'following their estimation to all '
     &,                 'functions with poorer self-consistency.'
         WRITE (istde,*) ' Revise this?'
         YES = GETYN ()
         IF (YES) ORTHST = .FALSE.
      ELSE
         WRITE (istde,*) 'Subshell radial wavefunctions will be '
     &,                 'Schmidt orthogonalised at the end of'
         WRITE (istde,*) 'each SCF cycle.   Revise this?'
         YES = GETYN ()
         IF (YES) ORTHST = .TRUE.
      ENDIF
                  !-------------------------------------------
                  ENDIF ! end of the _big_ IF
                  !-------------------------------------------
      RETURN
      END
