************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***          *****   *******  **   **  ****   *****    *****         ***
***         **   **  **       ***  **   **   **   **  **   **        ***
***         **       **       ** * **   **   **       **   **        ***
***         **  ***  ****     ** * **   **    *****   **   **        ***
***         **   **  **       ** * **   **        **  **   **        ***
***         **   **  **       **  ***   **   **   **  **   **        ***
***          *****   *******  **   **  ****   *****    *****         ***
***                                                                  ***
***   Package to generate the nuclear charge, geometry, mass, spin,  ***
***   and electromagnetic moment data file for the GRASP92 codes.    ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM GENISO
*                                                                      *
*   Generates the isotope data file for the GRASP92 suite of codes.    *
*                                                                      *
*   Call(s) to: GETCPR, GETYN, LENGTH, OPENFL.                         *
*                                                                      *
*   Written by Farid A. Parpia.           Last revision: 16 Oct 1994   *
*                                                                      *
*   Update: 2014-02-13 - Jon Grumer, Lund University, Sweden           *
*           Grid parameters are written to isodata so other routines   *
*           may read this file instead of asking the user (non-default *
*           grid parameters).                                          *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      LOGICAL GETYN,YES
*     
      ! JG (Lund, 2013) ---------------------------------------------------------
      INCLUDE 'parameters.def'       ! Grasp2k parameter file (src/lib/def/)
      CHARACTER*1 GRIDANSW,GRIDANSW2 ! Answers to the default grid questions
      DOUBLE PRECISION RNT, H, HP    ! Grid parameters to be included in isodata
      INTEGER N                      ! ...
      ! -------------------------------------------------------------------------
*
!XHH To be more implementation-independent, add 3 lines to each main 
!XHH program
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde

cb alpha constant from lib/lib92/setcon.f
cb AUMAMU from lib/lib92/setcon.f
cb
cb    DATA EMEAMU /5.48579903D-04/
cb   :     ALFAI  /137.0359895D 00/
      COMMON/DEF2/C
     :      /DEF11/FMTOAU,AUMAMU
cb
      ALFAI = 1.0D 00/C
cb EMEAMU: Electron mass in amu
      AUMAMU = EMEAMU
cb end alpha constant

      write(*,*)
      write(*,*) 'RNUCLEUS'
      write(*,*) 'This program defines nuclear data and the radial grid'
      write(*,*) 'Outputfile: isodata'
      write(*,*)
*
*   File  grasp92.iso  is FORMATTED
*
      DEFNAM = 'isodata'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      FILNAM = DEFNAM
*
      CALL OPENFL (22,FILNAM,FORM,STATUS,IERR)
*
      IF (IERR .NE. 0) THEN
         WRITE(istde,*) 'Error when opening isodata'
         STOP
      ENDIF
*
      WRITE(istde,*) 'Enter the atomic number:'

      READ *, Z
      WRITE (22,300) 'Atomic number:'
      WRITE (22,*) Z
*
      WRITE(istde,*) 'Enter the mass number (0 if the'
     &, ' nucleus is to be modelled as a point source:'

      READ *, A

      WRITE (22,300) 'Mass number (integer) :'
      WRITE (22,*) A

      IF (A .EQ. 0.0D 00) THEN

         CPARM = 0.0D 00
         APARM = 0.0D 00
      ELSE
         RRMS =  0.836D 00*A**(1.0D 00/3.0D 00)
     :          +0.570D 00
         WRITE(istde,*) 'The default root mean squared'
     &, ' radius is ',RRMS,' fm;'
         TPARM = 2.30D 00
         WRITE(istde,*) 'The default nuclear skin thickness'
     &, ' is   ',TPARM,' fm;'
         WRITE(istde,*) 'Revise these values?'
         YES = GETYN ()
         IF (YES) THEN
            WRITE(istde,*) 'Enter the root mean squared'
     &, ' radius of the nucleus (in fm):'
            READ *, RRMS
            WRITE(istde,*) 'Enter the skin thickness of'
     &, ' the nucleus (in fm):'
            READ *, TPARM
         ENDIF
         APARM = TPARM/(4.0D 00*LOG (3.0D 00))
         CALL GETCPR (RRMS,APARM,CPARM)
      ENDIF
      WRITE (22,300) 'Fermi distribution parameter a:'
      WRITE (22,*) APARM
      WRITE (22,300) 'Fermi distribution parameter c:'
      WRITE (22,*) CPARM
*
      WRITE(istde,*) 'Enter the mass of the neutral'
     &, ' atom (in amu) (0 if the nucleus is to be static):'
      READ *, AMAMU
      IF (AMAMU .NE. 0.0D 00) THEN
!         WRITE(istde,*) 'Enter your best estimate of the ground'
!     &, ' state energy of, the neutral atom (in Hartrees):'
!         READ *, EBIND
          EBIND = 0.D0

!XHH better use NINT, not INT. 
!         NENEU = INT (Z)
         NENEU = NINT (Z)

!         WRITE(istde,*) 'The number of electrons in the'
!     &, ' neutral atom is deduced to be',NENEU,';'
!         WRITE(istde,*) 'Revise this?'
!         YES = GETYN ()
!         IF (YES) THEN
!            WRITE(istde,*) 'Enter the number of electrons'
!     &, ' in the neutral atom:'
!            READ *, NENEU
!         ENDIF
         IF (EBIND .GT. 0.0D 00) EBIND = -EBIND
         EMNAMU = AMAMU-EMEAMU*DBLE (NENEU)
     :                 -EMEAMU*EBIND/ALFAI**2
      ELSE
         EMNAMU = 0.0D 00
      ENDIF
*
      WRITE (22,300) 'Mass of nucleus (in amu):'
      WRITE (22,*) EMNAMU
*
      WRITE(istde,*) 'Enter the nuclear spin quantum'
     &, ' number (I) (in units of h / 2 pi):'
      READ *, SQN
      WRITE (22,300) 'Nuclear spin (I) (in units of h / 2 pi):'
      WRITE (22,*) SQN
*
      WRITE(istde,*) 'Enter the nuclear dipole moment'
     &, ' (in nuclear magnetons):'
      READ *, DMOMNM
      WRITE (22,300) 'Nuclear dipole moment (in nuclear magnetons):'
      WRITE (22,*) DMOMNM
*
      WRITE(istde,*) 'Enter the nuclear quadrupole'
     &, ' moment (in barns):'
      READ *, QMOMB
      WRITE (22,300) 'Nuclear quadrupole moment (in barns):'
      WRITE (22,*) QMOMB
*     
*     Grid Parameters - Jon Grumer (Lund, 2013)
*      
*     Default values
*
*     If chosen to model nucleus as point source then
      IF (A .EQ. 0) THEN 
         WRITE(istde,*) 'You have chosen to model the nucleus as '
     &,                 'a point source!'
         RNT = EXP (-65.0D 00/16.0D 00) / real(Z)
         H   = 0.5D 00**4
         HP  = 0.d0
         N   = MIN (220,NNNP)
*     Otherwize set normal default values         
      ELSE
         RNT = 2.d-6
         H   = 5.d-2
         HP  = 0.d0
         N   = INT(NNNP) ! default number of grid points from parameters.def
      END IF
*
*     Print default values to screen
*
  101 WRITE(istde,*) '------------------------------------------------'
     &,              '-----------'
      WRITE(istde,*) 'The Grasp2K grid is defined by:'
      WRITE(istde,*) 
      WRITE(istde,*) 'R(i) = RNT*[ exp[ (i-1)*H ] - 1 ]'
      WRITE(istde,*) '  i  = 1, 2, ..., NNNP'
      WRITE(istde,*) 
      WRITE(istde,*) 'The default grid parameters are:'
      WRITE(istde,*) 
      WRITE(istde,*) 'RNT  (first grid point       ) = ', RNT
      WRITE(istde,*) 'H    (grid step-size         ) = ', H
!      WRITE(istde,*) 'HP   (related to linear grid ) = ', HP
      WRITE(istde,*) 'NNNP (max. no. of grid-points) = ', N
      WRITE(istde,*) '------------------------------------------------'
     &,              '-----------'
      WRITE(istde,*) 
      WRITE(istde,*) 'Do you want to revise these values (y/*)?'
*
      READ (*,'(A)') GRIDANSW
*
      IF (GRIDANSW.EQ.'y'.OR.GRIDANSW.EQ.'Y') THEN
*
         WRITE(istde,*)
         WRITE(istde,*) 'NNNP depends on H. A new value of'
         WRITE(istde,*) 'NNNP should satisfy the inequality       '
         WRITE(istde,*) 'NNNP_new >= (0.05/H_new)*590 '
         WRITE(istde,*)


         WRITE(istde,*) 'Do you want to revise RNT (y/*)?'
         READ (*,'(A)'), GRIDANSW2
         IF (GRIDANSW2.EQ.'y'.OR.GRIDANSW2.EQ.'Y') THEN
            WRITE(istde,*) 'Enter RNT (double precision real):'
            READ *, RNT
         ELSE
            WRITE(istde,*) 'Default value of RNT is set'
         END IF
*      
         WRITE(istde,*) 'Do you want to revise H (y/*)?'
         READ (*,'(A)'), GRIDANSW2
         IF (GRIDANSW2.EQ.'y'.OR.GRIDANSW2.EQ.'Y') THEN
            WRITE(istde,*) 'Enter H (double precision real):'
            READ *, H
         ELSE
            WRITE(istde,*) 'Default value of H is set'
         END IF
*      
!          WRITE(istde,*) 'Do you want to revise HP (y/*)?'
!          READ (*,'(A)'), GRIDANSW2
!          IF (GRIDANSW2.EQ.'y'.OR.GRIDANSW2.EQ.'Y') THEN
!             WRITE(istde,*) 'Enter HP (double precision real):'
!             READ *, HP
!          ELSE
!             WRITE(istde,*) 'Default value of HP is set'
!          END IF
*      
         WRITE(istde,*) 'Do you want to revise NNNP (y/*)?'
         READ (*,'(A)'), GRIDANSW2
         IF (GRIDANSW2.EQ.'y'.OR.GRIDANSW2.EQ.'Y') THEN
            WRITE(istde,'(a,i6,a)') 'Enter NNNP (integer <=', NNNP,' ):'
            READ *, N
            IF (N.le.(0.05/H)*590) THEN
               WRITE(*,*) 'Relation between NNNP and H not fullfilled'
               STOP
            END IF
         ELSE
            WRITE(istde,*) 'Default value of NNNP is set'
         END IF
*
      ELSE 
         WRITE(istde,*) 'Default grid parameter values are kept!' 
      END IF
*
*     Write grid parameters to isodata
*
      WRITE (22,300) 'RNT:'
      WRITE (22,*) RNT
*
      WRITE (22,300) 'H:'
      WRITE (22,*) H
*
      WRITE (22,300) 'HP:'
      WRITE (22,*) HP
*
      WRITE (22,300) 'NNNP:'
      WRITE (22,*) N
*
      CLOSE (22)
*
      WRITE(istde,*)
      WRITE(istde,*) 'Output file isodata successfully created. '
     &,              'Happy computing!'
      WRITE(istde,*)

      STOP
*
  300 FORMAT (A)
*
      END
