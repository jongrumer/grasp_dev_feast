!
!***********************************************************************
!                                                                      *
      MODULE blk_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  14:35:02   1/ 6/07  
!GG      INTEGER, DIMENSION(10) :: NCFBLK 
      INTEGER, DIMENSION(20) :: NCFBLK, NEVINBLK, NCFINBLK, TWO_J
      INTEGER :: NBLOCK 
      REAL(DOUBLE), DIMENSION(:), pointer :: IDXBLK 
      INTEGER :: NELECTOT, NCFTOT, NWTOT, NVECTOT, NVECSIZTOT, NBLOCK1 
      INTEGER :: NBLOCKI 
      INTEGER :: NBLOCKF 
      INTEGER, DIMENSION(10) :: NCFI 
      INTEGER, DIMENSION(10) :: NCFF 
      END MODULE blk_C 
!
!***********************************************************************
!                                                                      *
      MODULE cons_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:42:40   1/ 2/07  
      REAL(DOUBLE) :: ZERO = 0.0D00, &
                      HALF = 0.5D00, &
                      TENTH= 0.1D00, &
                      ONE  = 1.0D00, &
                      TWO  = 2.0D00, &
                      THREE= 3.0D00, &
                      FOUR = 4.0D00, &
                      SEVEN= 7.0D00, &
                      TEN  =10.0D00, &
                      EPS  = 1.0D-10 
      END MODULE cons_C 
!
!***********************************************************************
!                                                                      *
      MODULE debug_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  14:35:02   1/ 6/07  
      INTEGER :: IBUG1, IBUG2, IBUG3, IBUG4, IBUG5, IBUG6 
      LOGICAL, DIMENSION(5) :: LDBPA 
      LOGICAL, DIMENSION(5) :: LDBPG 
      LOGICAL, DIMENSION(30) :: LDBPR 
      REAL(DOUBLE) :: cutoff  ! used by bioscl
      END MODULE debug_C 
!
!***********************************************************************
!                                                                      *
      MODULE def_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  11:02:52   1/ 2/07  
      REAL(DOUBLE) :: TENMAX, EXPMAX, EXPMIN, PRECIS 
      REAL(DOUBLE) :: AUCM, AUEV, CCMS, FASI, FBSI 
      REAL(DOUBLE) :: FMTOAU, AUMAMU, B1
      INTEGER :: IONCTY, NELEC 
      REAL(DOUBLE) :: EMN, Z 
      INTEGER :: IONCTYFF, NELECFF 
      REAL(DOUBLE) :: EMNFF, ZFF 
      INTEGER :: IONCTYII, NELECII 
      REAL(DOUBLE) :: EMNII, ZII 
      INTEGER :: NELECR 
      REAL(DOUBLE) :: C 
      REAL(DOUBLE) :: EMPAM, RBCM 
      INTEGER :: NSCF, NSIC, NSOLV 
      REAL(DOUBLE) :: ACCY 
!     REAL(DOUBLE) :: PNTRWT, PWEIGH 
      REAL(DOUBLE), DIMENSION(:), pointer :: wt, weight
      INTEGER :: NCMIN, NCMAX 
!     REAL(DOUBLE) :: PCCMIN 
      INTEGER, DIMENSION(:), pointer :: iccmin
      REAL(DOUBLE) :: CVAC, PI 
      INTEGER, PARAMETER :: NNNP = 590 
      REAL(DOUBLE), DIMENSION(NNNP) :: DP, DQ 
      END MODULE def_C 
!
!***********************************************************************
!                                                                      *
      MODULE eigv_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07  
      REAL(DOUBLE) :: EAV, EAVFF, EAVII
      REAL(DOUBLE), DIMENSION(:), pointer :: EVAL, EVALFF, EVALII
      REAL(DOUBLE), DIMENSION(:), pointer :: EVEC, EVECFF, EVECII
      END MODULE eigv_C 
!
!***********************************************************************
!                                                                      *
      MODULE facts_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:46:52   2/14/04  
      INTEGER, PARAMETER :: MFACT = 500 
      REAL(DOUBLE), DIMENSION(MFACT) :: GAM 
      END MODULE facts_C 
!
!***********************************************************************
!                                                                      *
      MODULE iounit_C 
!                                                                      *
!***********************************************************************
!...Created by Pacific-Sierra Research 77to90  4.3E  06:27:59  12/28/06  
      INTEGER :: ISTDI=5, ISTDO=6, ISTDE=0 
      END MODULE iounit_C 
!
!***********************************************************************
!                                                                      *
      MODULE jj2lsj_C 
!                                                                      *
!     This module contains the (numerical) values of the subshell      *
!     terms in LS-coupling and in it are define some global variables  *
!     for jj2LSJ program.                                              *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  LONG, DOUBLE 
!
! Define some global data for the jj2LSJ transformation program
!
      type :: subshell_term_LS
         integer :: l_shell ! angular momentum l
         integer :: w       ! w
         integer :: Q       ! Subshell total quasispin 2*Q
         integer :: LL      ! Subshell total angular momentum 2*L
         integer :: S       ! Subshell Total sinp 2*S
      end type subshell_term_LS
!
      type :: LS_jj_me
!         sequence
         integer       :: w, Q, L, S, J, Nm, Qm, Jm, Qp, Jp
         integer       :: factor
         integer(LONG) :: nom, denom
      end type LS_jj_me
!
! Define the values of all LS-jj transformation coefficients
!
      integer, parameter :: LS_jj_number_p3 =  10, LS_jj_number_p4  =    9, &
                            LS_jj_number_p5 =   2, LS_jj_number_p6  =    1, &
                            LS_jj_number_d3 =  65, LS_jj_number_d4  =  166, &
                            LS_jj_number_d5 = 184, LS_jj_number_d6  =  166, &
                            LS_jj_number_d7 =  65, LS_jj_number_d8  =   19, &
                            LS_jj_number_d9 =   2, LS_jj_number_d10 =    1, &
                            LS_jj_number_f3 = 216, LS_jj_number_f4  = 1210, &
                            LS_jj_number_f5 =3799, LS_jj_number_f6  = 7313, &
                            LS_jj_number_f7 =8003
!
! Define all possible antisymetric subshell terms in LS-coupling
!
   type(subshell_term_LS), dimension(1:2), parameter ::   &
      term_LS_s =(/                                       &
      subshell_term_LS(0, 0, 0, 0, 1),                    &
      subshell_term_LS(0, 0, 1, 0, 0)  /)
   type(subshell_term_LS), dimension(1:6), parameter ::   &
      term_LS_p =(/                                       &
      subshell_term_LS(1, 0, 0, 0, 3),                    &
      subshell_term_LS(1, 0, 2, 2, 1),                    &
      subshell_term_LS(1, 0, 0, 4, 1),                    &
      subshell_term_LS(1, 0, 3, 0, 0),                    &
      subshell_term_LS(1, 0, 1, 2, 2),                    &
      subshell_term_LS(1, 0, 1, 4, 0)  /)
   type(subshell_term_LS), dimension(1:32), parameter ::  &
      term_LS_d =(/                                       &
      subshell_term_LS(2, 0, 0, 0, 5),                    &
      subshell_term_LS(2, 0, 0, 0, 1),                    &
      subshell_term_LS(2, 0, 2, 2, 3),                    &
      subshell_term_LS(2, 0, 2, 2, 1),                    &
      subshell_term_LS(2, 0, 4, 4, 1),                    &
      subshell_term_LS(2, 0, 2, 4, 1),                    &
      subshell_term_LS(2, 0, 0, 4, 3),                    &
      subshell_term_LS(2, 0, 0, 4, 1),                    &
      subshell_term_LS(2, 0, 2, 6, 3),                    &
      subshell_term_LS(2, 0, 2, 6, 1),                    &
      subshell_term_LS(2, 0, 0, 6, 1),                    &
      subshell_term_LS(2, 0, 2, 8, 1),                    &
      subshell_term_LS(2, 0, 0, 8, 3),                    &
      subshell_term_LS(2, 0, 0, 8, 1),                    &
      subshell_term_LS(2, 0, 2,10, 1),                    &
      subshell_term_LS(2, 0, 0,12, 1),                    &
      subshell_term_LS(2, 0, 5, 0, 0),                    &
      subshell_term_LS(2, 0, 1, 0, 0),                    &
      subshell_term_LS(2, 0, 3, 2, 2),                    &
      subshell_term_LS(2, 0, 1, 2, 2),                    &
      subshell_term_LS(2, 0, 1, 4, 4),                    &
      subshell_term_LS(2, 0, 1, 4, 2),                    &
      subshell_term_LS(2, 0, 3, 4, 0),                    &
      subshell_term_LS(2, 0, 1, 4, 0),                    &
      subshell_term_LS(2, 0, 3, 6, 2),                    &
      subshell_term_LS(2, 0, 1, 6, 2),                    &
      subshell_term_LS(2, 0, 1, 6, 0),                    &
      subshell_term_LS(2, 0, 1, 8, 2),                    &
      subshell_term_LS(2, 0, 3, 8, 0),                    &
      subshell_term_LS(2, 0, 1, 8, 0),                    &
      subshell_term_LS(2, 0, 1,10, 2),                    &
      subshell_term_LS(2, 0, 1,12, 0)  /)
   type(subshell_term_LS), dimension(1:238), parameter :: &
      term_LS_f =(/                                       &
      subshell_term_LS(3, 0, 0, 0, 7),                    &
      subshell_term_LS(3, 0, 2, 2, 5),                    &
      subshell_term_LS(3, 0, 0, 4, 5),                    &
      subshell_term_LS(3, 0, 2, 6, 5),                    &
      subshell_term_LS(3, 0, 0, 8, 5),                    &
      subshell_term_LS(3, 0, 2,10, 5),                    &
      subshell_term_LS(3, 0, 0,12, 5),                    &
      subshell_term_LS(3, 1, 4, 0, 3),                    &
      subshell_term_LS(3, 2, 0, 0, 3),                    &
      subshell_term_LS(3, 1, 2, 2, 3),                    &
      subshell_term_LS(3, 2, 2, 2, 3),                    &
      subshell_term_LS(3, 1, 4, 4, 3),                    &
      subshell_term_LS(3, 2, 2, 4, 3),                    &
      subshell_term_LS(3, 3, 2, 4, 3),                    &
      subshell_term_LS(3, 4, 0, 4, 3),                    &
      subshell_term_LS(3, 5, 0, 4, 3),                    &
      subshell_term_LS(3, 6, 0, 4, 3),                    &
      subshell_term_LS(3, 1, 4, 6, 3),                    &
      subshell_term_LS(3, 2, 2, 6, 3),                    &
      subshell_term_LS(3, 3, 2, 6, 3),                    &
      subshell_term_LS(3, 4, 2, 6, 3),                    &
      subshell_term_LS(3, 5, 0, 6, 3),                    &
      subshell_term_LS(3, 1, 4, 8, 3),                    &
      subshell_term_LS(3, 2, 2, 8, 3),                    &
      subshell_term_LS(3, 3, 2, 8, 3),                    &
      subshell_term_LS(3, 4, 2, 8, 3),                    &
      subshell_term_LS(3, 5, 0, 8, 3),                    &
      subshell_term_LS(3, 6, 0, 8, 3),                    &
      subshell_term_LS(3, 7, 0, 8, 3),                    &
      subshell_term_LS(3, 1, 2,10, 3),                    &
      subshell_term_LS(3, 2, 2,10, 3),                    &
      subshell_term_LS(3, 3, 2,10, 3),                    &
      subshell_term_LS(3, 4, 0,10, 3),                    &
      subshell_term_LS(3, 5, 0,10, 3),                    &
      subshell_term_LS(3, 1, 4,12, 3),                    &
      subshell_term_LS(3, 2, 2,12, 3),                    &
      subshell_term_LS(3, 3, 2,12, 3),                    &
      subshell_term_LS(3, 4, 0,12, 3),                    &
      subshell_term_LS(3, 5, 0,12, 3),                    &
      subshell_term_LS(3, 1, 2,14, 3),                    &
      subshell_term_LS(3, 2, 2,14, 3),                    &
      subshell_term_LS(3, 3, 0,14, 3),                    &
      subshell_term_LS(3, 1, 2,16, 3),                    &
      subshell_term_LS(3, 2, 0,16, 3),                    &
      subshell_term_LS(3, 3, 0,16, 3),                    &
      subshell_term_LS(3, 0, 2,18, 3),                    &
      subshell_term_LS(3, 0, 0,20, 3),                    &
      subshell_term_LS(3, 1, 0, 0, 1),                    &
      subshell_term_LS(3, 2, 0, 0, 1),                    &
      subshell_term_LS(3, 1, 4, 2, 1),                    &
      subshell_term_LS(3, 2, 2, 2, 1),                    &
      subshell_term_LS(3, 3, 2, 2, 1),                    &
      subshell_term_LS(3, 4, 2, 2, 1),                    &
      subshell_term_LS(3, 5, 0, 2, 1),                    &
      subshell_term_LS(3, 1, 4, 4, 1),                    &
      subshell_term_LS(3, 2, 4, 4, 1),                    &
      subshell_term_LS(3, 3, 2, 4, 1),                    &
      subshell_term_LS(3, 4, 2, 4, 1),                    &
      subshell_term_LS(3, 5, 2, 4, 1),                    &
      subshell_term_LS(3, 6, 0, 4, 1),                    &
      subshell_term_LS(3, 7, 0, 4, 1),                    &
      subshell_term_LS(3, 1, 6, 6, 1),                    &
      subshell_term_LS(3, 2, 4, 6, 1),                    &
      subshell_term_LS(3, 3, 2, 6, 1),                    &
      subshell_term_LS(3, 4, 2, 6, 1),                    &
      subshell_term_LS(3, 5, 2, 6, 1),                    &
      subshell_term_LS(3, 6, 2, 6, 1),                    &
      subshell_term_LS(3, 7, 2, 6, 1),                    &
      subshell_term_LS(3, 8, 0, 6, 1),                    &
      subshell_term_LS(3, 9, 0, 6, 1),                    &
      subshell_term_LS(3,10, 0, 6, 1),                    &
      subshell_term_LS(3, 1, 4, 8, 1),                    &
      subshell_term_LS(3, 2, 4, 8, 1),                    &
      subshell_term_LS(3, 3, 2, 8, 1),                    &
      subshell_term_LS(3, 4, 2, 8, 1),                    &
      subshell_term_LS(3, 5, 2, 8, 1),                    &
      subshell_term_LS(3, 6, 2, 8, 1),                    &
      subshell_term_LS(3, 7, 0, 8, 1),                    &
      subshell_term_LS(3, 8, 0, 8, 1),                    &
      subshell_term_LS(3, 9, 0, 8, 1),                    &
      subshell_term_LS(3,10, 0, 8, 1),                    &
      subshell_term_LS(3, 1, 4,10, 1),                    &
      subshell_term_LS(3, 2, 4,10, 1),                    &
      subshell_term_LS(3, 3, 2,10, 1),                    &
      subshell_term_LS(3, 4, 2,10, 1),                    &
      subshell_term_LS(3, 5, 2,10, 1),                    &
      subshell_term_LS(3, 6, 2,10, 1),                    &
      subshell_term_LS(3, 7, 2,10, 1),                    &
      subshell_term_LS(3, 8, 0,10, 1),                    &
      subshell_term_LS(3, 9, 0,10, 1),                    &
      subshell_term_LS(3, 1, 4,12, 1),                    &
      subshell_term_LS(3, 2, 2,12, 1),                    &
      subshell_term_LS(3, 3, 2,12, 1),                    &
      subshell_term_LS(3, 4, 2,12, 1),                    &
      subshell_term_LS(3, 5, 2,12, 1),                    &
      subshell_term_LS(3, 6, 0,12, 1),                    &
      subshell_term_LS(3, 7, 0,12, 1),                    &
      subshell_term_LS(3, 8, 0,12, 1),                    &
      subshell_term_LS(3, 9, 0,12, 1),                    &
      subshell_term_LS(3, 1, 4,14, 1),                    &
      subshell_term_LS(3, 2, 2,14, 1),                    &
      subshell_term_LS(3, 3, 2,14, 1),                    &
      subshell_term_LS(3, 4, 2,14, 1),                    &
      subshell_term_LS(3, 5, 2,14, 1),                    &
      subshell_term_LS(3, 6, 0,14, 1),                    &
      subshell_term_LS(3, 7, 0,14, 1),                    &
      subshell_term_LS(3, 1, 4,16, 1),                    &
      subshell_term_LS(3, 2, 2,16, 1),                    &
      subshell_term_LS(3, 3, 2,16, 1),                    &
      subshell_term_LS(3, 4, 0,16, 1),                    &
      subshell_term_LS(3, 5, 0,16, 1),                    &
      subshell_term_LS(3, 1, 2,18, 1),                    &
      subshell_term_LS(3, 2, 2,18, 1),                    &
      subshell_term_LS(3, 3, 0,18, 1),                    &
      subshell_term_LS(3, 4, 0,18, 1),                    &
      subshell_term_LS(3, 1, 2,20, 1),                    &
      subshell_term_LS(3, 2, 0,20, 1),                    &
      subshell_term_LS(3, 0, 2,22, 1),                    &
      subshell_term_LS(3, 0, 0,24, 1),                    &
      subshell_term_LS(3, 0, 1, 6, 6),                    &
      subshell_term_LS(3, 1, 3, 4, 4),                    &
      subshell_term_LS(3, 2, 1, 4, 4),                    &
      subshell_term_LS(3, 3, 1, 4, 4),                    &
      subshell_term_LS(3, 1, 3, 6, 4),                    &
      subshell_term_LS(3, 2, 1, 6, 4),                    &
      subshell_term_LS(3, 1, 3, 8, 4),                    &
      subshell_term_LS(3, 2, 1, 8, 4),                    &
      subshell_term_LS(3, 3, 1, 8, 4),                    &
      subshell_term_LS(3, 0, 1, 2, 4),                    &
      subshell_term_LS(3, 1, 1,10, 4),                    &
      subshell_term_LS(3, 2, 1,10, 4),                    &
      subshell_term_LS(3, 0, 3, 0, 4),                    &
      subshell_term_LS(3, 1, 3,12, 4),                    &
      subshell_term_LS(3, 2, 1,12, 4),                    &
      subshell_term_LS(3, 0, 1,14, 4),                    &
      subshell_term_LS(3, 0, 1,16, 4),                    &
      subshell_term_LS(3, 1, 5, 6, 2),                    &
      subshell_term_LS(3, 2, 3, 6, 2),                    &
      subshell_term_LS(3, 6, 1, 6, 2),                    &
      subshell_term_LS(3, 8, 1, 6, 2),                    &
      subshell_term_LS(3, 1, 3, 4, 2),                    &
      subshell_term_LS(3, 2, 3, 4, 2),                    &
      subshell_term_LS(3, 3, 1, 4, 2),                    &
      subshell_term_LS(3, 4, 1, 4, 2),                    &
      subshell_term_LS(3, 3, 3, 6, 2),                    &
      subshell_term_LS(3, 5, 1, 6, 2),                    &
      subshell_term_LS(3, 1, 3, 8, 2),                    &
      subshell_term_LS(3, 2, 3, 8, 2),                    &
      subshell_term_LS(3, 4, 1, 8, 2),                    &
      subshell_term_LS(3, 5, 1, 8, 2),                    &
      subshell_term_LS(3, 5, 1, 4, 2),                    &
      subshell_term_LS(3, 4, 3, 6, 2),                    &
      subshell_term_LS(3, 7, 1, 6, 2),                    &
      subshell_term_LS(3, 9, 1, 6, 2),                    &
      subshell_term_LS(3, 3, 3, 8, 2),                    &
      subshell_term_LS(3, 6, 1, 8, 2),                    &
      subshell_term_LS(3, 7, 1, 8, 2),                    &
      subshell_term_LS(3, 1, 5, 2, 2),                    &
      subshell_term_LS(3, 2, 3, 2, 2),                    &
      subshell_term_LS(3, 3, 3, 2, 2),                    &
      subshell_term_LS(3, 1, 5,10, 2),                    &
      subshell_term_LS(3, 2, 3,10, 2),                    &
      subshell_term_LS(3, 3, 3,10, 2),                    &
      subshell_term_LS(3, 4, 3,10, 2),                    &
      subshell_term_LS(3, 4, 1, 2, 2),                    &
      subshell_term_LS(3, 5, 1,10, 2),                    &
      subshell_term_LS(3, 6, 1,10, 2),                    &
      subshell_term_LS(3, 5, 1, 2, 2),                    &
      subshell_term_LS(3, 6, 1, 2, 2),                    &
      subshell_term_LS(3, 7, 1,10, 2),                    &
      subshell_term_LS(3, 8, 1,10, 2),                    &
      subshell_term_LS(3, 9, 1,10, 2),                    &
      subshell_term_LS(3, 1, 3,12, 2),                    &
      subshell_term_LS(3, 2, 3,12, 2),                    &
      subshell_term_LS(3, 3, 1,12, 2),                    &
      subshell_term_LS(3, 3, 1,12, 2),                    &
      subshell_term_LS(3, 5, 1,12, 2),                    &
      subshell_term_LS(3, 6, 1,12, 2),                    &
      subshell_term_LS(3, 1, 3,14, 2),                    &
      subshell_term_LS(3, 2, 3,14, 2),                    &
      subshell_term_LS(3, 3, 1,14, 2),                    &
      subshell_term_LS(3, 4, 1,14, 2),                    &
      subshell_term_LS(3, 5, 1,14, 2),                    &
      subshell_term_LS(3, 6, 1,14, 2),                    &
      subshell_term_LS(3, 1, 3,16, 2),                    &
      subshell_term_LS(3, 2, 1,16, 2),                    &
      subshell_term_LS(3, 3, 1,16, 2),                    &
      subshell_term_LS(3, 1, 3,18, 2),                    &
      subshell_term_LS(3, 2, 1,18, 2),                    &
      subshell_term_LS(3, 3, 1,18, 2),                    &
      subshell_term_LS(3, 0, 1,20, 2),                    &
      subshell_term_LS(3, 0, 1,22, 2),                    &
      subshell_term_LS(3, 2, 1, 6, 0),                    &
      subshell_term_LS(3, 3, 1, 6, 0),                    &
      subshell_term_LS(3, 4, 1, 6, 0),                    &
      subshell_term_LS(3, 1, 5, 4, 0),                    &
      subshell_term_LS(3, 2, 3, 4, 0),                    &
      subshell_term_LS(3, 3, 3, 4, 0),                    &
      subshell_term_LS(3, 1, 3, 6, 0),                    &
      subshell_term_LS(3, 1, 5, 8, 0),                    &
      subshell_term_LS(3, 2, 3, 8, 0),                    &
      subshell_term_LS(3, 3, 3, 8, 0),                    &
      subshell_term_LS(3, 5, 1, 4, 0),                    &
      subshell_term_LS(3, 5, 1, 8, 0),                    &
      subshell_term_LS(3, 6, 1, 4, 0),                    &
      subshell_term_LS(3, 6, 1, 8, 0),                    &
      subshell_term_LS(3, 7, 1, 8, 0),                    &
      subshell_term_LS(3, 8, 1, 8, 0),                    &
      subshell_term_LS(3, 4, 3, 4, 0),                    &
      subshell_term_LS(3, 4, 3, 8, 0),                    &
      subshell_term_LS(3, 1, 3,10, 0),                    &
      subshell_term_LS(3, 2, 3,10, 0),                    &
      subshell_term_LS(3, 0, 1, 2, 0),                    &
      subshell_term_LS(3, 3, 1,10, 0),                    &
      subshell_term_LS(3, 4, 1,10, 0),                    &
      subshell_term_LS(3, 1, 7, 0, 0),                    &
      subshell_term_LS(3, 1, 5,12, 0),                    &
      subshell_term_LS(3, 2, 3, 0, 0),                    &
      subshell_term_LS(3, 2, 3,12, 0),                    &
      subshell_term_LS(3, 3, 3,12, 0),                    &
      subshell_term_LS(3, 3, 1, 0, 0),                    &
      subshell_term_LS(3, 4, 1,12, 0),                    &
      subshell_term_LS(3, 5, 1,12, 0),                    &
      subshell_term_LS(3, 4, 1, 0, 0),                    &
      subshell_term_LS(3, 6, 1,12, 0),                    &
      subshell_term_LS(3, 7, 1,12, 0),                    &
      subshell_term_LS(3, 1, 3,14, 0),                    &
      subshell_term_LS(3, 2, 1,14, 0),                    &
      subshell_term_LS(3, 3, 1,14, 0),                    &
      subshell_term_LS(3, 1, 3,16, 0),                    &
      subshell_term_LS(3, 2, 3,16, 0),                    &
      subshell_term_LS(3, 3, 1,16, 0),                    &
      subshell_term_LS(3, 4, 1,16, 0),                    &
      subshell_term_LS(3, 1, 1,18, 0),                    &
      subshell_term_LS(3, 2, 1,18, 0),                    &
      subshell_term_LS(3, 1, 3,20, 0),                    &
      subshell_term_LS(3, 2, 1,20, 0),                    &
      subshell_term_LS(3, 0, 1,24, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_g1 =(/                                      &
      subshell_term_LS(4, 0, 8, 8, 1)  /)
   type(subshell_term_LS), dimension(1:9), parameter ::   &
      term_LS_g2 =(/                                      &
      subshell_term_LS(4, 0, 9, 0, 0),                    &
      subshell_term_LS(4, 0, 7, 2, 2),                    &
      subshell_term_LS(4, 0, 7, 4, 0),                    &
      subshell_term_LS(4, 0, 7, 6, 2),                    &
      subshell_term_LS(4, 0, 7, 8, 0),                    &
      subshell_term_LS(4, 0, 7,10, 2),                    &
      subshell_term_LS(4, 0, 7,12, 0),                    &
      subshell_term_LS(4, 0, 7,14, 2),                    &
      subshell_term_LS(4, 0, 7,16, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_h1 =(/                                      &
      subshell_term_LS(5, 0,10,10, 1)  /)
   type(subshell_term_LS), dimension(1:11), parameter ::  &
      term_LS_h2 =(/                                      &
      subshell_term_LS(5, 0,11, 0, 0),                    & 
      subshell_term_LS(5, 0, 9, 2, 2),                    & 
      subshell_term_LS(5, 0, 9, 4, 0),                    & 
      subshell_term_LS(5, 0, 9, 6, 2),                    & 
      subshell_term_LS(5, 0, 9, 8, 0),                    & 
      subshell_term_LS(5, 0, 9,10, 2),                    & 
      subshell_term_LS(5, 0, 9,12, 0),                    & 
      subshell_term_LS(5, 0, 9,14, 2),                    & 
      subshell_term_LS(5, 0, 9,16, 0),                    & 
      subshell_term_LS(5, 0, 9,18, 2),                    & 
      subshell_term_LS(5, 0, 9,20, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_i1 =(/                                      &
      subshell_term_LS(6, 0,12,12, 1)  /)
   type(subshell_term_LS), dimension(1:13), parameter ::  &
      term_LS_i2 =(/                                      &
      subshell_term_LS(6, 0,13, 0, 0),                    & 
      subshell_term_LS(6, 0,11, 2, 2),                    & 
      subshell_term_LS(6, 0,11, 4, 0),                    & 
      subshell_term_LS(6, 0,11, 6, 2),                    & 
      subshell_term_LS(6, 0,11, 8, 0),                    & 
      subshell_term_LS(6, 0,11,10, 2),                    & 
      subshell_term_LS(6, 0,11,12, 0),                    & 
      subshell_term_LS(6, 0,11,14, 2),                    & 
      subshell_term_LS(6, 0,11,16, 0),                    & 
      subshell_term_LS(6, 0,11,18, 2),                    & 
      subshell_term_LS(6, 0,11,20, 0),                    & 
      subshell_term_LS(6, 0,11,22, 2),                    & 
      subshell_term_LS(6, 0,11,24, 0)  /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_k1 =(/                                      &
      subshell_term_LS(7, 0,14,14, 1)  /)
   type(subshell_term_LS), dimension(1:15), parameter ::  &
      term_LS_k2 =(/                                      &
      subshell_term_LS(7, 0,15, 0, 0),                    & 
      subshell_term_LS(7, 0,13, 2, 2),                    & 
      subshell_term_LS(7, 0,13, 4, 0),                    & 
      subshell_term_LS(7, 0,13, 6, 2),                    & 
      subshell_term_LS(7, 0,13, 8, 0),                    & 
      subshell_term_LS(7, 0,13,10, 2),                    & 
      subshell_term_LS(7, 0,13,12, 0),                    & 
      subshell_term_LS(7, 0,13,14, 2),                    & 
      subshell_term_LS(7, 0,13,16, 0),                    & 
      subshell_term_LS(7, 0,13,18, 2),                    & 
      subshell_term_LS(7, 0,13,20, 0),                    & 
      subshell_term_LS(7, 0,13,22, 2),                    & 
      subshell_term_LS(7, 0,13,24, 0),                    & 
      subshell_term_LS(7, 0,13,26, 2),                    & 
      subshell_term_LS(7, 0,13,28, 0) /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_l1 =(/                                      &
      subshell_term_LS(8, 0,16,16, 1)  /)
   type(subshell_term_LS), dimension(1:17), parameter ::  &
      term_LS_l2 =(/                                      &
      subshell_term_LS(8, 0,17, 0, 0),                    & 
      subshell_term_LS(8, 0,15, 2, 2),                    & 
      subshell_term_LS(8, 0,15, 4, 0),                    & 
      subshell_term_LS(8, 0,15, 6, 2),                    & 
      subshell_term_LS(8, 0,15, 8, 0),                    & 
      subshell_term_LS(8, 0,15,10, 2),                    & 
      subshell_term_LS(8, 0,15,12, 0),                    & 
      subshell_term_LS(8, 0,15,14, 2),                    & 
      subshell_term_LS(8, 0,15,16, 0),                    & 
      subshell_term_LS(8, 0,15,18, 2),                    & 
      subshell_term_LS(8, 0,15,20, 0),                    & 
      subshell_term_LS(8, 0,15,22, 2),                    & 
      subshell_term_LS(8, 0,15,24, 0),                    & 
      subshell_term_LS(8, 0,15,26, 2),                    & 
      subshell_term_LS(8, 0,15,28, 0),                    & 
      subshell_term_LS(8, 0,15,30, 2),                    & 
      subshell_term_LS(8, 0,15,32, 0) /)
   type(subshell_term_LS), dimension(1:1), parameter ::   &
      term_LS_m1 =(/                                      &
      subshell_term_LS(9, 0,18,18, 1)  /)
   type(subshell_term_LS), dimension(1:19), parameter ::  &
      term_LS_m2 =(/                                      &
      subshell_term_LS(9, 0,19, 0, 0),                    & 
      subshell_term_LS(9, 0,17, 2, 2),                    & 
      subshell_term_LS(9, 0,17, 4, 0),                    & 
      subshell_term_LS(9, 0,17, 6, 2),                    & 
      subshell_term_LS(9, 0,17, 8, 0),                    & 
      subshell_term_LS(9, 0,17,10, 2),                    & 
      subshell_term_LS(9, 0,17,12, 0),                    & 
      subshell_term_LS(9, 0,17,14, 2),                    & 
      subshell_term_LS(9, 0,17,16, 0),                    & 
      subshell_term_LS(9, 0,17,18, 2),                    & 
      subshell_term_LS(9, 0,17,20, 0),                    & 
      subshell_term_LS(9, 0,17,22, 2),                    & 
      subshell_term_LS(9, 0,17,24, 0),                    & 
      subshell_term_LS(9, 0,17,26, 2),                    & 
      subshell_term_LS(9, 0,17,28, 0),                    & 
      subshell_term_LS(9, 0,17,30, 2),                    & 
      subshell_term_LS(9, 0,17,32, 0),                    & 
      subshell_term_LS(9, 0,17,34, 2),                    & 
      subshell_term_LS(9, 0,17,36, 0) /)
      END MODULE jj2lsj_C 
!
!***********************************************************************
!                                                                      *
      MODULE m_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06  
      INTEGER, PARAMETER :: NNNW = 120 
      INTEGER, DIMENSION(NNNW) :: NQ1, NQ2 
      INTEGER, DIMENSION(NNNW) :: JJC1, JJC2
      INTEGER, DIMENSION(3,NNNW) :: JJQ1, JJQ2 
      INTEGER, DIMENSION(NNNW) :: JLIST, KLIST 
      INTEGER :: NPEEL, NCORE 
      END MODULE m_C 
!
!***********************************************************************
!                                                                      *
      MODULE npar_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06  
      INTEGER :: NPARM 
      REAL(DOUBLE), DIMENSION(2) :: PARM 
      END MODULE npar_C 

      MODULE nsmdat_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  06:36:34  12/28/06  
      REAL(DOUBLE) :: SQN,DMOMNM,QMOMB
      REAL(DOUBLE) :: HFSI, HFSD, HFSQ 
      REAL(DOUBLE) :: SMSI, SMSD, SMSQ 
      END MODULE nsmdat_C 
!
!***********************************************************************
!                                                                      *
      MODULE orb_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE, BYTE 
!...Created by Pacific-Sierra Research 77to90  4.3E  08:57:22  12/25/06  

!     Formerly orb10_C and orb10r_C
      INTEGER, PARAMETER :: NNNW = 120 
      CHARACTER(LEN=2), DIMENSION(NNNW) :: NH 
      CHARACTER(LEN=2), DIMENSION(NNNW) :: NHR 
!     Formerly orb1_C
      REAL(DOUBLE), DIMENSION(NNNW) :: E, GAMA 
!     Formerly orb2_C and orb2r_C
      INTEGER :: NCF, NW, NCFR, NWR 
      INTEGER(BYTE), DIMENSION(:,:), pointer :: IQA
      INTEGER, DIMENSION(:,:), pointer :: IQAR
!     Formerly orb4_C and orb4r_C
      INTEGER, DIMENSION(NNNW) :: NP, NAK, NPR, NAKR 
!     Formerly orb5_c
      INTEGER, DIMENSION(NNNW) :: NKL, NKJ 
      END MODULE orb_C 
!
!***********************************************************************
!                                                                      *
      MODULE prnt_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE, BYTE
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07  
      INTEGER :: NVEC, NVECMX 
      INTEGER :: NVECFF, NVECMXFF 
      INTEGER :: NVECII, NVECMXII 
!     REAL(DOUBLE) :: PNIVEC 
!     REAL(DOUBLE) :: PNIVECFF 
      REAL(DOUBLE) :: PNIVECII 
      INTEGER, DIMENSION(:), pointer :: ivec, ivecff, ivecii
      END MODULE prnt_C 
!
!***********************************************************************
!                                                                      *
      MODULE stat_C 
!
!***********************************************************************
!                                                                      *
!     Formerly stat_C and statr_C
      USE vast_kind_param, ONLY:  DOUBLE, BYTE
!...Created by Pacific-Sierra Research 77to90  4.3E  07:21:55   1/ 6/07  
!     REAL :: PNJQSR, PJCUPR !  From statr
!     REAL(DOUBLE) :: PNTJQS, PNJCUP 
      INTEGER, DIMENSION(:,:,:), pointer :: jqsar
      INTEGER(BYTE), DIMENSION(:,:,:), pointer :: jqsa
      INTEGER, DIMENSION(:,:), pointer :: jcupar
      INTEGER(BYTE), DIMENSION(:,:), pointer :: jcupa
      END MODULE stat_C 
!
!***********************************************************************
!                                                                      *
      MODULE syma_C 
!                                                                      *
!***********************************************************************
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07  
!     REAL(DOUBLE) :: PIATJP, PIASPA 
!     REAL(DOUBLE) :: PIATJPFF, PIASPAFF 
!     REAL(DOUBLE) :: PIATJPII, PIASPAII 
      INTEGER, DIMENSION(:), pointer :: iatjpo, iaspar
      INTEGER, DIMENSION(:), pointer :: iatjpoff, iasparff
      INTEGER, DIMENSION(:), pointer :: iatjpoii, iasparii
!     .. I cannot tell whether that above are "zero" or "oh"
      END MODULE syma_C 
!
!***********************************************************************
!                                                                      *
      MODULE terms_C 
!                                                                      *
!***********************************************************************
!...Created by Pacific-Sierra Research 77to90  4.3E  06:16:25   2/14/04  
      implicit none
      INTEGER, DIMENSION(31) :: ITAB 
      INTEGER, DIMENSION(32) :: JTAB 
      INTEGER, DIMENSION(327) :: NTAB 
      INTEGER :: NROWS
      INTEGER, PRIVATE :: i
      DATA NROWS/ 31/
!
!   A row is defined by a subshell angular momentum and an occupation
!   number
!
!   Each entry ITAB gives the number of terms in a row
!
!   Each entry JTAB gives the starting location -1 of the first triad
!   in a row
!
!   Each triad in NTAB is (v,w,2J+1); here v is the seniority,
!   w resolves any degeneracy in the seniority scheme, and J is the
!   subshell total angular momentum
!
!   Empty subshell or full subshell
!
      DATA (ITAB(I),I=1,1)/ 1/
      DATA (JTAB(I),I=1,1)/ 0/
      DATA (NTAB(I),I=1,3)/ 0, 0, 1/
!
!   s, p-   (j = 1/2)
!
      DATA (ITAB(I),I=2,2)/ 1/
      DATA (JTAB(I),I=2,2)/ 3/
      DATA (NTAB(I),I=4,6)/ 1, 0, 2/
!
!   p, d-   (j = 3/2)
!
      DATA (ITAB(I),I=3,4)/ 1, 2/
      DATA (JTAB(I),I=3,4)/ 6, 9/
      DATA (NTAB(I),I=7,15)/ 1, 0, 4, 0, 0, 1, 2, 0, 5/
!
!  d, f-   (j = 5/2)
!
      DATA (ITAB(I),I=5,7)/ 1, 3, 3/
      DATA (JTAB(I),I=5,7)/ 15, 18, 27/
      DATA (NTAB(I),I=16,36)/ 1, 0, 6, 0, 0, 1, 2, 0, 5, 2, 0, 9, 1, 0, 6, 3, 0&
         , 4, 3, 0, 10/
!
!   f, g-   (j = 7/2)
!
      DATA (ITAB(I),I=8,11)/ 1, 4, 6, 8/
      DATA (JTAB(I),I=8,11)/ 36, 39, 51, 69/
      DATA (NTAB(I),I=37,93)/ 1, 0, 8, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, 1, &
         0, 8, 3, 0, 4, 3, 0, 6, 3, 0, 10, 3, 0, 12, 3, 0, 16, 0, 0, 1, 2, 0, 5&
         , 2, 0, 9, 2, 0, 13, 4, 0, 5, 4, 0, 9, 4, 0, 11, 4, 0, 17/
!
!   g, h-   (j = 9/2)
!
      DATA (ITAB(I),I=12,16)/ 1, 5, 10, 18, 20/
      DATA (JTAB(I),I=12,16)/ 93, 96, 111, 141, 195/
      DATA (NTAB(I),I=94,255)/ 1, 0, 10, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, 2&
         , 0, 17, 1, 0, 10, 3, 0, 4, 3, 0, 6, 3, 0, 8, 3, 0, 10, 3, 0, 12, 3, 0&
         , 14, 3, 0, 16, 3, 0, 18, 3, 0, 22, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, &
         13, 2, 0, 17, 4, 0, 1, 4, 0, 5, 4, 0, 7, 4, 0, 9, 4, 1, 9, 4, 0, 11, 4&
         , 0, 13, 4, 1, 13, 4, 0, 15, 4, 0, 17, 4, 0, 19, 4, 0, 21, 4, 0, 25, 1&
         , 0, 10, 3, 0, 4, 3, 0, 6, 3, 0, 8, 3, 0, 10, 3, 0, 12, 3, 0, 14, 3, 0&
         , 16, 3, 0, 18, 3, 0, 22, 5, 0, 2, 5, 0, 6, 5, 0, 8, 5, 0, 10, 5, 0, &
         12, 5, 0, 14, 5, 0, 16, 5, 0, 18, 5, 0, 20, 5, 0, 26/
!
!   h, i-   (j = 11/2)
!
!   h, i-   (j = 11/2)
!
!   First two rows only
!
      DATA (ITAB(I),I=17,18)/ 1, 6/
      DATA (JTAB(I),I=17,19)/ 255, 258, 277/
      DATA (NTAB(I),I=256,276)/ 1, 0, 12, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, &
         2, 0, 17, 2, 0, 21/
!
!   i, k-   (j = 13/2)
!
!   First two rows only
!
      DATA (ITAB(I),I=23,24)/ 1, 7/
      DATA (JTAB(I),I=23,25)/ 276, 279, 301/
      DATA (NTAB(I),I=277,300)/ 1, 0, 14, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, &
         2, 0, 17, 2, 0, 21, 2, 0, 25/
!
!   k, l-   (j = 15/2)
!
!   First two rows only
!
      DATA (ITAB(I),I=30,31)/ 1, 8/
      DATA (JTAB(I),I=30,32)/ 300, 303, 328/
      DATA (NTAB(I),I=301,327)/ 1, 0, 16, 0, 0, 1, 2, 0, 5, 2, 0, 9, 2, 0, 13, &
         2, 0, 17, 2, 0, 21, 2, 0, 25, 2, 0, 29/
      END MODULE terms_C 
