!***********************************************************************
!                                                                      *
      SUBROUTINE LODISO 
!                                                                      *
!   Loads the data from the  .iso  file.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 29 Sep 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:30:59   2/14/04  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE DEF_C 
      USE NPAR_C 
      USE NSMDAT_C, ONLY: SQN, DMOMNM, QMOMB 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: A, APARM, CPARM, EMNAMU 


!
!   Read and echo pertinent information from  .iso  file
!
!   Atomic number
!
      READ (22, *) Z 
!
!   Nuclear geometry
!
      READ (22, *) 
      READ (22, *) A 
      READ (22, *) 
      READ (22, *) APARM 
      READ (22, *) 
      READ (22, *) CPARM 
!
      IF (A /= 0.D0) THEN 
         NPARM = 2 
         PARM(1) = CPARM*FMTOAU 
         PARM(2) = APARM*FMTOAU 
      ELSE 
         NPARM = 0 
      ENDIF 
!
!   Nuclear mass
!
      READ (22, *) 
      READ (22, *) EMNAMU 
!
      IF (EMNAMU /= 0.D0) THEN 
         EMN = EMNAMU/AUMAMU 
      ELSE 
         EMN = 0.D0 
      ENDIF 
!
!   Nuclear spin and moments
!
      READ (22, *) 
      READ (22, *) SQN 
      READ (22, *) 
      READ (22, *) DMOMNM 
      READ (22, *) 
      READ (22, *) QMOMB 
 
      RETURN  
      END SUBROUTINE LODISO 
