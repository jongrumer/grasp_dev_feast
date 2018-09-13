      SUBROUTINE EDENSITYFIT(XVEC,YVEC,Z,PAR,NRNUC,P,F,RHO,RES)
!     Fits polynomial b1 + b2r^2 b3r^3 + b4r^4 to (r,rho) electron density data using least squares method
      IMPLICIT NONE
      DOUBLE PRECISION :: M(3,3), RM(3), MI(3,3), PM(3),MDET
      DOUBLE PRECISION :: X(600), Y(600), C(3,3), B(3), CI(3,3), P(4)
      DOUBLE  PRECISION :: XVEC(600), YVEC(600), RHO, Z
      DOUBLE PRECISION :: CDET, AU2FM, W(600), NORM, RES, F(4)
      DOUBLE PRECISION :: PI, CONST, A0, A1, A2
      DOUBLE PRECISION :: DR(4), FDSUM, PAR(2), FO90
      
      INTEGER :: I, NMIN, NMAX, NRNUC, NR

      PI = 3.14159265358979d0
      AU2FM = 52917.72083d0
      CONST = 27.2113834d0*AU2FM*1000.0d0

      ! 90% FALL OFF RADIUS FO90 DETERMINED - R <= FO90 DEFINES WITHIN NUCLEUS
      ! PAR(1) = 50% FALL OFF RADIUS C
      ! PAR(2) = SKIN THICKNESS A. 10% TO 90% FALL OFF DISTANCE T = 4.0*LN(3)*A
      FO90 = PAR(1)+2.0d0*LOG(3.0d0)*PAR(2)

      X(:) = XVEC(:)
      Y(:) = YVEC(:)

      ! DETERMINE FIRST POINT BEYOND FO90 IN GRID CALLED NR
      DO I=1,200
         IF(X(I+1).GT.FO90.AND.X(I).LT.FO90) THEN
            NR = I+1 
         END IF
      END DO

      ! Y(NR) AND X(NR) ARE RECALCULATED SO THAT X(NR) = FO90
      ! Y(NR) IS RECALCULATED TROUGH LINEAR INTERPOLATION BETWEEN Y(NR-1) AND Y(NR) 
      Y(NR) = Y(NR-1) + 
     :     (Y(NR)-Y(NR-1))/(X(NR)-X(NR-1))*(FO90-X(NR-1))
      X(NR) = FO90

      ! SINCE FIRST FEW DATA POINTS IN DENSITY ARE NOT RELIABLE DUE TO DIVISION WITH SMALL
      ! R^2 VALUES (STAGGERING IS SEEN), NMIN IS THE FIRST POINT TO BE USED IN THE SUBSEQUENT LEAST SQUARES FIT
      NMIN = 8

      NMAX = NR          ! SETS NMAX TO NR
      NRNUC = NR


! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
! TO DETERMINE RHO:
! BELOW PM(1) + PM(2)*X(I)**2 + PM(3)*X(I)**4 IS FITTED TROUGH DATA POINTS
! (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)), (X(NMIN+2),Y(NMIN+2)).
! THE DENSITY IN THE FIRST NMIN-1 POINTS IS THEN REPLACED BY EXTRAPOLATING 
! THE POLYNOMIAL ABOVE.
! ESPECIALLY WE HAVE: RHO(0) = PM(1)
! (THIS EXTAPOLATION IS USED SINCE DATA POINTS SMALLER THAN NMIN ARE NOT RELIABLE )

      RM(1) = Y(NMIN)
      RM(2) = Y(NMIN+1)
      RM(3) = Y(NMIN+2)

      M(1,1) = 1.0d0
      M(1,2) = X(NMIN)**2.0d0
      M(1,3) = X(NMIN)**4.0d0
      M(2,1) = 1.0d0
      M(2,2) = X(NMIN+1)**2.0d0
      M(2,3) = X(NMIN+1)**4.0d0
      M(3,1) = 1.d0
      M(3,2) = X(NMIN+2)**2.0d0
      M(3,3) = X(NMIN+2)**4.0d0

      MDET = M(1,1)*M(2,2)*M(3,3)-M(1,1)*M(2,3)*M(3,2)-
     :     M(1,2)*M(2,1)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+
     :     M(1,3)*M(2,1)*M(3,2)-M(1,3)*M(2,2)*M(3,1)

! Determine inverse of C matrix
      MI(1,1) = (M(2,2)*M(3,3)-M(2,3)*M(3,2))/MDET
      MI(1,2) = (M(1,3)*M(3,2)-M(1,2)*M(3,3))/MDET
      MI(1,3) = (M(1,2)*M(2,3)-M(1,3)*M(2,2))/MDET

      MI(2,1) = (M(2,3)*M(3,1)-M(2,1)*M(3,3))/MDET
      MI(2,2) = (M(1,1)*M(3,3)-M(1,3)*M(3,1))/MDET
      MI(2,3) = (M(1,3)*M(2,1)-M(1,1)*M(2,3))/MDET

      MI(3,1) = (M(2,1)*M(3,2)-M(2,2)*M(3,1))/MDET
      MI(3,2) = (M(1,2)*M(3,1)-M(1,1)*M(3,2))/MDET
      MI(3,3) = (M(1,1)*M(2,2)-M(1,2)*M(2,1))/MDET

      ! Determine parameters
      PM(1) = MI(1,1)*RM(1)+MI(1,2)*RM(2)+MI(1,3)*RM(3)
      PM(2) = MI(2,1)*RM(1)+MI(2,2)*RM(2)+MI(2,3)*RM(3)
      PM(3) = MI(3,1)*RM(1)+MI(3,2)*RM(2)+MI(3,3)*RM(3)

      ! Finallay RHO(0) in au^{-3} is determined
      RHO = PM(1)

      ! REPLACE THE DENSITY IN THE FIRST NMIN-1 POINTS WITH
      ! VALUES GIVEN BY THE EXTRAPOLATED FITTED POLYNOMIAL ABOVE
      DO I=1,NMIN-1
         YVEC(I) = PM(1) + PM(2)*X(I)**2.0d0 + PM(3)*X(I)**4.0d0
      END DO

      ! ------------------------------------------------------------------------------
      ! ------------------------------------------------------------------------------
      ! START LEAST SQUARE FIT PROCEDURE FOR DATA POINTS 
      ! (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)),...,(X(N),Y(N))
      ! THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
      ! P(2)*X(I)**2 + P(3)*X(I)**3 + P(4)*X(I)**4      
 
      NORM = -(Y(NMAX)-RHO)/AU2FM**3.0d0
      DO I=NMIN,NMAX
         X(I) = AU2FM*X(I)
         Y(I) = (Y(I)-RHO)/AU2FM**3.0d0
         Y(I) = Y(I)/NORM
         W(I) = X(I)                        ! WEIGHTS SET TO R(I) TO COMPENSATE FOR EXPONENTIAL GRID DENSITY
         !write(*,*) 'W(I): ', W(I)
      END DO
!     wrIte(*,*) 'rmax: ', X(N)

!     DetermINe B_l aNd C_{kl} matrIxelemeNts
      B(:) = 0.0d0
      C(:,:) = 0.0d0
      DO I=NMIN,NMAX
         B(1) = B(1) + Y(I)*X(I)**2.0d0*W(I)
         B(2) = B(2) + Y(I)*X(I)**3.0d0*W(I)
         B(3) = B(3) + Y(I)*X(I)**4.0d0*W(I)
         C(1,1) = C(1,1) + X(I)**4.0d0*W(I)
         C(2,2) = C(2,2) + X(I)**6.0d0*W(I)
         C(3,3) = C(3,3) + X(I)**8.0d0*W(I)
         C(1,2) = C(1,2) + X(I)**5.0d0*W(I)
         C(1,3) = C(1,3) + X(I)**6.0d0*W(I)
         C(2,3) = C(2,3) + X(I)**7.0d0*W(I)
      END DO
      C(2,1) = C(1,2)
      C(3,1) = C(1,3)
      C(3,2) = C(2,3)

!     DetermINe determINaNt
      CDET = C(1,1)*C(2,2)*C(3,3)-C(1,1)*C(2,3)*C(3,2)-
     :     C(1,2)*C(2,1)*C(3,3)+C(1,2)*C(2,3)*C(3,1)+
     :     C(1,3)*C(2,1)*C(3,2)-C(1,3)*C(2,2)*C(3,1)

!     DetermINe INverse of C matrIx
      CI(1,1) = (C(2,2)*C(3,3)-C(2,3)*C(3,2))/CDET
      CI(1,2) = (C(1,3)*C(3,2)-C(1,2)*C(3,3))/CDET
      CI(1,3) = (C(1,2)*C(2,3)-C(1,3)*C(2,2))/CDET

      CI(2,1) = (C(2,3)*C(3,1)-C(2,1)*C(3,3))/CDET
      CI(2,2) = (C(1,1)*C(3,3)-C(1,3)*C(3,1))/CDET
      CI(2,3) = (C(1,3)*C(2,1)-C(1,1)*C(2,3))/CDET

      CI(3,1) = (C(2,1)*C(3,2)-C(2,2)*C(3,1))/CDET
      CI(3,2) = (C(1,2)*C(3,1)-C(1,1)*C(3,2))/CDET
      CI(3,3) = (C(1,1)*C(2,2)-C(1,2)*C(2,1))/CDET

!     DetermINe fIttINg parameters
      P(2) = CI(1,1)*B(1)+CI(1,2)*B(2)+CI(1,3)*B(3)
      P(3) = CI(2,1)*B(1)+CI(2,2)*B(2)+CI(2,3)*B(3)
      P(4) = CI(3,1)*B(1)+CI(3,2)*B(2)+CI(3,3)*B(3)

      P(1) = RHO/AU2FM**3.0d0
      P(2) = P(2)*NORM
      P(3) = P(3)*NORM
      P(4) = P(4)*NORM
 
      F(1) = 2.0d0*PI/3.0d0*Z*P(1)*CONST
      F(2) = 2.0d0*PI/10.0d0*Z*P(2)*CONST
      F(3) = 2.0d0*PI/15.0d0*Z*P(3)*CONST
      F(4) = 2.0d0*PI/21.0d0*Z*P(4)*CONST

!      write(*,*) 'F(1): ', F(1)
!      write(*,*) 'F(2): ', F(2)
!      write(*,*) 'F(3): ', F(3)
!      write(*,*) 'F(4): ', F(4)

      write(*,*)
      write(*,*) '        r [fm]                    rho(r)'
      write(*,*) '---------------------------------------------------'

      RES = 0.0d0
      DO I=NMIN,NMAX
         Y(I) = Y(I)*AU2FM**3.0d0*NORM+RHO
      RES = RES + (Y(I)-RHO
     :        -AU2FM**3.0d0*(P(2)*X(I)**2.0d0+P(3)*X(I)**3.0d0
     :        +P(4)*X(I)**4.0d0))**2.0d0
!     RES = RES + abs(Y(I)
!    :        -RHO-AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0))
      write(*,*)X(I),Y(I),Y(I)-RHO-AU2FM**3.0d0*(P(2)*X(I)**2.0d0
     :     +P(3)*X(I)**3.0d0+P(4)*X(I)**4.0d0)
      END DO

      RES = sqrt(RES/(NMAX-NMIN+1))
      RES = RES/RHO*1000.0d0  ! In per mille of RHO

      ! RESULTING PARAMETERS RETURNED BY SUBROUTINE
      P(1) = RHO                        ! IN AU^-3
      P(2) = P(2)*AU2FM**3.0d0          ! IN AU^-3*FM^-2
      P(3) = P(3)*AU2FM**3.0d0          ! IN AU^-3*FM^-3
      P(4) = P(4)*AU2FM**3.0d0          ! IN AU^-3*FM^-4

      write(*,*) 
      write(*,*) 'RES: ', RES
      write(*,*) 'RHO: ', RHO

!     PURELY FOR TESTING PUPOSES
      DR(1) = 1.290943024d0
      DR(2) = 79.46276275d0
      DR(3) = 590.1056868d0
      DR(4) = 4375.679306d0

      write(*,*) 
      FDSUM = 0.0d0
      DO I=1,4 
         write(*,*) I, F(I)*DR(I) 
         FDSUM = FDSUM + F(I)*DR(I) 
      END DO
      write(*,*) FDSUM
      write(*,*) 

      RETURN
      END SUBROUTINE EDENSITYFIT
      
