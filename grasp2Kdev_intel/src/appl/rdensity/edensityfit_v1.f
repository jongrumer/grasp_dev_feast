      SUBROUTINE EDENSITYFIT(XVEC,YVEC,Z,PAR,NPARFIT,
     :     DRMS,P,F,RHO,RES,NRNUC)
!     Fits polynomial b1 + b2r^2 b3r^3 + b4r^4 to (r,rho) electron density data using least squares method
      IMPLICIT NONE
      DOUBLE PRECISION :: M(3,3), RM(3), MI(3,3), PM(3),MDET
      DOUBLE PRECISION :: X(1000), Y(1000), C(3,3), B(3), CI(3,3), P(5)
      DOUBLE  PRECISION :: XVEC(1000), YVEC(1000), RHO, Z, DRMS
      DOUBLE PRECISION :: CDET, AU2FM, W(1000), NORM, RES, F(5)
      DOUBLE PRECISION :: PI, CONST, A0, A1, A2
      DOUBLE PRECISION :: DR(4), FDSUM, PAR(2), FO90
      DOUBLE PRECISION :: DR2(4),PARF(2),PARF2(2),R2,DR4,DR6

      DOUBLE PRECISION :: DY(1000),DF2(3),DE
      
      INTEGER :: I, NMIN, NMAX, NRNUC, NR, NPARFIT

      PI = 3.14159265358979d0
      AU2FM = 52917.72083d0
      CONST = 27.2113834d0*AU2FM*1000.0d0

      DR2(1) = DRMS

      ! 90% FALL OFF RADIUS FO90 DETERMINED - R <= FO90 DEFINES WITHIN NUCLEUS
      ! PAR(1) = 50% FALL OFF RADIUS C
      ! PAR(2) = SKIN THICKNESS A. 10% TO 90% FALL OFF DISTANCE T = 4.0*LN(3)*A
      FO90 = PAR(1)+2.0d0*LOG(3.0d0)*PAR(2)
      PARF(1) = PAR(1)*AU2FM
      PARF(2) = PAR(2)*AU2FM

      X(:) = XVEC(:)
      Y(:) = YVEC(:)

      ! DETERMINE FIRST POINT BEYOND FO90 IN GRID CALLED NR
      DO I=1,1000
         IF(X(I+1).GT.FO90.AND.X(I).LT.FO90) THEN
            NR = I+1
            EXIT
         END IF
      END DO
      !write(*,*) X(NR)

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
         !YVEC(I) = PM(1) + PM(2)*X(I)**2.0d0 + PM(3)*X(I)**4.0d0
      END DO

      ! ------------------------------------------------------------------------------
      ! ------------------------------------------------------------------------------
      ! START LEAST SQUARE FIT PROCEDURE FOR DATA POINTS 
      ! (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)),...,(X(N),Y(N))
 
      NORM = -(Y(NMAX)-RHO)/AU2FM**3.0d0

      DO I=NMIN,NMAX
         X(I) = AU2FM*X(I)
         Y(I) = (Y(I)-RHO)/AU2FM**3.0d0
         Y(I) = Y(I)/NORM
         W(I) = X(I)                        ! WEIGHTS SET TO R(I) TO COMPENSATE FOR EXPONENTIAL GRID DENSITY
         !write(*,*) 'W(I): ', W(I)
         !write(*,*) X(I), Y(I)
      END DO
!     wrIte(*,*) 'rmax: ', X(N)


      IF(NPARFIT.EQ.4) THEN
      ! THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
      ! P(2)*X(I)**2 + P(3)*X(I)**3 + P(4)*X(I)**4      
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
     :        C(1,2)*C(2,1)*C(3,3)+C(1,2)*C(2,3)*C(3,1)+
     :        C(1,3)*C(2,1)*C(3,2)-C(1,3)*C(2,2)*C(3,1)
         
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
         
!     write(*,*) 'F(1): ', F(1)
!     write(*,*) 'F(2): ', F(2)
!     write(*,*) 'F(3): ', F(3)
!     write(*,*) 'F(4): ', F(4)
         
         RES = 0.0d0
         DO I=NMIN,NMAX
            Y(I) = Y(I)*AU2FM**3.0d0*NORM+RHO
            RES = RES + (Y(I)-RHO
     :           -AU2FM**3.0d0*(P(2)*X(I)**2.0d0+P(3)*X(I)**3.0d0
     :           +P(4)*X(I)**4.0d0))**2.0d0
!            write(*,*)X(I),Y(I),Y(I)-RHO-AU2FM**3.0d0*(P(2)*X(I)**2.0d0
!     :           +P(3)*X(I)**3.0d0+P(4)*X(I)**4.0d0)
         END DO
         
         RES = sqrt(RES/(NMAX-NMIN+1))
         RES = RES/RHO*1000.0d0 ! In per mille of RHO
         
! RESULTING PARAMETERS RETURNED BY SUBROUTINE
         P(1) = RHO             ! IN AU^-3
         P(2) = P(2)*AU2FM**3.0d0 ! IN AU^-3*FM^-2
         P(3) = P(3)*AU2FM**3.0d0 ! IN AU^-3*FM^-3
         P(4) = P(4)*AU2FM**3.0d0 ! IN AU^-3*FM^-4
         
         write(*,*) 
         write(*,'(a6,f15.5,a10)') 'RES: ', RES, 'PER MILLE'
         write(*,'(a6,f15.2,a6)') 'RHO: ', RHO, 'AU^-3'
         
!     PURELY FOR TESTING PUPOSES
         R2 = 3.d0/5.d0*PARF(1)**2.d0 
     :        + 7.d0/5.d0*PI**2.d0*PARF(2)**2.d0 + DR2(1)
         PARF2(1) =  SQRT(5.d0/3.d0*(R2 
     :        - 7.d0/5.d0*(PI*PARF(2))**2.d0))
         DR2(2) = 3.d0/7.d0*(PARF2(1)**4.d0 - PARF(1)**4.d0)
     :        + 18.d0/7.d0*(PI*PARF(2))**2.d0*(PARF2(1)**2.d0 - 
     :        PARF(1)**2.d0)
         DR2(3) = 3.d0/8.d0*(PARF2(1)**5.d0 - PARF(1)**5.d0)
     :        + 25.d0/8.d0*(PI*PARF(2))**2.d0*(PARF2(1)**3.d0 - 
     :        PARF(1)**3.d0)
     :        + 73.d0/8.d0*(PI*PARF(2))**4.d0*(PARF2(1) - PARF(1))
     :        + 51.d0/8.d0*(PI*PARF(2))**6.d0*
     :        (1.d0/PARF2(1) - 1.d0/PARF(1))
     :        - 16.d0/5.d0*(PI*PARF(2))**8.d0*(1.d0/PARF2(1)**3.d0 - 
     :        1.d0/PARF(1)**3.d0) 
     :        + 16.d0/5.d0*(PI*PARF(2))**10.d0*(1.d0/PARF2(1)**5.d0 - 
     :        1.d0/PARF(1)**5.d0) 
     :        - 16.d0/5.d0*(PI*PARF(2))**12.d0*(1.d0/PARF2(1)**7.d0 - 
     :        1.d0/PARF(1)**7.d0) 
         DR2(4) = 3.d0/9.d0*(PARF2(1)**6.d0 - PARF(1)**6.d0)
     :        + 11.d0/3.d0*(PI*PARF(2))**2.d0*(PARF2(1)**4.d0 - 
     :        PARF(1)**4.d0) + 239.d0/15.d0*(PI*PARF(2))**4.d0*
     :        (PARF2(1)**2.d0 - PARF(1)**2.d0)

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
         write(*,*) 
         FDSUM = 0.0d0
         DO I=1,4 
            write(*,*) I, F(I)*DR2(I) 
            FDSUM = FDSUM + F(I)*DR2(I) 
         END DO
         F(5) = FDSUM
         write(*,*) FDSUM
         write(*,*) 
         
      ELSE IF (NPARFIT.EQ.3) THEN
      ! THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
      ! P(2)*X(I)**2 + P(3)*X(I)**4      

         !     DetermINe B_l aNd C_{kl} matrIxelemeNts
         B(:) = 0.0
         C(:,:) = 0.0
         DO I=NMIN,NMAX
            B(1) = B(1) + Y(I)*X(I)**2.0*W(I)
            B(2) = B(2) + Y(I)*X(I)**4.0*W(I)
            C(1,1) = C(1,1) + X(I)**4.0*W(I)
            C(2,2) = C(2,2) + X(I)**8.0*W(I)
            C(1,2) = C(1,2) + X(I)**6.0*W(I)
         END DO
         C(2,1) = C(1,2)
         
!     DetermINe determINaNt
         CDET = C(1,1)*C(2,2)-C(1,2)*C(2,1)
         
!     DetermINe INverse of C matrIx
         CI(1,1) = C(2,2)/CDET
         CI(1,2) = -C(1,2)/CDET
         CI(2,1) = -C(2,1)/CDET
         CI(2,2) = C(1,1)/CDET
         
!     DetermINe fIttINg parameters
         P(2) = CI(1,1)*B(1)+CI(1,2)*B(2)
         P(3) = CI(2,1)*B(1)+CI(2,2)*B(2)
         
         P(1) = RHO/AU2FM**3.0
         P(2) = P(2)*NORM
         P(3) = P(3)*NORM
         
         F(1) = 2.0*PI/3.0*Z*P(1)*CONST
         F(2) = PI/5.0*Z*P(2)*CONST
         F(3) = 2.0*PI/21.0*Z*P(3)*CONST

         RES = 0.0
         C(:,:) = 0.0
         DO I=NMIN,NMAX
            Y(I) = Y(I)*AU2FM**3.0*NORM+RHO
            RES = RES + (Y(I)-RHO
     :           -AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0))**2.0
!     RES = RES + abs(Y(I)
!     :        -RHO-AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0))
      !write(*,*)X(I),Y(I),RHO+AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0)
       DY(I) = RHO+AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0);
       DY(I) = ABS(Y(I) - DY(I))
       DY(I) = DY(I)/(NORM*AU2FM**3.0)
       DY(I) = 1/DY(I)**2.0
       !DY(I) = 1.0/(ABS(Y(I) - DY(I)))**2.0
       !DY(I) = 5.0
       !write(*,*)X(I),Y(I),DY(I)
            C(1,1) = C(1,1) + X(I)**4.0*DY(I)
            C(2,2) = C(2,2) + X(I)**8.0*DY(I)
            C(1,2) = C(1,2) + X(I)**6.0*DY(I)
         END DO
         C(2,1) = C(1,2)

!     DetermINe determINaNt
         CDET = C(1,1)*C(2,2)-C(1,2)*C(2,1)
         
!     DetermINe INverse of C matrIx
         CI(1,1) = C(2,2)/CDET
         CI(1,2) = -C(1,2)/CDET
         CI(2,1) = -C(2,1)/CDET
         CI(2,2) = C(1,1)/CDET

         write(*,*)
         write(*,*) 'F(2): ', F(2)
         write(*,*) 'Error in F(2): ',sqrt(CI(1,1))*PI/5.0*Z*CONST*
     :        NORM
         DF2(1) = CI(1,1)*(PI/5.0*Z*CONST*NORM)**2.0
         !write(*,*) sqrt(CI(2,2))
         write(*,*) 'F(3): ', F(3)
         write(*,*) 'Error in F(3): ',sqrt(CI(2,2))*2.0*PI/21.0*Z*CONST*
     :        NORM
         DF2(2) = CI(2,2)*(2.0*PI/21.0*Z*CONST*NORM)**2.0
         write(*,*) CI(1,2)
         write(*,*) CI(2,1)
         DF2(3) = CI(1,2)*(2.0*PI/21.0*Z*CONST*NORM)*
     :        (PI/5.0*Z*CONST*NORM)
         !write(*,*) sqrt(CI(2,2))*2.0*PI/21.0*Z*P(3)*CONST*NORM


         RES = sqrt(RES/(NR-NMIN+1))
         RES = RES/RHO*1000.0   ! In per mille of RHO(0)
!     RES = RES/(N-NMIN+1)
!     wrIte(*,*) 'RESIdual: ', RES
         
! RESULTING PARAMETERS RETURNED BY SUBROUTINE
         P(1) = RHO             ! IN AU^-3
         P(2) = P(2)*AU2FM**3.0d0 ! IN AU^-3*FM^-2
         P(3) = P(3)*AU2FM**3.0d0 ! IN AU^-3*FM^-4
         
         write(*,*) 
         write(*,'(a6,f15.5,a10)') 'RES: ', RES, 'PER MILLE'
         write(*,'(a6,f15.2,a6)') 'RHO: ', RHO, 'AU^-3'
         
!     PURELY FOR TESTING PUPOSES
         R2 = 3.d0/5.d0*PARF(1)**2.d0 
     :        + 7.d0/5.d0*PI**2.d0*PARF(2)**2.d0 + DR2(1)
         PARF2(1) =  SQRT(5.d0/3.d0*(R2 
     :        - 7.d0/5.d0*(PI*PARF(2))**2.d0))
         DR2(2) = 3.d0/7.d0*(PARF2(1)**4.d0 - PARF(1)**4.d0)
     :        + 18.d0/7.d0*(PI*PARF(2))**2.d0*(PARF2(1)**2.d0 - 
     :        PARF(1)**2.d0)
         DR2(3) = 3.d0/9.d0*(PARF2(1)**6.d0 - PARF(1)**6.d0)
     :        + 11.d0/3.d0*(PI*PARF(2))**2.d0*(PARF2(1)**4.d0 - 
     :        PARF(1)**4.d0) + 239.d0/15.d0*(PI*PARF(2))**4.d0*
     :        (PARF2(1)**2.d0 - PARF(1)**2.d0)

         DR(1) = 1.290943024d0
         DR(2) = 79.46276275d0
         DR(3) = 4375.679306d0
         write(*,*) 
         FDSUM = 0.0d0
         DO I=1,3 
            write(*,*) I, F(I)*DR(I) 
            FDSUM = FDSUM + F(I)*DR(I) 
         END DO
         write(*,*) FDSUM
         write(*,*) 

         FDSUM = 0.0d0
         DO I=1,3 
            write(*,*) I, F(I)*DR2(I) 
            FDSUM = FDSUM + F(I)*DR2(I) 
         END DO
         DE = DR(2)**2.*DF2(1) + DR(3)**2.*DF2(2)
         !DE = DE + 2.*DR(2)*DR(3)*DF2(3)
         DE = sqrt(DE)
         F(5) = FDSUM
         write(*,*) FDSUM, DE
         write(*,*) 


      END IF
      write(*,*) 'RHO(c)(RHO(0)', Y(NR)/RHO
      
      RETURN

      END SUBROUTINE EDENSITYFIT
      
