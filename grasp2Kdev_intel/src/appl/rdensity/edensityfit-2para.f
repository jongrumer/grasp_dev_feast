      SUBROUTINE EDENSITYFIT(XVEC,YVEC,Z,NRNUC,P,F,RHO,RES)
!     Fits polynomial b1 + b2r^2 b3r^4 to (r,rho) electron density data using least squares method
      IMPLICIT NONE
      DOUBLE PRECISION :: M(3,3), RM(3), MI(3,3), PM(3),MDET
      DOUBLE PRECISION :: X(600), Y(600), C(2,2), B(2), CI(2,2), P(3)
      DOUBLE  PRECISION :: XVEC(600), YVEC(600), RHO, Z
      DOUBLE PRECISION :: CDET, AU2FM, W(600), NORM, RES, F(3)
      DOUBLE PRECISION :: PI, CONST, A0, A1, A2
      INTEGER :: I, N, NMIN, NRNUC
      
      PI = 3.14159265358979
      AU2FM = 52917.72083
      CONST = 27.2113834*AU2FM*1000.0
      
      NMIN = 8
      N = NRNUC
      
      X(:) = XVEC(:)
      Y(:) = YVEC(:)

! BELOW PM(1) + PM(2)*X(I)**2 + PM(3)*X(I)**4 IS FITTED TO DATA POINTS
! (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)), (X(NMIN+2),Y(NMIN+2)).
! THE DENSITY IN THE FIRST NMIN-1 POINTS IS THEN REPLACED BY EXTRAPOLATING 
! THE POLYNOMIAL ABOVE.
! ESPECIALLY WE HAVE: RHO(0) = PM(1)
! (THIS EXTAPOLATION IS USED SINCE DATA POINTS SMALLER THAN NMIN ARE NOT RELIABLE )

      RM(1) = Y(NMIN)
      RM(2) = Y(NMIN+1)
      RM(3) = Y(NMIN+2)

      M(1,1) = 1d0
      M(1,2) = X(NMIN)**2
      M(1,3) = X(NMIN)**4
      M(2,1) = 1d0
      M(2,2) = X(NMIN+1)**2
      M(2,3) = X(NMIN+1)**4
      M(3,1) = 1d0
      M(3,2) = X(NMIN+2)**2
      M(3,3) = X(NMIN+2)**4

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

      RHO = PM(1)
!      write(*,*) 'RHO: ', RHO

! REPLACE THE DENSITY IN THE FIRST NMIN-1 POINTS WITH
! VALUES GIVEN BY THE FITTED POLYNOMIAL ABOVE
      DO I=1,NMIN-1
         YVEC(I) = PM(1) + PM(2)*X(I)**2 + PM(3)*X(I)**4
      END DO

! START LEAST SQUARE FIT PROCEDURE FOR DATA POINTS 
! (X(NMIN),Y(NMIN)), (X(NMIN+1),Y(NMIN+1)),...,(X(N),Y(N))
! THE DATA POINTS, SUBTRACTED SO THAT Y(I) = Y(I)-RHO, ARE FITTED TO POLYNOMIAL:
! P(2)*X(I)**2 + P(3)*X(I)**4      
      NORM = -(Y(N)-RHO)/AU2FM**3.0
      DO I=NMIN,N
         X(I) = AU2FM*X(I)
         Y(I) = (Y(I)-RHO)/AU2FM**3.0
         Y(I) = Y(I)/NORM
         W(I) = 1.0
      END DO
!     wrIte(*,*) 'rmax: ', X(N)

!     DetermINe B_l aNd C_{kl} matrIxelemeNts
      B(:) = 0.0
      C(:,:) = 0.0
      DO I=NMIN,N
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
      DO I=NMIN,N
         Y(I) = Y(I)*AU2FM**3.0*NORM+RHO
      RES = RES + (Y(I)-RHO
     :        -AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0))**2.0
!     RES = RES + abs(Y(I)
!    :        -RHO-AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0))
!      write(*,*)X(I),Y(I),RHO+AU2FM**3.0*(P(2)*X(I)**2.0+P(3)*X(I)**4.0)
      END DO
      RES = sqrt(RES/(N-NMIN+1))
      RES = RES/RHO*1000.0  ! In per mille of RHO(0)
!     RES = RES/(N-NMIN+1)
!     wrIte(*,*) 'RESIdual: ', RES

      RETURN
      END SUBROUTINE EDENSITYFIT
      
