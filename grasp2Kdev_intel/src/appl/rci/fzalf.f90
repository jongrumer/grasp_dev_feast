   !
   function fzalf(n,kappa,Z)
   ! function relci_qed_F(n,kappa,Z)                             result(F)
   !--------------------------------------------------------------------
   ! Estimates the function  F (Z*\alpha) by using an interpolation of 
   ! tabulated data from Mohr (1983) or a series expansion
   ! from S Klarsfeld and A Maquet (1973).
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: n, kappa
      real(kind=8), intent(in)         :: Z
      real(kind=8)                     :: F
      !
      if (n <= 2) then
         F = relci_qed_F_Mohr(n,kappa,z)
      else if (3 <= n   .and.   n <= 7) then
         F = relci_qed_F_Mohr_Kim(n,kappa,z)
      else
         F = zero
      end if
         if(iabs(kappa).ge.3) CALL KLAMAQ (N,KAPPA,Z,F)
!     F = relci_qed_F_Klarsfeld(n,kappa,z)
      fzalf=f
!      print *, n,kappa,z,f,'n,kappa,z,f'
      !
   end function fzalf      
   !
   !
   !
   !
   function relci_qed_F_Mohr(n,kappa,Z)                        result(F)
   !--------------------------------------------------------------------
   ! Computes the function  F (Z*\alpha) for the  1s  2s  2p-  2p  
   ! symmetries by interpolating in, or extrapolating from, the table 
   ! due to  P J Mohr. See  P J Mohr, At Data Nucl Data Tables 29 
   ! (1983) 453. 
   ! This procedure is adapted from RCI92 of GRASP92, written
   ! by Farid A Parpia, to the Fortran 95 standard.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: n, kappa
      real(kind=8), intent(in)         :: Z
      real(kind=8)                     :: F, value
      !
      ! Number of data points
      integer, parameter       :: numval = 12
      real(kind=8), parameter :: accy   = 1.0e-3
      !
      ! 1s data
      real(kind=8), dimension(numval), parameter :: val1s =     &
         (/   10.3168D+0,  4.6540D+0,   3.2460D+0,   2.5519D+0,  &
               2.1351D+0,  1.8644D+0,   1.6838D+0,   1.5675D+0,  &
               1.5032D+0,  1.4880D+0,   1.5317D+0,   1.6614D+0  /)
      ! 2s data
      real(kind=8), dimension(numval), parameter :: val2s =     &
         (/   10.5468D+0,  4.8930D+0,   3.5063D+0,   2.8391D+0,  &
               2.4550D+0,  2.2244D+0,   2.0948D+0,   2.0435D+0,  &
               2.0650D+0,  2.1690D+0,   2.3870D+0,   2.7980D+0  /)
      ! 2p- data
      real(kind=8), dimension(numval), parameter :: val2p1 =    &
         (/   -0.1264D+0, -0.1145D+0,  -0.0922D+0,  -0.0641D+0,  &
              -0.0308D+0,  0.0082D+0,   0.0549D+0,   0.1129D+0,  &
               0.1884D+0,  0.2934D+0,   0.4530D+0,   0.7250D+0  /)
      ! 2p data
      real(kind=8), dimension(numval), parameter :: val2p3 =    &
         (/    0.1235D+0,  0.1303D+0,   0.1436D+0,   0.1604D+0,  &
               0.1794D+0,  0.1999D+0,   0.2215D+0,   0.2440D+0,  &
               0.2671D+0,  0.2906D+0,   0.3141D+0,   0.3367D+0  /)
      ! Z data
      real(kind=8), dimension(numval), parameter :: arg =       &
         (/    1.0D+0,    10.0D+0,     20.0D+0,     30.0D+0,     &
              40.0D+0,    50.0D+0,     60.0D+0,     70.0D+0,     &
              80.0D+0,    90.0D+0,    100.0D+0,    110.0D+0     /)
      !
      ! Interpolate or issue error message as appropriate
      if (n == 1) then
         select case(kappa)
         case(-1)
            call interp(arg,val1s,numval,z,value,accy)
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr(): program stop A."
         end select
      else if (n == 2) then
         select case(kappa)
         case(-1)
            call interp(arg,val2s, numval,z,value,accy)
         case(1)
            call interp(arg,val2p1,numval,z,value,accy)
         case(-2)
            call interp(arg,val2p3,numval,z,value,accy)
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr(): program stop B."
         end select
      else
         print *, "Principal quantum number, ",n,"should be either 1 or 2."
         stop     "relci_qed_F_Mohr(): program stop B."
      end if
      !
      F = value
      relci_qed_F_Mohr=F
      !
   end function relci_qed_F_Mohr
   !
   !
   function relci_qed_F_Mohr_Kim(n,kappa,Z)                    result(F)
   !--------------------------------------------------------------------
   ! Computes the function  F (Z*\alpha) for the  1s  2s  2p-  2p  
   ! symmetries by interpolating in, or extrapolating from, the table 
   ! due to  P J Mohr and Y-K Kim. See  P J Mohr and Y-K Kim, 
   ! Phys. Rev. A45 (1992) 2723.
   !
   ! Since no values are given for Z = 1, these values are estimated by
   ! extrapolation.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: n, kappa
      real(kind=8), intent(in)         :: Z
      real(kind=8)                     :: F, value
      !
      ! Number of data points
      integer, parameter       :: numval = 12
      real(kind=8), parameter :: accy   = 1.0e-3
      !
      ! 3s data
      real(kind=8), dimension(numval), parameter :: val3s =     &
         (/   10.5000D+0,  4.9524D+0,   3.5633D+0,   2.8940D+0,  &
               2.5083D+0,  2.2757D+0,   2.1431D+0,   2.0874D+0,  &
               2.1018D+0,  2.1935D+0,   2.3897D+0,   2.7609D+0  /)
      ! 3p- data
      real(kind=8), dimension(numval), parameter :: val3p1 =    &
         (/   -0.1300D+0, -0.1021D+0,  -0.0760D+0,  -0.0430D+0,  &
               0.0041D+0,  0.0414D+0,   0.0956D+0,   0.1623D+0,  &
               0.2483D+0,  0.3660D+0,   0.5408D+0,   0.8322D+0  /)
      ! 3p data
      real(kind=8), dimension(numval), parameter :: val3p3 =    &
         (/    0.1250D+0,  0.1421D+0,   0.1572D+0,   0.1761D+0,  &
               0.1977D+0,  0.2214D+0,   0.2470D+0,   0.2745D+0,  &
               0.3038D+0,  0.3350D+0,   0.3679D+0,   0.4020D+0  /)
      ! 3d- data
      real(kind=8), dimension(numval), parameter :: val3d3 =    &
         (/   -0.0440D+0, -0.0428D+0,  -0.0420D+0,  -0.0410D+0,  &
              -0.0396D+0, -0.0378D+0,  -0.0353D+0,  -0.0321D+0,  &
              -0.0279D+0, -0.0225D+0,  -0.0154D+0,  -0.0062D+0  /)
      !
      ! 4s data
      real(kind=8), dimension(numval), parameter :: val4s =     &
         (/   10.5000D+0,  4.9749D+0,   3.5834D+0,   2.9110D+0,  &
               2.5215D+0,  2.2842D+0,   2.1455D+0,   2.0814D+0,  &
               2.0840D+0,  2.1582D+0,   2.3262D+0,   2.6484D+0  /)
      ! 4p- data
      real(kind=8), dimension(numval), parameter :: val4p1 =    &
         (/   -0.1200D+0, -0.0963D+0,  -0.0690D+0,  -0.0344D+0,  &
               0.0064D+0,  0.0538D+0,   0.1098D+0,   0.1780D+0,  &
               0.2649D+0,  0.3819D+0,   0.5525D+0,   0.8311D+0  /)
      ! 4p data
      real(kind=8), dimension(numval), parameter :: val4p3 =    &
         (/    0.1250D+0,  0.1477D+0,   0.1630D+0,   0.1827D+0,  &
               0.2052D+0,  0.2299D+0,   0.2568D+0,   0.2858D+0,  &
               0.3170D+0,  0.3507D+0,   0.3868D+0,   0.4247D+0  /)
      ! 4d- data
      real(kind=8), dimension(numval), parameter :: val4d3 =    &
         (/   -0.0410D+0, -0.0403D+0,  -0.0399D+0,  -0.0387D+0,  &
              -0.0371D+0, -0.0348D+0,  -0.0317D+0,  -0.0276D+0,  &
              -0.0222D+0, -0.0149D+0,  -0.0053D+0,   0.0074D+0  /)
      !
      ! 5s data
      real(kind=8), dimension(numval), parameter :: val5s =     &
         (/   10.5000D+0,  4.9858D+0,   3.5923D+0,   2.9173D+0,  &
               2.5246D+0,  2.2833D+0,   2.1395D+0,   2.0686D+0,  &
               2.0619D+0,  2.1225D+0,   2.2696D+0,   2.5566D+0  /)
      ! 5p- data
      real(kind=8), dimension(numval), parameter :: val5p1 =    &
         (/   -0.1200D+0, -0.0933D+0,  -0.0652D+0,  -0.0299D+0,  &
               0.0116D+0,  0.0597D+0,   0.1161D+0,   0.1843D+0,  &
               0.2703D+0,  0.3848D+0,   0.5497D+0,   0.8150D+0  /)
      ! 5p data
      real(kind=8), dimension(numval), parameter :: val5p3 =    &
         (/    0.1300D+0,  0.1502D+0,   0.1662D+0,   0.1861D+0,  &
               0.2089D+0,  0.2341D+0,   0.2614D+0,   0.2910D+0,  &
               0.3229D+0,  0.3574D+0,   0.3946D+0,   0.4338D+0  /)
      ! 5d- data
      real(kind=8), dimension(numval), parameter :: val5d3 =    &
         (/   -0.0405D+0, -0.0396D+0,  -0.0387D+0,  -0.0374D+0,  &
              -0.0356D+0, -0.0331D+0,  -0.0297D+0,  -0.0252D+0,  &
              -0.0190D+0, -0.0108D+0,   0.0001D+0,   0.0145D+0  /)
      ! Z data
      real(kind=8), dimension(numval), parameter :: arg =       &
         (/    1.0D+0,    10.0D+0,     20.0D+0,     30.0D+0,     &
              40.0D+0,    50.0D+0,     60.0D+0,     70.0D+0,     &
              80.0D+0,    90.0D+0,    100.0D+0,    110.0D+0     /)
      !
      ! Interpolate or issue error message as appropriate
      if (n == 3) then
         select case(kappa)
         case(-1)
            call interp(arg,val3s,numval,z,value,accy)
         case(1)
            call interp(arg,val3p1,numval,z,value,accy)
         case(-2)
            call interp(arg,val3p3,numval,z,value,accy)
         case(2)
            call interp(arg,val3d3,numval,z,value,accy)
         case(-3)
            call interp(arg,val3d3,numval,z,value,accy)
            value = 0.9D+0 * value
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop A."
         end select
      else if (n == 4) then
         select case(kappa)
         case(-1)
            call interp(arg,val4s,numval,z,value,accy)
         case(1)
            call interp(arg,val4p1,numval,z,value,accy)
         case(-2)
            call interp(arg,val4p3,numval,z,value,accy)
         case(2)
            call interp(arg,val4d3,numval,z,value,accy)
         case(-3)
            call interp(arg,val4d3,numval,z,value,accy)
            value = 0.9D+0 * value
            value = zero
         case(3)
            value = zero
         case(-4)
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop B."
         end select
      else if (5 <= n  .and.   n <= 7) then
         select case(kappa)
         case(-1)
            call interp(arg,val5s,numval,z,value,accy)
         case(1)
            call interp(arg,val5p1,numval,z,value,accy)
         case(-2)
            call interp(arg,val5p3,numval,z,value,accy)
         case(2)
            call interp(arg,val5d3,numval,z,value,accy)
         case(-3)
            call interp(arg,val4d3,numval,z,value,accy)
            value = 0.9D+0 * value
            value = zero
         case(3:6)
            value = zero
         case(-7:-4)
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop C."
         end select
      else
         print *, "Principal quantum number, ",n,"should be in the interval "//&
                  "3 <= n <= 7."
         stop     "relci_qed_F_Mohr(): program stop D."
      end if
      !
      F = value
      relci_qed_F_Mohr_Kim=F
      !
   end function relci_qed_F_Mohr_Kim
   !
   !
   function relci_qed_F_Klarsfeld(n,kappa,Z)                   result(F)
   !--------------------------------------------------------------------
   ! Estimates the function  F (Z*\alpha) by using a series expansion
   ! from S Klarsfeld and A Maquet, Physics Letters  43B (1973) 201,
   ! Eqs (1) and (2) and the table of Bethe logarithms. The 
   ! vacuum-polarization contribution in Eq (2) is omitted. 
   ! This procedure is adapted from RCI92 of GRASP92, written
   ! by Farid A Parpia, to the Fortran 95 standard.
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: n, kappa
      real(kind=8), intent(in)         :: Z
      real(kind=8)                     :: F
      !
      real(kind=8), dimension(36), parameter :: bethe = &
         (/ 2.9841285D+0,   2.8117699D+0,  -0.0300167D+0,   2.7676636D+0, &
           -0.0381902D+0,  -0.0052321D+0,   2.7498118D+0,  -0.0419549D+0, &
           -0.0067409D+0,  -0.0017337D+0,   2.7408237D+0,  -0.0440347D+0, &
           -0.0076008D+0,  -0.0022022D+0,  -0.0007721D+0,   2.7356642D+0, &
           -0.0453122D+0,  -0.0081472D+0,  -0.0025022D+0,  -0.0009628D+0, &
           -0.0004079D+0,   2.7324291D+0,  -0.0461552D+0,  -0.0085192D+0, &
           -0.0027091D+0,  -0.0010945D+0,  -0.0004997D+0,  -0.0002409D+0, &
            2.7302673D+0,  -0.0467413D+0,  -0.0087850D+0,  -0.0028591D+0, &
           -0.0011904D+0,  -0.0005665D+0,  -0.0002904D+0,  -0.0001539D+0 /) 
      !
      real(kind=8), parameter :: C401 = 11.0D+0/24.0D+0,C = 137.036D+0,     &
                               C402 = 3.0D+0/8.0D+0, ovlfac = 4.0D+0/3.0D+0
      !
      integer       :: l, loc,one=1
      real(kind=8) :: bethel, factor, term
      !
      ! Ensure that the principal quantum number is in range
      if (n < 1   .or.   n > 8) then
         print *, "Principal quantum number,",n,", should be in the range 1-8."
         stop     "relci_qed_F_Klarsfeld(): program stop A."
      end if
      !
!      l = angular_momentum_l(kappa)
      if(kappa > 0) then
         l=kappa
      else
         l=iabs(kappa+1)
      endif
      print *, kappa,l,'kappa,l'
      if (l > n-1) then
         print *, "Kappa = ",kappa," is out of range for n = ",n,"."
         stop     "relci_qed_F_Klarsfeld(): program stop B."
      end if
      !
      ! Find the appropriate entry in the table
      loc    = (n*n-n)/2+l+1
      bethel = bethe(loc)
      !
      ! Determine the quantity in square brackets in eq.(1) of
      ! Klarsfeld and Maquet
      term = -bethel
      !
      if (kappa > 0) then
         term = term - c402 / (l*(l+l+one))
      else
         term = term + c402 / ((l+one)*(l+l+one))
         if (kappa == -1) then
            factor = log (Z/c)
            factor = - (factor + factor)
            term   = term + factor + c401
         end if
      end if
      !
!     print *, ovlfac, term,'ovlfac, term'
      F = ovlfac * term
      !
      relci_qed_F_Klarsfeld = F
!     print *, relci_qed_F_Klarsfeld,F,'relci_qed_F_Klarsfeld,F'
   
   end function relci_qed_F_Klarsfeld
