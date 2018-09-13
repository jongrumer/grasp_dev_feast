************************************************************************
*                                                                      *
      SUBROUTINE VAC2
*                                                                      *
*   This routine sets up the second-order vacuum polarization poten-   *
*   tial using  equations (1) and (4) of  L Wayne Fullerton and  G A   *
*   Rinker, Jr,  Phys Rev A  13 (1976) 1283-1287.  The  potential is   *
*   accumulated  in  array  TB(I), I = 1, ..., N  which is in COMMON   *
*   block /TATB/ .                                                     *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*               [RCI92]: FUNK.                                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF9/CVAC,PI
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NCDIST/ZDIST(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
*   Overall initialization
*
      EPSI = PRECIS*PRECIS
      TWOCV = CVAC+CVAC
*
*   Potential for a point nucleus: equation (1)
*   (this is also the asymptotoc form for a finite nucleus)
*
      FACTOR = -(2.0D 00*Z)/(3.0D 00*PI*CVAC)
*
      TB(1) = 0.0D 00
*
      I = 1
    1 I = I+1
*
      RI = R(I)
      X = TWOCV*RI
      TBI = (FACTOR/RI)*FUNK (X,1)
*
      IF (ABS (TBI) .GE. EPSI) THEN
         TB(I) = TBI
         IF (I .LT. N) GOTO 1
      ELSE
         DO 2 K = I,N
            TB(K) = 0.0D 00
    2    CONTINUE
      ENDIF
*
*   Potential for a finite nucleus: equation (4)
*
      IF (NPARM .EQ. 2) THEN
*
         FACTOR = -2.0D 00/(3.0D 00*CVAC**2)
*
*   Set up integrand
*
         TB(1) = 0.0D 00
*
         K = 1
    3    K = K+1
*
         RK = R(K)
         XK = TWOCV*RK
*
         TA(1) = 0.0D 00
         DO 4 I = 2,MTP
            XI = TWOCV*R(I)
            XM = ABS (XK-XI)
            XP = XK+XI
            TA(I) = ( FUNK (XM,0)-FUNK (XP,0) )*ZDIST(I)
    4    CONTINUE
*
         CALL QUAD (X)
*
         X = X*FACTOR/RK
*
*   Get out of loop if the asymptotic value has been attained
*
         IF (ABS (X) .GE. EPSI) THEN
            IF (ABS ((X-TB(K))/X) .GT. 1.0D-05) THEN
               TB(K) = X
               IF (K .LT. N) GOTO 3
            ENDIF
         ENDIF
*
      ENDIF
*
      RETURN
      END
