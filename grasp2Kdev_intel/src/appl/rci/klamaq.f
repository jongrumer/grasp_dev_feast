************************************************************************
*                                                                      *
      SUBROUTINE KLAMAQ (N,KAPPA,Z,FZALFA)
*                                                                      *
*   The function  F (Z*\alpha)  is estimated here. We use the series   *
*   expansion given by  Eqs (1) and (2) and the table of Bethe loga-   *
*   rithms in  S Klarsfeld and A Maquet, Physics Letters  43B (1973)   *
*   201. The vacuum-polarization contribution in Eq (2) is omitted.    *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL FIRST
*
      DIMENSION BETHE(36)
*
      COMMON/DEF2/C
     :      /DEF9/CVAC,PI
*
*----------------------------------------------------------------------*
*
      DATA BETHE/ 2.9841285D 00,
     :            2.8117699D 00,-0.0300167D 00,
     :            2.7676636D 00,-0.0381902D 00,
     :              -0.0052321D 00,
     :            2.7498118D 00,-0.0419549D 00,
     :              -0.0067409D 00,-0.0017337D 00,
     :            2.7408237D 00,-0.0440347D 00,
     :              -0.0076008D 00,-0.0022022D 00,
     :                 -0.0007721D 00,
     :            2.7356642D 00,-0.0453122D 00,
     :              -0.0081472D 00,-0.0025022D 00,
     :                 -0.0009628D 00,-0.0004079D 00,
     :            2.7324291D 00,-0.0461552D 00,
     :              -0.0085192D 00,-0.0027091D 00,
     :                 -0.0010945D 00,-0.0004997D 00,
     :                    -0.0002409D 00,
     :            2.7302673D 00,-0.0467413D 00,
     :              -0.0087850D 00,-0.0028591D 00,
     :                 -0.0011904D 00,-0.0005665D 00,
     :                   -0.0002904D 00,-0.0001539D 00/
*
*----------------------------------------------------------------------*
*
      DATA FIRST/.TRUE./
*
      DATA C401/0.0D 00/,C402/0.0D 00/,
     :     OVLFAC/0.0D 00/
*
*   Set up the constants
*
      IF (FIRST) THEN
*
         C401 = 11.0D 00/24.0D 00
         C402 = 3.0D 00/8.0D 00
         OVLFAC = 4.0D 00/3.0D 00
*
         FIRST = .FALSE.
*
      ENDIF
*
*   Ensure that the principal quantum number is in range
*
      IF ((N .LT. 1) .OR. (N .GT. 8)) THEN
         WRITE (*,300)
         WRITE (*,301) N
         STOP
      ENDIF
*
*   Determine the azimuthal quantum number
*
      IF     (KAPPA .GT. 0) THEN
         L =  KAPPA
      ELSEIF (KAPPA .EQ. 0) THEN
         WRITE (*,300)
         WRITE (*,302)
         STOP
      ELSE
         L = -KAPPA-1
      ENDIF
*
*   Ensure that the azimuthal quantum number is in range
*
      IF (L .GT. N-1) THEN
         WRITE (*,300)
         WRITE (*,303) KAPPA,N
         STOP
      ENDIF
*
*   Find the appropriate entry in the table
*
      LOC = (N*N-N)/2+L+1
      BETHEL = BETHE(LOC)
*
*   Determine the quantity in square brackets in eq. (1) of
*   Klarsfeld and Maquet
*
      TERM = -BETHEL
*
      IF (KAPPA .GT. 0) THEN
         TERM = TERM-C402/DBLE (L*(L+L+1))
      ELSE
         TERM = TERM+C402/DBLE ((L+1)*(L+L+1))
         IF (KAPPA .EQ. -1) THEN
            ZALFA = Z/C
            FACTOR = LOG (ZALFA)
            FACTOR = -(FACTOR+FACTOR)
            TERM = TERM+FACTOR+C401
         ENDIF
      ENDIF
*
      FZALFA = OVLFAC*TERM
*
      RETURN
*
  300 FORMAT ('KLAMAQ:')
  301 FORMAT (' Principal quantum number, ',1I2,
     :        ', should be in the range  1--8.')
  302 FORMAT (' Kappa is  0 .')
  303 FORMAT (' Kappa, ',1I3,', is out of range for n, ',1I2,'.')
*
      END
