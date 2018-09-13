************************************************************************
*                                                                      *
      FUNCTION FZALF (N,KAPPA,Z)
*                                                                      *
*   An estimate of the function  F (Z*\alpha) is computed here.        *
*                                                                      *
*   Call(s) to: [RCI92]: KLAMAQ, MOHR.                                 *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1990   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      IF (N .LE. 2) THEN
         CALL MOHR (N,KAPPA,Z,VALUE)
      ELSE
         IF ( (KAPPA .EQ. -1) .OR.
     :        (KAPPA .EQ.  1) .OR.
     :        (KAPPA .EQ. -2) ) THEN
            NEFF = 2
            CALL MOHR (NEFF,KAPPA,Z,VALUE)
         ELSE
            CALL KLAMAQ (N,KAPPA,Z,VALUE)
         ENDIF
      ENDIF
*
      FZALF = VALUE
*
      RETURN
      END
