************************************************************************
*                                                                      *
      FUNCTION DSUBRS (EOL,I,J,jblock)
*                                                                      *
*   The coefficients d   for  I = r, J = s  are calculated here.       *
*                     rs                                               *
*                                                                      *
*                         NCMIN                                        *
*                          Sum     w     c          c                  *
*                         t = 1     t     r Gamma    s Gamma           *
*                                                t          t          *
*    (E)OL:        d   = -------------------------------------         *
*                   rs             NCMIN                               *
*                                   Sum     w                          *
*                                  t = 1     t                         *
*                                                                      *
*                                                                      *
*                                           /  NCF                     *
*    (E)AL:           d   =  delta     w   /   Sum   w                 *
*                      rs         rs    r /   t = 1   t                *
*                                                                      *
*                                                                      *
*   Written by Farid A Parpia               Last update: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
*   This routine assumes (and does not check) that I and J belong
*   the same block and the job is done within the block specified
*   by parameter jblock
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      POINTER (PWEIGH,WEIGHDUMMY)
      LOGICAL EOL
*
      POINTER (PNTRWT,WT(*))
      POINTER (PNEVEC,EVEC(*))
*
      COMMON/DEF5/PNTRWT,PWEIGH
     :      /EIGVEC/PNEVEC
*
      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      POINTER (pnevblk, nevblk(*))
      POINTER (pncmaxblk, ncmaxblk(*))
      COMMON/hblock2/pnevblk, pncmaxblk

      POINTER (pncfpast, ncfpast(*))
      POINTER (pncminpast, ncminpast(*))
      POINTER (pnevecpast, nevecpast(*))
      COMMON/pos/pncfpast,pncminpast,pnevecpast,ncftot,nvecsiz

      IF (EOL) THEN
         DSUBRS = 0.D0
         DO K = 1, nevblk(jblock)
            kcmin = k + ncminpast(jblock)
            ik = I + nevecpast(jblock) + (K-1)*NCFblk(jblock)
            jk = J + nevecpast(jblock) + (K-1)*NCFblk(jblock)
            DSUBRS = DSUBRS + EVEC(ik) * EVEC(jk) * WT(kcmin)
         ENDDO
      ELSEIF (I .EQ. J) THEN
         DSUBRS = WT(I)
      ELSE
         DSUBRS = 0.D0
      ENDIF

      RETURN
      END
