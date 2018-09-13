      SUBROUTINE screenpar (ncore)
      IMPLICIT REAL*8           (a-h,o-z)
************************************************************************
*
* Purpose:
*   Compute hydrogenic screen parameters
*
* Input:
*   ncore
*   in the common - nw, nkj()
*
* Output:
*   in the common - sigma()
*
************************************************************************
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      POINTER (PNTRIQ,IQAdum)
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     &      /hydpar/sigma(NNNW)

       !...Core orbitals
      nelectron = 0
      DO i = 1, ncore
         sigma(i) = nelectron + (nkj(i) + 1) / 2
         nelectron = nelectron + nkj(i) + 1
      ENDDO

       !...Peel orbitals
      DO i = ncore + 1, nw
         sigma(i) = nelectron
      ENDDO

      RETURN
      END
