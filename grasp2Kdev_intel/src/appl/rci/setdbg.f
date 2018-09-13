************************************************************************
*                                                                      *
      SUBROUTINE SETDBG
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*                                                                      *
*   Written by Farid A Parpia               Last update: 21 Oct 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL LDBPA,LDBPG,LDBPR

      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEBUGR/LDBPR(30)

      DO 1 I = 1,5
         LDBPA(I) = .FALSE.
         LDBPG(I) = .FALSE.
    1 CONTINUE

      DO 3 I = 1,30
         LDBPR(I) = .FALSE.
    3 CONTINUE

      RETURN
      END
