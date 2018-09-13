************************************************************************
*                                                                      *
      SUBROUTINE POSNFL (NUNIT,NREC)
*                                                                      *
*   This routine positions standard GRASP92 files at the end of NREC   *
*   records following the file header.                                 *
*   No subroutines called.                                             *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      REWIND NUNIT
*
      DO 1 I = 1,1+NREC
         READ (NUNIT)
    1 CONTINUE
*
      RETURN
      END
