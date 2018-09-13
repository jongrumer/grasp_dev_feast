************************************************************************
*                                                                      *
      SUBROUTINE MRGCSL(NAME)
*                                                                      *
*   Entry routine for merging two csl lists                            *
*                                                                      *
************************************************************************
*
      implicit double precision (a-h,o-z)
      CHARACTER*24 NAME(2)

      COMMON/DEBUGA/LDBPA(5)
     :      /DEF1/EMN,IONCTY,NELEC,Z

      PRINT *
      PRINT *, 'MRGCSL: Execution begins ...'
*
*   Load the first  .csl  file
*
      CALL LDCSL1 (NCORER,NAME(1))
*
*   Load the second  .csl  file
*
      CALL LDCSL2 (NCORE,NAME(2))
*
*   Merge the two  .csl  lists, observe that there may be doublets 
*   among the CSF's
*
      CALL MERG12 (NAME,NCORER,NCORE)

      RETURN
      END
