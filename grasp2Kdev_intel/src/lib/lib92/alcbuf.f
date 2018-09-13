*  CFF -- testing the replacement of iunpck in jqs
*
************************************************************************
*                                                                      *
      SUBROUTINE ALCBUF (MODE)
*                                                                      *
*   The arrays in COMMON/BUFFER/ are allocated (MODE = 1), realloca-   *
*   ted (MODE = 2), and deallocated (MODE = 3) in this routine.        *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*                                                                      *
*   Written by Farid A. Parpia.             Last update: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PLABEL,PCOEFF
      POINTER (PLABEL,LABELDUMMY)
      POINTER (PCOEFF,COEFFDUMMY)
*
      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
*
      IF (MODE .EQ. 1) THEN
*
*   Allocation
*
c        NBDIM = 1
         NBDIM = 10
         CALL ALLOC (PLABEL,6*NBDIM,4)
         CALL ALLOC (PCOEFF,  NBDIM,8)
*
      ELSEIF (MODE .EQ. 2) THEN
*
*   Reallocation to double storage
*
         NEWSIZ = NBDIM+NBDIM
         CALL RALC2D (PLABEL,6,NBDIM,6,NEWSIZ,4)
         CALL RALLOC (PCOEFF,NBDIM,NEWSIZ,8)
         NBDIM = NEWSIZ
*
      ELSEIF (MODE .EQ. 3) THEN
*
*   Deallocation
*
         CALL DALLOC (PLABEL)
         CALL DALLOC (PCOEFF)
*
      ELSE
*
*   Argument error
*
         PRINT *, 'ALCBUF: Invalid argument MODE = ',MODE
         STOP
*
      ENDIF
*
      RETURN
*
      END
