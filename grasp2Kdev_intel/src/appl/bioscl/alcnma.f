************************************************************************
*                                                                      *
      SUBROUTINE ALCNMA (PIPTR,PISLDR,PISLDR1,PXSLDR,NMDIM,IMODE)
*                                                                      *
*   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
*   and deallocates (IMODE = 3) storage for  certain arrays that are   *
*   local to SUBROUTINE TRSORT.                                        *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PIPTR,PISLDR,PXSLDR
      POINTER (PIPTR,IPTRDUMMY)  
      POINTER (PISLDR,ISLDRDUMMY)
      POINTER (PISLDR1,ISLDR1DUMMY)
      POINTER (PXSLDR,XSLDRDUMMY)
*
      IF (IMODE .EQ. 1) THEN
*
*   Initial array dimension
*
         NMDIM = 64
*
*   Allocate storage for arrays
*
	 CALL ALLOC (PIPTR ,NMDIM,4)
	 CALL ALLOC (PISLDR,NMDIM,4)
	 CALL ALLOC (PISLDR1,NMDIM,4)
	 CALL ALLOC (PXSLDR,NMDIM,8)
*
      ELSEIF (IMODE .EQ. 2) THEN
*
*   Double the allocation of storage for the arrays
*
         NEWSIZ = 2*NMDIM
*
         CALL RALLOC (PIPTR ,NMDIM,NEWSIZ,4)
         CALL RALLOC (PISLDR,NMDIM,NEWSIZ,4)
         CALL RALLOC (PISLDR1,NMDIM,NEWSIZ,4)
         CALL RALLOC (PXSLDR,NMDIM,NEWSIZ,8)
*
         NMDIM = NEWSIZ
*
      ELSEIF (IMODE .EQ. 3) THEN
*
*   Deallocate the storage for the arrays
*
	 CALL DALLOC (PIPTR )
	 CALL DALLOC (PISLDR)
	 CALL DALLOC (PISLDR1)
	 CALL DALLOC (PXSLDR)
*
      ELSE
*
         PRINT *, 'ALCNMA: Invalid argument IMODE = ',IMODE
         STOP
*
      ENDIF
*
      RETURN
      END
