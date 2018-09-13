************************************************************************
*                                                                      *
      SUBROUTINE ALCNTA (PISLDR,PISLDR1,PXSLDR,NTDIM,IMODE)
*                                                                      *
*   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
*   and deallocates (IMODE = 3) storage for  certain arrays that are   *
*   used in OSCL.                                                      *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*                                                                      *
*   Farid A. Parpia.                      Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PISLDR,PXSLDR
      POINTER (PISLDR,ISLDRDUMMY)
      POINTER (PISLDR1,ISLDR1DUMMY)
      POINTER (PXSLDR,XSLDRDUMMY)
*
      IF (IMODE .EQ. 1) THEN
*
*   Initial array dimension
*
         NTDIM = 1
*
*   Allocate storage for arrays
*
	 CALL ALLOC (PISLDR,NTDIM,4)
	 CALL ALLOC (PISLDR1,NTDIM,4)
	 CALL ALLOC (PXSLDR,NTDIM,8)
*
      ELSEIF (IMODE .EQ. 2) THEN
*
*   Double the allocation of storage for the arrays
*
         NEWSIZ = 2*NTDIM
*
	 CALL RALLOC (PISLDR,NTDIM,NEWSIZ,4)
	 CALL RALLOC (PISLDR1,NTDIM,NEWSIZ,4)
	 CALL RALLOC (PXSLDR,NTDIM,NEWSIZ,8)
*
         NTDIM = NEWSIZ
*
      ELSEIF (IMODE .EQ. 3) THEN
*
*   Deallocate the storage for the arrays
*
	 CALL DALLOC (PISLDR)
	 CALL DALLOC (PISLDR1)
	 CALL DALLOC (PXSLDR)
*
      ELSE
*
         PRINT *, 'ALCNTA: Invalid argument IMODE = ',IMODE
         STOP
*
      ENDIF
*
      RETURN
      END
