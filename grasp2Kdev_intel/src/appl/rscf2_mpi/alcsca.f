************************************************************************
*                                                                      *
      SUBROUTINE ALCSCA (PNTI,PNTR,IDIM,IMODE)
*                                                                      *
*   This subprogram allocates (IMODE = 1), reallocates  (IMODE = 2),   *
*   and deallocates (IMODE = 3) storage for arrays in COMMON/SCF2/.    *
*
*     imode = 1 then idim = 64 (out)
*     imode = 2 then idim = 2*idim (inout)
*     imode = 3 then dalloc
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*                                                                      *
*   Farid A. Parpia.                      Last revision: 17 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTI,PNTR
      POINTER (PNTI,IDUMMY)
      POINTER (PNTR,RDUMMY)
*
      IF (IMODE .EQ. 1) THEN
*
*   Initial array dimension
*
         IDIM = 64
*
*   Allocate storage for arrays
*
	 CALL ALLOC (PNTI,IDIM,4)
	 CALL ALLOC (PNTR,IDIM,8)
*
      ELSEIF (IMODE .EQ. 2) THEN
*
*   Double the allocation of storage for the arrays
*
         NEWSIZ = 2*IDIM
*
         CALL RALLOC (PNTI,IDIM,NEWSIZ,4)
         CALL RALLOC (PNTR,IDIM,NEWSIZ,8)
*
         IDIM = NEWSIZ
*
      ELSEIF (IMODE .EQ. 3) THEN
*
         CALL DALLOC (PNTI)
         CALL DALLOC (PNTR)
*
      ELSE
*
         PRINT *, 'ALCSCA: Invalid argument IMODE = ',IMODE
         STOP
*
      ENDIF
*
      RETURN
      END
