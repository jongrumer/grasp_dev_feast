************************************************************************
*                                                                      *
      SUBROUTINE ALCNSA (PNTJJA,PNTJJB,PNTHB1,PNTHB2,
     :                   PNTHC1,PNTHC2,PNTHM1,PNTHM2,
     :                   PNTLAB,PNNPTR,NSDIM,IMODE)
*                                                                      *
*   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
*   and deallocates (IMODE = 3) storage for  certain arrays that are   *
*   used in 12_oscl.                                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
*                                                                      *
*   Farid A. Parpia.                      Last revision: 28 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTJJA,PNTJJB,PNTHB1,PNTHB2,PNTHC1,PNTHC2,PNTHM1,PNTHM2,
Cww     :        PNTLAB,PNNPTR
      POINTER (PNTJJA,JJADUMMY)
      POINTER (PNTJJB,JJBDUMMY)
      POINTER (PNTHB1,HB1DUMMY)
      POINTER (PNTHB2,HB2DUMMY)
      POINTER (PNTHC1,HC1DUMMY)
      POINTER (PNTHC2,HC2DUMMY)
      POINTER (PNTHM1,HM1DUMMY)
      POINTER (PNTHM2,HM2DUMMY)
      POINTER (PNTLAB,LABDUMMY)
      POINTER (PNNPTR,NPTRDUMMY)
*
      IF (IMODE .EQ. 1) THEN
*
*   Initial array dimension
*
         NSDIM = 1
*
*   Allocate storage for arrays
*
	 CALL ALLOC (PNTJJA,NSDIM,4)
	 CALL ALLOC (PNTJJB,NSDIM,4)
	 CALL ALLOC (PNTHB1,NSDIM,8)
	 CALL ALLOC (PNTHB2,NSDIM,8)
	 CALL ALLOC (PNTHC1,NSDIM,8)
	 CALL ALLOC (PNTHC2,NSDIM,8)
	 CALL ALLOC (PNTHM1,NSDIM,8)
	 CALL ALLOC (PNTHM2,NSDIM,8)
	 CALL ALLOC (PNTLAB,NSDIM,4)
	 CALL ALLOC (PNNPTR,NSDIM,4)
*
      ELSEIF (IMODE .EQ. 2) THEN
*
*   Double the allocation of storage for the arrays
*
         NEWSIZ = 2*NSDIM
*
	 CALL RALLOC (PNTJJA,NSDIM,NEWSIZ,4)
	 CALL RALLOC (PNTJJB,NSDIM,NEWSIZ,4)
	 CALL RALLOC (PNTHB1,NSDIM,NEWSIZ,8)
	 CALL RALLOC (PNTHB2,NSDIM,NEWSIZ,8)
	 CALL RALLOC (PNTHC1,NSDIM,NEWSIZ,8)
	 CALL RALLOC (PNTHC2,NSDIM,NEWSIZ,8)
	 CALL RALLOC (PNTHM1,NSDIM,NEWSIZ,8)
	 CALL RALLOC (PNTHM2,NSDIM,NEWSIZ,8)
	 CALL RALLOC (PNTLAB,NSDIM,NEWSIZ,4)
	 CALL RALLOC (PNNPTR,NSDIM,NEWSIZ,4)
*
         NSDIM = NEWSIZ
*
      ELSEIF (IMODE .EQ. 3) THEN
*
*   Deallocate the storage for the arrays
*
	 CALL DALLOC (PNTJJA)
	 CALL DALLOC (PNTJJB)
	 CALL DALLOC (PNTHB1)
	 CALL DALLOC (PNTHB2)
	 CALL DALLOC (PNTHC1)
	 CALL DALLOC (PNTHC2)
	 CALL DALLOC (PNTHM1)
	 CALL DALLOC (PNTHM2)
	 CALL DALLOC (PNTLAB)
	 CALL DALLOC (PNNPTR)
*
      ELSE
*
         PRINT *, 'ALCNSA: Invalid argument IMODE = ',IMODE
         STOP
*
      ENDIF
*
      RETURN
      END
