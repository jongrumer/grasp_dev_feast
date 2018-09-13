************************************************************************
*                                                                      *
      SUBROUTINE ALCLLA (PIBEG,PILAB,PILAST,PILEFT,PIPTCS,
     :                   PIRIGH,PLBLIN,LLDIM,IMODE)
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
Cww      INTEGER PIBEG,PILAB,PILAST,PILEFT,PIPTCS,PIRIGH,PLBLIN
      POINTER (PIBEG,IBEGDUMMY)
      POINTER (PILAB,ILABDUMMY)
      POINTER (PILAST,ILASTDUMMY)
      POINTER (PILEFT,ILEFTDUMMY)
      POINTER (PIPTCS,IPTCSDUMMY)
      POINTER (PIRIGH,IRIGHDUMMY)                                             
      POINTER (PLBLIN,LBLINDUMMY)
*
      IF (IMODE .EQ. 1) THEN
*
*   Initial array dimension
*
         LLDIM = 64
*
*   Allocate storage for arrays
*
	 CALL ALLOC (PIBEG ,LLDIM,4)
	 CALL ALLOC (PILAB ,LLDIM,4)
	 CALL ALLOC (PILAST,LLDIM,4)
	 CALL ALLOC (PILEFT,LLDIM,4)
	 CALL ALLOC (PIPTCS,LLDIM,4)
	 CALL ALLOC (PIRIGH,LLDIM,4)
	 CALL ALLOC (PLBLIN,LLDIM,4)
*
      ELSEIF (IMODE .EQ. 2) THEN
*
*   Double the allocation of storage for the arrays
*
         NEWSIZ = 2*LLDIM
*
         CALL RALLOC (PIBEG ,LLDIM,NEWSIZ,4)
         CALL RALLOC (PILAB ,LLDIM,NEWSIZ,4)
         CALL RALLOC (PILAST,LLDIM,NEWSIZ,4)
         CALL RALLOC (PILEFT,LLDIM,NEWSIZ,4)
         CALL RALLOC (PIPTCS,LLDIM,NEWSIZ,4)
         CALL RALLOC (PIRIGH,LLDIM,NEWSIZ,4)
         CALL RALLOC (PLBLIN,LLDIM,NEWSIZ,4)
*
         LLDIM = NEWSIZ
*
      ELSEIF (IMODE .EQ. 3) THEN
*
*   Deallocate the storage for the arrays
*
	 CALL DALLOC (PIBEG )
	 CALL DALLOC (PILAB )
	 CALL DALLOC (PILAST)
	 CALL DALLOC (PILEFT)
	 CALL DALLOC (PIPTCS)
	 CALL DALLOC (PIRIGH)
	 CALL DALLOC (PLBLIN)
*
      ELSE
*
         PRINT *, 'ALCLLA: Invalid argument IMODE = ',IMODE
         STOP
*
      ENDIF
*
      RETURN
      END
