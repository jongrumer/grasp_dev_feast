************************************************************************
*                                                                      *
      SUBROUTINE RALLOC (PTR,OLDSIZ,NEWSIZ,LENGTH)
*                                                                      *
*   This  subprogram reallocates  the memory for the pointee of PTR.   *
*   OLDSIZ is the number of storage locations of length LENGTH (4 or   *
*   8) bytes that are associated with  PTR upon entry. NEWSIZ is the   *
*   number of storage locations associated with  PTR  on exit.  This   *
*   version is specific to  SUN FORTRAN 77  (Sun  OS) or  Cray CFT77   *
*   (UNICOS)  depending on the action of the preprocessor PREPRO.      *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This revision by Farid A. Parpia.     Last revision: 17 Sep 1992   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER OLDSIZ
      COMMON/iounit/istdi,istdo,istde
*
*      INTEGER*8 PTR
      POINTER (PTR,ptrdummy)
*
      POINTER (IPTR,IA(*)),(ITPTR,ITA(*))
      POINTER (RPTR,RA(*)),(RTPTR,RTA(*))
*
      IF (NEWSIZ*LENGTH .LE. 0) THEN
	 WRITE(istde,*) 'RALLOC: Invalid memory request:'
	 WRITE(istde,*) ' NEWSIZ = ',NEWSIZ,', LENGTH = ',LENGTH,'.'
	 STOP
      ENDIF
*
*   M is the number of array elements that are to be copied into the
*   new array
*
      M = MIN (OLDSIZ,NEWSIZ)
*
      IF (LENGTH .EQ. 8) THEN
*
*   Arrays of length 8 bytes are assumed to be REAL*8          ;
*   these are treated in the first branch
*
	 RTPTR = PTR
	 CALL ALLOC (RPTR,NEWSIZ,LENGTH)
*
	 DO 1 I = 1,M
	    RA(I) = RTA(I)
    1    CONTINUE
*
*   Reset PTR so that it is associated with the new array
*
         PTR = RPTR
*
*   Deallocate the array that used to be associated with PTR
*
         CALL DALLOC (RTPTR)
*
      ELSEIF (LENGTH .EQ. 4) THEN
*
*   Integer arrays are treated in the second branch; the procedure is
*   essentially the same as in the first branch
*
*   Associate ITA with PTR
*
	 ITPTR = PTR
*
*   Associate IA with IPTR
*
	 CALL ALLOC (IPTR,NEWSIZ,LENGTH)
*
*   Copy the contents of ITA into IA
*
	 DO 2 I = 1,M
	    IA(I) = ITA(I)
    2    CONTINUE
*
*   Reset PTR so that it is associated with the new array
*
         PTR = IPTR
*
*   Release ITA
*
         CALL DALLOC (ITPTR)
*
      ELSE
*
         WRITE(istde,*) 'RALLOC: Invalid argument LENGTH = ',LENGTH,'.'
         STOP
*
      ENDIF
*
      RETURN
      END
