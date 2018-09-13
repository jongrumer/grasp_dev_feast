************************************************************************
*                                                                      *
      SUBROUTINE RALC2D (PTR,ICDIMO,IRDIMO,ICDIMN,IRDIMN,LENGTH)
*                                                                      *
*   This  subprogram reallocates  the memory for the pointee of PTR:   *
*   This is assumed to be a  two-dimensional  array of data types of   *
*   length LENGTH (4 or 8) bytes.  ICDIMO is the column dimension on   *
*   entry;  ICDIMN is that on exit.  IRDIMO is the row  dimension on   *
*   entry; IRDIMN is that on exit.  This version is specific to Cray   *
*   CFT77  (UNICOS),  IBM XL FORTRAN (AIX/6000),  or SUN  FORTRAN 77   *
*   (Sun OS), depending on the action of the preprocessor PREPRO.      *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
*                                                                      *
*   Written by Farid A. Parpia.           Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
Cww      INTEGER PTR
      POINTER (PTR,PTRDUMMY)
*
      INTEGER*8 NBYTE
      POINTER (IPTR,IA(ICDIMN,*))
      POINTER (ITPTR,ITA(ICDIMO,*))
      POINTER (RPTR,RA(ICDIMN,*))
      POINTER (RTPTR,RTA(ICDIMO,*))
      COMMON/iounit/istdi,istdo,istde
*
*     IF (ICDIMN*IRDIMN*LENGTH .LE. 0) THEN
      NBYTE=ICDIMN
      NBYTE=NBYTE*IRDIMN
      NBYTE=NBYTE*LENGTH
      IF (NBYTE .LE. 0) THEN
*
*   Trap an invalid memory request
*
	 WRITE(istde,*) 'RALC2D: Invalid memory request:'
	 WRITE(istde,*) ' ICDIMN = ',ICDIMN,', IRDIMN = ',IRDIMN,
     :                                ', LENGTH = ',LENGTH,'.'
	 STOP
*
      ENDIF
*
*   How many columns and rows are to be copied?
*
      ICMIN = MIN (ICDIMO,ICDIMN)
      IRMIN = MIN (IRDIMO,IRDIMN)
*
      IF (LENGTH .EQ. 8) THEN
*
*   Arrays of length 8 bytes are assumed to be REAL*8         ;
*   these are treated in the first branch
*
	 RTPTR = PTR
	 CALL ALLOC (RPTR,ICDIMN*IRDIMN,LENGTH)
*
         DO 2 J = 1,IRMIN
	    DO 1 I = 1,ICMIN
	       RA(I,J) = RTA(I,J)
    1       CONTINUE
    2    CONTINUE
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
	 CALL ALLOC (IPTR,ICDIMN*IRDIMN,LENGTH)
*
*   Copy the contents of ITA into IA
*
         DO 4 J = 1,IRMIN
	    DO 3 I = 1,ICMIN
	       IA(I,J) = ITA(I,J)
    3       CONTINUE
    4    CONTINUE
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
         WRITE(istde,*) ' RALC2D: Invalid argument LENGTH = ',LENGTH,'.'
         STOP
*
      ENDIF
*
      RETURN
      END
