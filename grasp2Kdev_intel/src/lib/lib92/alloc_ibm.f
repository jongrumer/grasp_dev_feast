************************************************************************
*                                                                      *
      SUBROUTINE ALLOC (PTR,NMLOCS,LENGTH)
*                                                                      *
*   This  routine allocates  NMLOCS memory locations, each of length   *
*   LENGTH bytes to the pointee of the pointer PTR.  This version is   *
*   specific to  Sun FORTRAN 77 (Sun OS) or CFT77 (UNICOS) depending   *
*   ing upon the action of the preprocessor  PREPRO.  LENGTH is  not   *
*   required for the  CFT77  version --- the  8-byte word  length is   *
*   assumed common to all data types; it has been retained to ensure   *
*   compatibility with the code that calls this subprogram.            *
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This version by Farid A. Parpia.      Last revision: 17 Sep 1992   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*                                                                      *
************************************************************************
*
c
c
cbieron DEC pointer
c
c      INTEGER PTR                 ! DEC and CRAY require ?
      pointer (PTR,ptrdummy)
c
*
*   Compute the size of the memory region in bytes; this is retained for
*   calling compatibility and debugging purposes in the Cray version
*
      NBYTES = NMLOCS*LENGTH
*
      IF (NBYTES. LE. 0) THEN
         PRINT *, 'ALLOC: Invalid memory request:'
         PRINT *, ' PTR = ',PTR,', NMLOCS = ',NMLOCS,
     :               ', LENGTH = ',LENGTH,'.'
         STOP
*
      ELSE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         PTR = MALLOC (NBYTES)                ! SUN & DEC
!         PTR = MALLOC(%VAL(NBYTES))           ! IBM
!         IF (PTR .EQ. 0) THEN
!
!         IERR = MZMALLOC(PTR,NMLOCS)          ! TSS
!         IABORT = 0
!         CALL HPALLOC(PTR,NMLOCS,IERR,IABORT) ! CRAY 
!         IF (ierr.ne.0) THEN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Choose one from above

         PTR = MALLOC(%VAL(NBYTES))           ! IBM
         IF (PTR .EQ. 0) THEN

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            PRINT *, 'ALLOC: Unable to allocate memory:'
            PRINT *, ' PTR = ',PTR,', NMLOCS =',NMLOCS,
     :                  ', LENGTH = ',LENGTH,
     :                  '.'
            STOP
         ENDIF
      ENDIF
      RETURN
      END
