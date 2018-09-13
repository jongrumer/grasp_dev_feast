      subroutine   daxpyi   ( nz, a, x, indx, y )
C
C     ==================================================================
C     ==================================================================
C     ====  DAXPYI -- INDEXED REAL*8           ELEMENTARY           ====
C     ====            VECTOR OPERATION                              ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DAXPYI ADDS A REAL*8           SCALAR MULTIPLE OF 
C             A REAL*8           SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         TO  
C             A REAL*8           VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST BE 
C         DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         A       DOUBLE      SCALAR MULTIPLIER OF  X.
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C
C     UPDATED ...
C
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ON OUTPUT
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C         
      REAL*8              Y (*), X (*), A
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( A .EQ. 0.0D0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I))  = Y(INDX(I)) + A * X(I)
   10 CONTINUE
C
      RETURN
      END
