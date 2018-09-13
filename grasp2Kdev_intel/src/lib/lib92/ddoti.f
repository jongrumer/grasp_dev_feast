      REAL*8           function   ddoti   ( nz, x, indx, y )
C
C     ==================================================================
C     ==================================================================
C     ====  DDOTI -- REAL*8           INDEXED DOT PRODUCT           ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DDOTI COMPUTES THE VECTOR INNER PRODUCT OF 
C             A REAL*8           SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A REAL*8           VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         DDOTI   DOUBLE      REAL*8           FUNCTION VALUE EQUAL TO
C                             THE VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  DDOTI IS SET TO ZERO.
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
      REAL*8              X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      DDOTI = 0.0D0
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          DDOTI = DDOTI  +  X(I) * Y(INDX(I))
   10 CONTINUE
C
      RETURN
      END
