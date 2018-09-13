      PROGRAM plotmcdf
      IMPLICIT NONE

      INTEGER, PARAMETER:: NPTS0=500

      DOUBLE PRECISION pg(NPTS0), qg(NPTS0), rg(NPTS0), pgg(NPTS0)
      DOUBLE PRECISION energy, a0
      CHARACTER        title*6
      INTEGER          np, lp, jp, nn, laky, ll, jj, npts, j

      PRINT *, '****** PLOTMCDF ******'
      PRINT *, ' Convert one orbital to ASCII form for plotting'
      PRINT *
  123 PRINT *, ' Enter the n,  l,  J(1 or -1) '
      READ *, np, lp, jp
      IF (np. LE. lp .OR. lp .LT. 0 .OR. 
     &    .NOT. (jp .EQ. 1 .OR. jp .EQ. -1) ) THEN
         PRINT *, ' Inputs not acceptable, re-do ...'
         GOTO 123
      ENDIF

      OPEN(3,FILE='mcdf.w',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(4,FILE='mcdf.w.dat',STATUS='UNKNOWN')
      READ(3) title
      IF (title .NE. 'G92RWF') THEN   ! Extra safety
         PRINT *, 'title = ', title, 'does not match G92RWF'
         STOP
      ENDIF

      DO
         READ(3, END = 20) nn, laky, energy, npts
         PRINT *, nn, laky, energy, npts

         IF (laky .GT. 0) THEN
            ll = laky
            jj = -1
         ELSEIF (laky .LE. -1) THEN
            ll = -laky - 1
            jj = 1
         ELSE
            WRITE(*,*)'Unexpected case in reading mcdf.w'
            STOP
         ENDIF

         IF (npts .GT. NPTS0) THEN
            WRITE(*,*) 'df2hf: npts .GT. NPTS0'
            STOP
         ENDIF

         READ(3) a0, (pg(j), j=1,npts), (qg(j), j=1,npts)
         READ(3) (rg(j), j=1,npts)

         IF (nn .EQ. np .AND. ll .EQ. lp .AND. jj .EQ. jp) THEN
            DO j = 1, npts
               WRITE (4, '(3E15.5)') rg(j), pg(j), qg(j)
            ENDDO
            EXIT
         ENDIF
      ENDDO
   20 CONTINUE

      STOP
      END
