      PROGRAM RWFNPLOT
      IMPLICIT NONE

      INTEGER, PARAMETER:: NPTS0=5000

      DOUBLE PRECISION pg(NPTS0), qg(NPTS0), rg(NPTS0), pgg(NPTS0)
      DOUBLE PRECISION rg2(NPTS0)
      DOUBLE PRECISION energy, a0
      CHARACTER        title*6, orbl*4, nnstr*2, new*3
      CHARACTER        name*100
      INTEGER          np, lp, jp, nn, laky, ll, jj, npts, j

      write(*,*) 'RWFNPLOT'
      write(*,*) 'Program to generate GNU Octave plot file' 
      write(*,*) 'Input file:  name.w'
      write(*,*) 'Output file: octave_name.m'
      write(*,*)
      write(*,*) 'To plot orbital: press enter'
      write(*,*) 'To remove orbital: type "d" or "D" and press enter'
      write(*,*)
      write(*,*) '                             Jorgen Ekman Aug 2014'
      write(*,*)

      WRITE(*,*) 'Name of state:'
      READ(*,*) name
      OPEN(3,FILE=trim(name)//'.w',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(4,FILE='octave_'//trim(name)//'.m',STATUS='UNKNOWN')
      write(4,*) '%     r                   P(r)                   Q(r)'
      write(4,*) 'clf'
      READ(3) title
      IF (title .NE. 'G92RWF') THEN   ! Extra safety
         PRINT *, 'title = ', title, 'does not match G92RWF'
         STOP
      ENDIF
      
      DO
         READ(3, END = 20) nn, laky, energy, npts
         write(nnstr,'(I2)') nn
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

         if(jj.eq.1) then
            if(ll.eq.0)  orbl = nnstr//'s '
            if(ll.eq.1)  orbl = nnstr//'p '
            if(ll.eq.2)  orbl = nnstr//'d '
            if(ll.eq.3)  orbl = nnstr//'f '
            if(ll.eq.4)  orbl = nnstr//'g '
            if(ll.eq.5)  orbl = nnstr//'h '
            if(ll.eq.6)  orbl = nnstr//'i '
            if(ll.eq.7)  orbl = nnstr//'j '
            if(ll.eq.8)  orbl = nnstr//'k '
            if(ll.eq.9)  orbl = nnstr//'l '
         else
            if(ll.eq.1)  orbl = nnstr//'p-'
            if(ll.eq.2)  orbl = nnstr//'d-'
            if(ll.eq.3)  orbl = nnstr//'f-'
            if(ll.eq.4)  orbl = nnstr//'g-'
            if(ll.eq.5)  orbl = nnstr//'h-'
            if(ll.eq.6)  orbl = nnstr//'i-'
            if(ll.eq.7)  orbl = nnstr//'j-'
            if(ll.eq.8)  orbl = nnstr//'k-'
            if(ll.eq.9)  orbl = nnstr//'l-'
         end if

         IF (npts .GT. NPTS0) THEN
            WRITE(*,*) 'df2hf: npts .GT. NPTS0'
            STOP
         ENDIF

         READ(3) a0, (pg(j), j=1,npts), (qg(j), j=1,npts)
         READ(3) (rg(j), j=1,npts)

         WRITE(*,'(2x,A,A,$)') orbl,' = '
         READ(*,'(A)') new
         IF ( NEW .NE. 'd  ' .AND. NEW .NE. 'D  ' ) THEN
            write(4,*) 'P = ['
            DO j = 1, npts
               rg2(j) = sqrt(rg(j))
               WRITE (4, '(3D20.10)') rg2(j), pg(j), qg(j)
            ENDDO
            write(4,*) '];'
            write(4,*) 'plot(P(:,1), P(:,2), ";',orbl,';")'
            write(4,*) 'hold all'            
         ENDIF
      ENDDO
   20 CONTINUE
      write(4,*) 'xlabel ("sqrt(r)", "fontsize", 12)' 
      write(4,*) 'ylabel ("P(r)", "fontsize", 12)' 
      write(4,*) 'grid on'
      PRINT *, ' FINISHED .....'
      STOP
      END
