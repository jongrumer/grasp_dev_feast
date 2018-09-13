************************************************************************
*                                                                      *
      PROGRAM TRANSTABLE 
*                                                                      *
*   This program reads the output from biotra2 together with spectral  *
*   designation and energies to produce a latex table                  * 
*                                                                      *
*   Written by Per Jonsson, Malmo University, November 2010            *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      PARAMETER (NTRANS = 100000)
      INTEGER ILAB(NTRANS),GFPRINT,IU(NTRANS),IL(NTRANS),INDX(NTRANS)
      CHARACTER*80 ROW,FILENAME1,FILENAME2,LATEX(NTRANS),ULATEX(NTRANS)
      CHARACTER*80 LLATEX(NTRANS)
      CHARACTER*80 LATEXFILE1,LATEXFILE2
      CHARACTER*2 FU(NTRANS),FL(NTRANS),GAUGE(NTRANS),FLAB(NTRANS)
      CHARACTER*2 EM(NTRANS)
      CHARACTER*4 JU(NTRANS),JL(NTRANS),JLAB(NTRANS),PU(NTRANS)
      CHARACTER*4 PL(NTRANS)
      CHARACTER*1 ANS,PAR
      CHARACTER*14 FILE(NTRANS),FILE1,FILE2,FILEU(NTRANS),FILEL(NTRANS)
      CHARACTER*14 DUMMY
      CHARACTER*26 ROW1
      REAL*8 UENERGY(NTRANS),LENERGY(NTRANS),ENERGY(NTRANS)
      REAL*8 S(NTRANS),GF(NTRANS),DELTAE(NTRANS),R(NTRANS)
      REAL*8 DELTAA(NTRANS),AL(NTRANS),AV(NTRANS),TL(NTRANS),TV(NTRANS)
*
      WRITE(*,*) 
      WRITE(*,*) ' RTABTRANS2    '
      WRITE(*,*) ' This program reads energy label data and transition'     
      WRITE(*,*) ' data and creates LaTeX tables of transition data   '
      WRITE(*,*) ' and lifetime data.                                 '
      WRITE(*,*) ' Energy label data is given in the file energylabel '
      WRITE(*,*) ' created by the rtabtrans1 program           '
      WRITE(*,*) ' Transition data file can be conctenated *.t or *.ct'
      WRITE(*,*) ' files. '
      WRITE(*,*) ' Two LaTeX files are created: one with transition   '
      WRITE(*,*) ' data and one with lifetimes                        '
      WRITE(*,*) ' Input files: energylabel, transitiondatafile       '
      WRITE(*,*) ' Output files: transitiontable.tex,            '
      WRITE(*,*) '               lifetimetable.tex,            '

      WRITE(*,*)    
      WRITE(*,*) ' Give the name of the transition data file '
      READ(*,*) FILENAME1

!      WRITE(*,*) ' Give name of the energylabel file '
!      READ(*,*) FILENAME2
      FILENAME2 = 'energylabel'

!      WRITE(*,*) ' Give name of the output LaTeX transition file '
!      READ(*,*) LATEXFILE1
      LATEXFILE1 = 'transitiontable.tex'

!      WRITE(*,*) ' Give name of the output LaTeX lifetimes file '
!      READ(*,*) LATEXFILE2 
      LATEXFILE2 = 'lifetimetable.tex'
      
      WRITE(*,*) ' Give cut-off for printing A values '
   	READ(*,*) CUTOFF
      WRITE(*,*) ' Transition data wave length sorted? '
      READ(*,'(A)') ANS
      IF ((ANS.EQ.'y').OR.(ANS.EQ.'Y')) THEN
         NSORT = 1
      ELSE
         NSORT = 0
      END IF
      WRITE(*,*)

      OPEN(36,FILE=FILENAME1,FORM='FORMATTED',STATUS = 'OLD')
      OPEN(37,FILE=FILENAME2,FORM='FORMATTED',STATUS = 'OLD')
      OPEN(38,FILE=LATEXFILE1,FORM='FORMATTED',STATUS = 'UNKNOWN')
      OPEN(39,FILE=LATEXFILE2,FORM='FORMATTED',STATUS = 'UNKNOWN')

*---- Read energylabelfile --------------------

      DO I = 1,12
         READ(37,*)
      END DO

      I = 1
      NLEN = 0
      DO
        READ(37,200,END=91) NO,ILAB(I),JLAB(I),PAR,ENERGY(I),
     &        DENER,FILE(I),LATEX(I) 
        IF (LEN_TRIM(LATEX(I)).GT.NLEN) NLEN = LEN_TRIM(LATEX(I))  
        JLAB(I) = ADJUSTR(JLAB(I))
        I = I + 1
      END DO
      
   91 CONTINUE
      N = I - 1

*---- Check that there are no two quantum labels that are the same

      NQEQUAL = 0
      DO I = 1,N
         DO K = I+1,N
            IF (LATEX(I).EQ.LATEX(K)) THEN
               WRITE(*,*) 'Quantum labels for states',I,K,'equal'
               NQEQUAL = NQEQUAL + 1
            END IF
         END DO
      END DO
      IF (NQEQUAl.GT.0) THEN
         WRITE(*,*) 'Some quantum states have equal labels'
         WRITE(*,*) 'Update energylabel file and rerun'
         STOP
      END IF


*---- Read transition files ---------------------------------------

      K = 1   
      DO
   98   CONTINUE

*---- Start reading the file. Find out if transition between configurations
*     or within configuration

        READ(36,600) ROW1
        IF (LEN(TRIM(ROW1)).GT.22) THEN
          N12 = 2
          READ(36,500) FILE1
          READ(36,500) FILE2
    	  ELSE
          N12 = 1
          READ(36,501) FILE1
        END IF

*---- Loop over multipolarities in the file

        DO
  101     CONTINUE
 
*---- Find out if electric or magnetic case

          DO I = 1,3
            READ(36,'(A)') ROW
          END DO

          IF (ROW(2:9).EQ.'Electric') THEN
            NEM = 1
          ELSE
            NEM = 2
          END IF
   	    IF (ROW(16:16).EQ.'1') THEN
	         NEMO = 1
	       ELSE
	         NEMO = 2
	       END IF

*---- Continue to read until the rate information comes

          DO I = 1,4
            READ(36,'(A)') ROW
          END DO

*---- Start reading the rates for the found polarity (electric or magnetic)
          DO
            READ(36,300) FU(K),IU(K),JU(K),PU(K),FL(K),IL(K),JL(K),
     :	                PL(K),DELTAE(K),GAUGE(K),AV(K),GF(K),S(K)

*---- If magnetic then save rate both in AV(K) and AL(K) 
*     If electric AL(K) will be overwritten as the next line is read

            AL(K) = AV(K)
            R(K) = 0.D0

            IF (NEM.EQ.1) THEN
              READ(36,301) GAUGE(K),AL(K),GF(K),S(K)
              R(K) = dabs(AL(K)-AV(K))/maxval([AL(K),AV(K)])
*              write(*,*) R(K)
*              write(*,'(F11.2)') R(K) 
*              write(*,'(F11.3)') R(K) 
            END IF
            IF (FU(K).EQ.'f1') THEN
              FILEU(K) = FILE1 
              FILEL(K) = FILE2
	    ELSE IF (FU(K).EQ.'f2') THEN
	      FILEU(K) = FILE2
              FILEL(K) = FILE1
	    ELSE
	      FILEU(K) = FILE1
              FILEL(K) = FILE1
	    END IF 

*-----Set polarity for case

            IF ((NEM.EQ.1).AND.(NEMO.EQ.1)) THEN
		        EM(K) = 'E1'
	         ELSE IF ((NEM.EQ.1).AND.(NEMO.EQ.2)) THEN
              EM(K) = 'E2'
            ELSE IF ((NEM.EQ.2).AND.(NEMO.EQ.1)) THEN
              EM(K) = 'M1'
	         ELSE
	           EM(K) = 'M2'
	         END IF

*---- Try and read another line
*     Four cases occur: (1) end of file in which case we close the preset file
*                       (2) we can read a line and continue reading in the normal way
*                       (3) we read a blank line which indicates that we have one
*                            more polarity to read
*                       (4) we can read a file starting with 'Transition' in which case
*                           we have more to read

            READ(36,700,END=99) ROW
            K = K + 1
            IF (K.EQ.NTRANS) THEN
               WRITE(*,*) 'Too many transitions'
               WRITE(*,*) 'Increase NTRANS and recompile'
               STOP
            END IF
            BACKSPACE 36

*---- More polarity

            IF (LEN(TRIM(ROW)).EQ.0) THEN
			     GOTO 101
	         END IF

*---- More transitions but within new files

            IF (LEN(TRIM(ROW)).LT.30) THEN
			     GOTO 98
	         END IF

*---- If not the above cases just continue with another row in the normal way

          END DO
	     END DO
      END DO

   99 CONTINUE

*---- Write header to LaTeX table ----------------------------------

      write(38,'(A)') '\documentclass[10pt]{article}'
      write(38,'(A)') '\usepackage{longtable}'
      write(38,'(A)') '\begin{document}'      
	   WRITE(38,'(A)') '\begin{longtable}{lllrllll} \hline'
	   WRITE(38,'(A)') 'Upper & Lower & EM & ',
     :   '$\Delta E$ (cm$^{-1}$) & $\lambda$ (\AA) & ',
     :   '$A$ (s$^{-1}$) & $gf$ & $dT$ 
     :         \\ \hline'          

*---- Sort energies --------------------------------------------

      INDX = 0

      DO J = 1,N
        DO I = 1,N
          DO L = 1,K
            IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.
     :         (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J)).
     :          AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
                DELTAA(L) = 1.0E+8/DELTAE(L)
	         END IF              
          END DO
        END DO
      END DO	    

      call HPSORT(K,DELTAA,indx)


      IF (NSORT.EQ.0) THEN

*---- Write transition data to LaTeX file. ----------------------
*     Order according to the order in energy file

        DO J = 1,N
          DO I = 1,N
            DO L = 1,K
              IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.
     :          (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J)).
     :          AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
	             IF (AL(L).GT.CUTOFF) THEN
                  DELTAA(L) = 1.0E+8/DELTAE(L)
                  WRITE(38,410) LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN),
     :              ' & ',EM(L),' & ',INT(DELTAE(L)),' &',
     :              DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
	             END IF
	           END IF              
            END DO
          END DO
        END DO	    

      ELSE

*---- Write transition data to LaTeX file. ----------------------
*     Order according to wave length

        DO LL = 1,K
          L = INDX(LL)
          DO J = 1,N
            DO I = 1,N
              IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.
     :          (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J)).
     :          AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
	             IF (AL(L).GT.CUTOFF) THEN
                  DELTAA(L) = 1.0E+8/DELTAE(L)
                  WRITE(38,410) LATEX(J)(1:NLEN),' & ',LATEX(I)(1:NLEN),
     :              ' & ',EM(L),' & ',INT(DELTAE(L)),' &',
     :              DELTAA(L),' &',AL(L),' &',GF(L),' &',0.1*R(L),'\\'
	             END IF
	           END IF              
            END DO
          END DO
        END DO

      END IF

      write(38,'(A)') '\hline\\'                                                      
      write(38,'(A)') '\caption{Transition data}'
      write(38,'(A)') '\end{longtable}'
      write(38,'(A)') '\end{document}'

*---- Write header to LaTeX table ----------------------------------

      write(39,'(A)') '\documentclass[10pt]{article}'
      write(39,'(A)') '\usepackage{longtable}'
      write(39,'(A)') '\begin{document}'      
	   WRITE(39,'(A)') '\begin{longtable}{lll} \hline'
      WRITE(39,'(A)') 'State & $\tau_l$ & $\tau_v$  \\ \hline'

*---- Compute and print lifetimes for all the states ---------------

      DO J = 1,N
        TL(J) = 0.D0
        TV(J) = 0.D0
        DO I = 1,N
          DO L = 1,K
            IF ((FL(L) == FL(L)).AND.(IL(L) == ILAB(I)).
     :        AND.(JL(L) == JLAB(I)).AND.(FILEL(L) == FILE(I)).AND.
     :          (FU(L) == FU(L)).AND.(IU(L) == ILAB(J)).
     :        AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
              TL(J) = TL(J) + AL(L)
              TV(J) = TV(J) + AV(L)
            END IF
          END DO
        END DO
        IF (TL(J).GT.0.D0) THEN
          WRITE(39,420) LATEX(J)(1:NLEN),' & ',1.D0/TL(J),
     :                                   ' & ',1.D0/TV(J),'\\'
        END IF
      END DO 

      write(39,'(A)') '\hline\\'                                                      
      write(39,'(A)') '\caption{Lifetimes in s$^{-1}$.}'
      write(39,'(A)') '\end{longtable}'
      write(39,'(A)') '\end{document}'

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) ' Program finished. The transition tables in latex'
      WRITE(*,*) ' have been written to file '

 200  FORMAT(2I3,1X,A4,1x,A1,2X,F14.7,F12.2,2X,A14,6X,A)      
 300  FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,1D12.6,A2,1P,
     :     D13.5,2D12.4)
 301  FORMAT(41X,A2,1P,D13.5,2D12.4)
 410  FORMAT(6A,I11,A,F20.6,A,1P,E11.3,A,E11.3,A,F11.3,A)
 420  FORMAT(2A,1P,E11.4,A,E11.4,A)
 500  FORMAT(6X,A14)
 501  FORMAT(5X,A14)
 600  FORMAT(A26)
 700  FORMAT(A80)
      
      END PROGRAM TRANSTABLE 
      
      SUBROUTINE HPSORT(N,RA,IND)
      double precision RA(N)
      integer IND(N)
      L=N/2+1
      IR=N
      do I=1,N
         IND(I) = I
      end do

  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10    continue
      if(L > 1)then
         L=L-1
         RRA=RA(L)
         IRRA = IND(L)  ! je
      else
         RRA=RA(IR)
         IRRA = IND(IR) ! je
         RA(IR)=RA(1)
         IND(IR)=IND(1)  !je
         IR=IR-1
         if(IR.eq.1)then
            RA(1)=RRA
            IND(1)=IRRA !je
            return
         end if
      end if
      I=L
      J=L+L
20    if(J.le.IR)then
         if(J < IR)then
            if(RA(J) < RA(J+1))  J=J+1
         end if
         if(RRA < RA(J))then
            RA(I)=RA(J)
            IND(I)=IND(J) !je
            I=J; J=J+J
         else
            J=IR+1
         end if
     
         goto 20
      end if
      RA(I)=RRA
      IND(I)=IRRA
      goto 10
      END
