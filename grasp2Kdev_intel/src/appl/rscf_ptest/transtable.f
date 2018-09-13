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

*  Example in /home/per/graspruns/B-like/MnXXI/trans on GRODAN

      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      INTEGER ILAB(100000),GFPRINT,IU(100000),IL(100000),DEBUG
      CHARACTER*80 ROW,FILENAME1,FILENAME2,LATEX(100000),ULATEX(100000)
      CHARACTER*80 LLATEX(100000)
      CHARACTER*80 LATEXFILE1,LATEXFILE2
      CHARACTER*2 FU(100000),FL(100000),GAUGE(100000),FLAB(100000)
      CHARACTER*2 EM(100000)
      CHARACTER*4 JU(100000),JL(100000),JLAB(100000),PU(100000)
      CHARACTER*4 PL(100000)
      CHARACTER*1 ANS
      CHARACTER*14 FILE(100000),FILE1,FILE2,FILEU(100000),FILEL(100000)
      CHARACTER*14 DUMMY
      CHARACTER*26 ROW1
      REAL*8 UENERGY(100000),LENERGY(100000),ENERGY(100000),A(100000)
      REAL*8 S(100000),GF(100000),DELTAE(100000),T(100000),R(100000)
      REAL*8 DELTAA(100000)
*
      WRITE(*,*) 
      WRITE(*,*) ' Welcome to the transition table program that     '
      WRITE(*,*) ' converts transition tables to LaTeX              '     
      WRITE(*,*) 
      WRITE(*,*) ' Begin and concatenate all transition files to    '
      WRITE(*,*) ' one file referred to as THE transition file      '

      WRITE(*,*)
      WRITE(*,*) ' Debug (0/1) ? '
      READ(*,*) DEBUG

      WRITE(*,*)    
      WRITE(*,*) ' Give the name of the transition file '
      READ(*,*) FILENAME1

      WRITE(*,*) ' Give name of the energylabel file '
      READ(*,*) FILENAME2

      WRITE(*,*) ' Give name of the output LaTeX transition file '
      READ(*,*) LATEXFILE1

      WRITE(*,*) ' Give name of the output LaTeX lifetimes file '
      READ(*,*) LATEXFILE2 
      
      WRITE(*,*) ' Give cut-off for printing A values '
	READ(*,*) CUTOFF
      WRITE(*,*)

      WRITE(*,*) 'CUTOFF: ', CUTOFF

      OPEN(36,FILE=FILENAME1,FORM='FORMATTED',STATUS = 'OLD')
      OPEN(37,FILE=FILENAME2,FORM='FORMATTED',STATUS = 'OLD')
      OPEN(38,FILE=LATEXFILE1,FORM='FORMATTED',STATUS = 'UNKNOWN')
      OPEN(39,FILE=LATEXFILE2,FORM='FORMATTED',STATUS = 'UNKNOWN')

*---- Read file with energy labels and energies --------------------
      I = 1 
      DO 
        READ(37,200,END=91) FILE(I),ILAB(I),JLAB(I),LATEX(I)
        READ(37,*) ENERGY(I)
        IF (DEBUG.EQ.1) THEN
          WRITE(*,200) FILE(I),ILAB(I),JLAB(I),LATEX(I)
	    WRITE(*,*) ENERGY(I)
        END IF
        I = I + 1
      END DO
      
   91 CONTINUE
      N = I - 1
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
          READ(36,500) FILE1
        END IF

*---- Loop over multipolarities in the file

        DO
  101     CONTINUE
 
*---- Find out if electric or magnetic case

          DO I = 1,3
            READ(36,'(A)') ROW
          END DO
          IF (DEBUG.EQ.1) WRITE(*,'(A)') ROW

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
            IF (DEBUG.EQ.1) WRITE(*,'(A)') ROW
          END DO

*---- Start reading the rates for the found polarity (electric or magnetic)
          DO
            READ(36,300) FU(K),IU(K),JU(K),PU(K),FL(K),IL(K),JL(K),
     :	    PL(K),DELTAE(K),GAUGE(K),A(K),GF(K),S(K)
            R(K) = 1.D0
	    IF (NEM.EQ.1) THEN
              A0 = A(K)
              READ(36,301) GAUGE(K),A(K),GF(K),S(K)
              R(K) = A(K)/A0
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

            IF (DEBUG.EQ.1) THEN
              WRITE(*,*) 'K',K
              WRITE(*,*) FU(K),IU(K),JU(K),PU(K),FL(K),IL(K),JL(K),
     :          PL(K),DELTAE(K),GAUGE(K),A(K),GF(K),S(K),R(K)
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
            BACKSPACE 36

*---- More polarity

            IF (LEN(TRIM(ROW)).EQ.0) THEN
		    WRITE(*,*) 'More multipolarities'
			GOTO 101
	      END IF

*---- More transitions but within new files

            IF (LEN(TRIM(ROW)).LT.30) THEN
		    WRITE(*,*) 'Transitions in new files'
			GOTO 98
	      END IF

*---- If not the above cases just continue with another row in the normal way

          END DO
	  END DO
      END DO

   99 CONTINUE

*---- Write header to LaTeX table ----------------------------------

      WRITE(38,'(A)') '\begin{table}'
      WRITE(38,'(A)') '\caption{Unscaled transition values}'
	WRITE(38,'(A)') '\begin{tabular}{llrllll} \hline'
	WRITE(38,'(A)') 'Upper & Lower & ',
     :        '$\Delta E_{calc}$ & $\lambda$ (nm)& $A$ & $gf$ & $R$ 
     :         \\ \hline'          

*---- Write transition data to LaTeX file. ----------------------
*     Order according to the order in energy file

      WRITE(*,*) 'Ordering (0/1)'
      READ(*,*) NORDER

      IF (NORDER.EQ.0) THEN

        DO J = 1,N
          DO I = 1,N
            DO L = 1,K
              IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.
     :          (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J)).
     :          AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
	          IF (A(L).GT.CUTOFF) THEN
                     DELTAA(L) = 1.0E+7/DELTAE(L)
                     WRITE(38,410) TRIM(LATEX(J)),' & ',TRIM(LATEX(I)),
     :                    ' & ',INT(DELTAE(L)),' &',
     :                 DELTAA(L),' &',A(L),' &',GF(L),' &',R(L)/10,'\\'
	          END IF
	        END IF              
            END DO
          END DO
        END DO	    

      ELSE
        DO I = 1,N
          DO J = 1,N
            DO L = 1,K
              IF ((IL(L) == ILAB(I)).AND.(JL(L) == JLAB(I)).AND.
     :	      (FILEL(L) == FILE(I)).AND.(IU(L) == ILAB(J)).
     :          AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
                DELTAEEXP = ENERGY(J)-ENERGY(I)
                IF (A(L).GT.CUTOFF) THEN
			    WRITE(38,410) TRIM(LATEX(J)),' & ',TRIM(LATEX(I)),
     :			' & ',INT(DELTAEEXP),' &',
     :            EM(L),' &',GF(L),' &',A(L),' &',R(L)/10,'\\'
	          END IF
              END IF
            END DO
          END DO
        END DO

      END IF

      WRITE(38,*) '\end{tabular}'
      WRITE(38,*) '\end{table}'

*---- Write header to LaTeX table ----------------------------------

      WRITE(39,'(A)') '\begin{table}'
      WRITE(39,'(A)') '\caption{Lifetimes}'
      WRITE(39,'(A)') '\begin{tabular}{ll} \hline'
      WRITE(39,'(A)') 'State & $\tau$  \\ \hline'

*---- Compute and print lifetimes for all the states ---------------

      DO J = 1,N
        T(J) = 0.D0
        DO I = 1,N
          DO L = 1,K
            IF ((FL(L) == FL(L)).AND.(IL(L) == ILAB(I)).
     :        AND.(JL(L) == JLAB(I)).AND.(FILEL(L) == FILE(I)).AND.
     :          (FU(L) == FU(L)).AND.(IU(L) == ILAB(J)).
     :        AND.(JU(L) == JLAB(J)).AND.(FILEU(L) == FILE(J))) THEN
              T(J) = T(J) + A(L)
            END IF
          END DO
        END DO
        IF (T(J).GT.0.D0) THEN
          WRITE(39,420) TRIM(LATEX(J)),' & ',1.D0/T(J),'\\'
        END IF
      END DO 

      WRITE(39,*) '\end{tabular}'
      WRITE(39,*) '\end{table}'

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) ' Program finished. The transition tables in latex'
      WRITE(*,*) ' have been written to file '
      
  200 FORMAT(A14,I2,A4,A)
  300 FORMAT(1X,A2,I3,1X,2A4,A2,I3,1X,2A4,1D12.6,A2,1P,
     :     D13.5,2D12.4)
 301  FORMAT(41X,A2,1P,D13.5,2D12.4)
 410  FORMAT(4A,I11,A,F20.6,A,1P,E11.3,A,E11.3,A,F11.2,A)
 420  FORMAT(2A,E11.4,A)
 500  FORMAT(6X,A14)
 600  FORMAT(A26)
 700  FORMAT(A80)
      
      END PROGRAM TRANSTABLE 
      
