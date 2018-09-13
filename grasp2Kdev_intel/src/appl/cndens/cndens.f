************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***        *****   **   **  ******   *******  **   **   *****        ***
***       **   **  ***  **  **   **  **       ***  **  **   **       ***
***       **       ***  **  **   **  **       ***  **  **            ***
***       **       ** * **  **   **  ****     ** * **   *****        ***
***       **       **  ***  **   **  **       **  ***       **       ***
***       **   **  **   **  **   **  **       **   **  **   **       ***
***        *****   **   **  ******   *******  **   **   *****        ***
***                                                                  ***
***    Program for `condensing' Configuration Symmetry List Files    ***
***                         GRASP92 Version                          ***
***                         Dynamic Storage                          ***
***                                                                  ***
***   ============================================================   ***
***   Copyright (c) 1992 by F A Parpia, I P Grant, and C F Fischer   ***
***   ============================================================   ***
***   All rights reserved.  No part of this software or its accom-   ***
***   panying documentation may be reproduced, stored in a retrie-   ***
***   val system,  or transmitted,  in any form or  by any  means,   ***
***   electronic, mechanical,  photocopying,  recording, or other-   ***
***   wise, without the prior written permission of the authors.     ***
***                                                                  ***
***                           Disclaimer                             ***
***                           ==========                             ***
***   The  authors make  no warranties,  express or implied,  that   ***
***   this software or its  accompanying documentation are free of   ***
***   error or that  they will meet your requirements for any par-   ***
***   ticular application.  Neither this software nor its accompa-   ***
***   nying documentation should be relied  upon for solving prob-   ***
***   lems if an incorrect solution could result in injury or loss   ***
***   of, or damage to, property. If you use this software or  its   ***
***   accompanying documentation,  you do so entirely  at your own   ***
***   risk;  the authors disclaim all liability for direct or con-   ***
***   sequential damage.                                             ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM CNDENS
*                                                                      *
*   Entry routine for CNDENS. Controls the entire computation.         *
*                                                                      *
*   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC.                        *
*               [MRGCSL]: CHKPLT, CUTOUT, SETDBG.                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN, YES
      CHARACTER*24 NAME
      COMMON/DEFAULT/NDEF
      PARAMETER (nblk0 = 20)
      CHARACTER*8 idblk(nblk0)
      common /BLK/NBLOCK,NCFBLK(30)

      write(*,*)
      write(*,*) 'RMIXCONDENS'
      write(*,*) 'Condens name.c and name.(c)m by deleting CSFs below'
      write(*,*) 'a given cut-off in all eigenvectors'
      write(*,*) 'Input files:  name.c, name.(c)m'
      write(*,*) 'Output files: name_cond.c, name_cond.(c)m '
      write(*,*)

      PRINT *, 'Default settings?'
      YES = GETYN ()
      IF (YES) THEN
        NDEF = 0
      ELSE
        NDEF = 1
      ENDIF

   10 PRINT *, 'Name of state'
      READ(*,'(A)') NAME
      print *, name
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
        PRINT *, 'Names may not start with a blank'
        GOTO 10
      ENDIF
      PRINT *, 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
        NCI = 0
      ELSE
        NCI = 1
      ENDIF
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Load the  .csl  file
*
      CALL SETCSLB(NAME,NCORE)
c     CALL SETCSL(TRIM(NAME)//'.c',NCORE,NBLK0,IDBLK)
c     CALL SETCSLA(NAME,NCORE)
c     CALL SETCSLA(NAME)
*
*   Load the  .mix  file
*
c     IF (NCI.EQ.0) THEN
c       CALL GETMIXC(NAME)
c     ELSE
c       CALL GETMIXA(NAME)
c     ENDIF
*
*   `Condense the  .csl  list, eliminating CSFs with sufficiently
*   small weights
*
      CALL CUTOUT(NAME,NCI)
*
*   Print completion message
*
      PRINT *
      PRINT *, 'RMIXCONDENS: Execution complete.'
*
      STOP
      END
