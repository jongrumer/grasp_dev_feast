**********************************************************************************
**********************************************************************************
**                                                                              ** 
**  **  ** *****  ****   ******* ***** ***** **     **       *       **     **  **
**  **  ** **    **  **      **  **    **    ***   ***      ***      ***    **  **
**  **  ** **    **         **   **    **    ** * * **     ** **     ** *   **  **
**  ****** ****   ****     **    ****  ****  **  *  **    **   **    **  *  **  **
**  **  ** **        **   **     **    **    **     **   *********   **   * **  **
**  **  ** **    **  **  **      **    **    **     **  **       **  **    ***  **
**  **  ** **     ****  *******  ***** ***** **     ** **         ** **     **  **
**                                                                              **
**                                                                              **
**                  Relativistic Hyperfine and Zeeman Program                   **
**                               GRASP2K Version                                **
**                               Dynamic Storage                                **
**                                                                              **
**  =========================================================================   **
**  Copyright (c) 2006 by M Andersson, P Jonsson                                **
**  =========================================================================   **
**   All rights reserved.  No part of this software or its accompanying docu-   **
**   mentation may be reproduced,  stored in a retrieval system, or transmit-   **
**   ted, in any form or by any  means, electronic, mechanical, photocopying,   **
**   recording,  or otherwise,  without the  prior written  permission of the   **
**   authors.                                                                   **
**                                                                              **
**                                Disclaimer                                    **
**                                ==========                                    **
**   The authors make no warranties,  express or implied,  that this software   **
**   or its accompanying  documentation are  free of error or that  they will   **
**   meet your  requirements  for any  particular  application.  Neither this   **
**   software  nor its accompanying  documentation should  be relied upon for   **
**   solving problems if an incorrect solution could result in injury or loss   **
**   of, or damage to, property. If you use this software or its accompanying   **
**   documentation, you do so entirely at your own risk; the authors disclaim   **
**   all liability for direct or consequential damage.                          **
**                                                                              **
**********************************************************************************
**********************************************************************************
*                                                                                *
      PROGRAM HFSZEEMAN06
*                                                                                *
*   Entry routine for HFSZEEMAN. Controls the entire computation.                *
*                                                                                *
*   Call(s) to: [LIB92]: GETMIXA, GETMIXC, SETCSL, SETMC, SETCON.                *
*               [HFSZEEMAN06]: CHKPLT, GETHFD, GETMIXBLOCK, HFSZEEMAN, SETDBG,   *
*                               SETSUM, STRSUM.                                  *
*               [NJGRAF]: FACTT.                                                 *
*                                                                                *
*      M. Andersson and P Jönsson                               2006             *
*                                                                                *
**********************************************************************************
*
      LOGICAL GETYN, YES
      CHARACTER*24 NAME
      COMMON/DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde

      WRITE (istde,*)
      WRITE (istde,*) 'HFSZEEMAN: Execution begins ...'

      WRITE (istde,*)
      WRITE (istde,*) 'Default settings?'
      YES = GETYN ()
      WRITE (istde,*)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 WRITE (istde,*) 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         WRITE (istde,*) 'Names may not start with a blank'
         GOTO 10
      ENDIF
      WRITE (istde,*)
      WRITE (istde,*) 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
      WRITE (istde,*)
      WRITE (istde,*) 'Calculate off-diagonal matrix elements? '
      YES = GETYN ()
      IF (YES) THEN
         NOFFD = 0
      ELSE
         NOFFD = 1
      ENDIF
*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT
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
*   Set up the physical constants
*
      CALL SETCON
*
*   Open the  .sum  file
*
      CALL SETSUM(NAME,NCI,NOFFD)
*
*   Open, check, load data from, and close, the  .csl  file
*
      CALL SETCSLA(NAME,ncore_not_used)
*
*   Get the remaining information
*
      CALL GETHFD(NAME)
*
*   Get the eigenvectors
*
      WRITE(istde,*) 'Block format?'
      YES = GETYN ()
      WRITE (istde,*)
      IF (YES) THEN
         CALL GETMIXBLOCK(NAME,NCI)
      ELSE
         IF (NCI.EQ.0) THEN
            CALL GETMIXC(NAME)
         ELSE
            CALL GETMIXA(NAME)
         ENDIF
      ENDIF
*
*   Append a summary of the inputs to the  .sum  file
*
      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the HFSZEEMAN calculation
*
      CALL HFSZEEMAN(NOFFD)
*
*   Print completion message
*
      WRITE (istde,*)
      WRITE (istde,*) 'HFSZEEMAN: Execution complete.'
*
      STOP
      END
