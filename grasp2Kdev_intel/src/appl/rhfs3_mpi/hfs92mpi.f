************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***           **   **  *******   *****    *****    *****             ***
***           **   **  **       **   **  **   **  **   **            ***
***           **   **  **       **       **   **      **             ***
***           *******  ****      *****    *****      **              ***
***           **   **  **            **      **     **               ***
***           **   **  **       **   **     **     **                ***
***           **   **  **        *****    **      *******            ***
***                                                                  ***
***            Relativistic Hyperfine Structure Program              ***
***                         GRASP92 Version                          ***
***                         Dynamic Storage                          ***
***                                                                  ***
***   ============================================================   ***
***   Copyright (c) 1995 by P Jonsson, F A Parpia, and C F Fischer   ***
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
      PROGRAM HFS92mpi
*                                                                      *
*   Entry routine for HFS92. Controls the entire computation.          *
*                                                                      *
*   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
*               [HFS92]: CHKPLT, GETHFD, HFS, SETDBG, SETSUM,          *
*                        STRSUM.                                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
cb adapted for MPI by Jacek Bieron
cb 14 Apr 2008
*                                                                      *
*   Modified by Gediminas Gaigalas for new spin-angular integrations.  *
*                                                   Vilnius May 2012   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN, YES
      CHARACTER*24 NAME
      CHARACTER*26 name_rwf
      COMMON/DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr, lenhost
      COMMON /mpi/ myid, nprocs, ierr
      CHARACTER host*(MPI_MAX_PROCESSOR_NAME), idstring*3

! timing
      INTEGER ncount1

! cpath uses
      INTEGER       lenperm, lentmp
      CHARACTER*128 startdir, permdir, tmpdir
!-----------------------------------------------------------------------

cbn workaround lam 7.1.1 i/o redirection bug
cbn
      CHARACTER*128 currentdir
      INTEGER    lencurrentdir
      Logical inputFileExists
cbn

*=======================================================================
*  Start mpi --- get processor info: myid, nprocs, host name; and print
*=======================================================================        
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,
     &                     startdir, permdir, tmpdir, 'hfs92MPI')
      WRITE (idstring, '(I3.3)') myid

      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)                                                

*=======================================================================        

      if (myid .eq. 0) then

      print*, ' tmpdir = ', tmpdir(1:lentmp)

cbn workaround lam 7.1.1 i/o redirection bug
cbn
        call sys_getwd(currentdir,ierr)
        lencurrentdir = len_trim(currentdir)
        call sys_chdir(permdir,lenperm,ierr)
         WRITE (istde,*) ' permdir = ', permdir(1:lenperm)

!       INQUIRE(file = 'hfs.input', exist=inputFileExists)
!       if(inputFileExists) then
!         OPEN (5, FILE ='hfs.input', STATUS = 'OLD', IOSTAT = ierr)
!         print*, ' inputFileExists = ', inputFileExists 
!       else
!        WRITE (istde,*) ' no hfs.input '
!       endif
c        OPEN (5, FILE ="xinput", STATUS = 'UNKNOWN', IOSTAT = ierr)

      WRITE (istde,*)
      WRITE (istde,*) 'hfs92MPI: Execution begins ...'

      WRITE (istde,*)
      WRITE (istde,*) 'Default settings ?'
cb     WRITE (istde,'(A)',ADVANCE='NO') 'Default settings?  (y/n) '
cb     print*, ' before GETYN = '
cb    YES = .false.
cb    print*, ' YES = ', YES
      YES = GETYN ()
cb     print*, ' after GETYN = '
cb     print*, ' YES = ', YES
      WRITE (istde,*)
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 WRITE (istde,*) 'Name of state'
      READ(*,'(A)') NAME
cb
      WRITE (istde,*) 'Name of state = ', NAME

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
      endif   !myid=0

      call sys_chdir(tmpdir,lentmp,ierr)
c     CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c     CALL MPI_Bcast (NCI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

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
      if (myid .EQ. 0) then

      call sys_chdir(permdir,lenperm,ierr)

      CALL SETSUM (NAME,NCI)

      endif   !myid=0
*
*   Open, check, load data from, and close, the  .csl  file
*
      CALL lodcslmpijb (NAME,ncore_not_used)

c     print *, ' after lodcslmpijb '
c
*   Open, check, load data from, and close the  .iso  file
*
      CALL SETISOmpi ('isodata')
*
c     print *, ' after SETISOmpi '
c
*   Get the remaining information
*
      CALL GETHFDmpi (NAME)
*
c     print *, ' after GETHFDmpi '
c
*   Load the radial wavefunctions
*
cbieron
cb
      NAME_len = len_trim(NAME)
      name_rwf = NAME(1:NAME_len)//'.w'

cb      WRITE (istde,*) name_rwf

      CALL SETRWFmpi (name_rwf)
*
c     print *, ' after SETRWFmpi ,myid = ', myid
c
*   Get the eigenvectors
*
      CALL readmixmpi (NAME,NCI)
*
c     print *, ' after readmixmpi '
c
*   Append a summary of the inputs to the  .sum  file
*
      if (myid .EQ. 0) then

      CALL STRSUM

      endif  !myid = 0
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the HFS calculation
*
      call sys_chdir(tmpdir,lentmp,ierr)

      CALL HFS(host)
*
*   Print completion message
*
      IF (myid .EQ. 0) THEN
        WRITE (istde,*)
        WRITE (istde,*) 'hfs92MPI:: Execution complete.'
      ENDIF

*=======================================================================
*  Execution finished; Statistics output
*=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost,
     &                     ncount1, 'hfs92MPI')
*
      STOP
      END
