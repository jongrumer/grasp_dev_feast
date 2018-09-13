program rmixreduce
  implicit none
  integer :: i,j,jj,k,kk,l,ios,err
  integer :: nelec, ncftot, nw, nvectot, nvecsize, nblock
  integer :: nb, nevblk(100), iatjp, iaspa
  integer :: ivec
  integer, allocatable :: checkcsf(:,:), ind(:,:), ind2(:,:), ncfblk(:)
  integer, allocatable :: ncfblockout(:)

  double precision :: eav,eval
  double precision :: dc2evec, totc2evecblk, c2evecblklim
  double precision, allocatable :: evec(:,:,:),c2evec(:,:)

  character(len=100) :: state, file1, file2
  character(len=600) :: header(5), string
  character(len=6) :: G92MIX
  character(len=1) :: ciflag
  character(len=100), allocatable  :: conf(:,:), coupling(:,:), spin(:,:)

  write(*,*) 'WELCOME TO PROGRAM RMIXREDUCE'
  write(*,*) 'Input files: <state>.(c)m, <state>.c'
  write(*,*) 'Reduced CSF list is written to rcsf.out'

  write(*,*) 
  write(*,*) 'Give name of the state: '
  read(*,*) state
  write(*,*) 'Expansion coefficients resulting from CI calculation (y/n)? '
  read(*,*) ciflag
  write(*,*) 'Fraction of total wave function [0-1] to be included in reduced list: '
  read(*,*) c2evecblklim

  l = index(state,' ')
  if(ciflag.eq.'y') then
     file1 = state(1:l-1)//'.cm'
  else
     file1 = state(1:l-1)//'.m'
  end if
  file2 = state(1:l-1)//'.c'

  ! Open and read data from mix file <state>.(c)m
  open (7,file=file1,status='old',form='unformatted')
  read (7,iostat=ios) G92MIX
  !write(*,*) G92MIX
  read (7) nelec, ncftot, nw, nvectot, nvecsize, nblock
  !write(*,*) '  nelec   = ', nelec
  !write(*,*) '  ncftot  = ', ncftot
  !write(*,*) '  nw      = ', nw
  !write(*,*) '  nblock  = ', nblock
  !write(*,*)

  ! Allocate various arrays
  allocate( conf(nblock,ncftot) )
  allocate( coupling(nblock,ncftot) )
  allocate( spin(nblock,ncftot) )
  allocate( ncfblockout(nblock) )
  allocate( ncfblk(ncftot) )
  allocate( ind(nblock,ncftot), ind2(nblock,ncftot) )
  allocate( checkcsf(nblock,ncftot) )
  allocate( evec(nblock,100,ncftot), c2evec(nblock,ncftot) )

  ! Continue to read data from mixing file <state>.(c)m
  ! Mixing coefficients are stored in 3D array: evec(block,eig,csf)
  evec(:,:,:) = 0.d0
  write(*,*)
  write(*,*) 'Block data read from mixing file'
  write(*,*) '        block        ncf         nev        2j+1          parity'
  do i=1, nblock
     READ (7,end=98) nb, ncfblk(i), nevblk(i), iatjp, iaspa
     write(*,*) nb, ncfblk(i), nevblk(i), iatjp, iaspa
     if(nevblk(i).gt.0) then
!        DO j = nvecpat + 1, nvecpat + nevblk(i)
!           ! ivec(i)   = ivec(i) + ncfpat ! serial # of the state
!           iatjpo(j) = iatjp
!           iaspar(j) = iaspa
!        ENDDO

        read (7) (ivec, j = 1,nevblk(i))
        !write(*,*) ivec(nvecpat+j)

        read (7) eav, (eval, j = 1, nevblk(i))
        !write(*,*) eav

        read (7) ((evec(i,k,j), j = 1, ncfblk(i)), k = 1, nevblk(i))
        !write(*,*) evec(nvecsizpat + ncfpat+j + (k-1)*ncftot)
     end if
     !nvecpat = nvecpat + nevblk(i)
     !ncfpat = ncfpat + ncfblk(i)
     !nvecsizpat = nvecsizpat + nevblk(i)*ncftot
  end do
98 continue
  close(7)

  ! For each CSF, calculate the sum of square of expansion coefficents for eigenvalues
  ! Block divided
  c2evec(:,:) = 0.0
  !write(*,*)
  do i=1,nblock
     !do j=1,2
        do k=1,ncfblk(i)
           do j=1,nevblk(i)
              c2evec(i,k) = c2evec(i,k)+evec(i,j,k)**2.0
           end do
           c2evec(i,k) = c2evec(i,k)/nevblk(i)
           !write(*,*) i, k, evec(i,1:nevblk(i),k),c2evec(i,k)
        end do
     !end do
  end do

  ! Sort square of expansions coefficients (c2evec) and index table (ind)
  ! Compute sum of square of expansion coefficients for each block
  ! and "flag" CSF:s that contribute to user defined fraction of total wave functions
  !write(*,*)
  checkcsf(:,:) = 0
  ind2(:,:) = 0
  do i = 1, nblock
     call HPSORT(ncfblk(i),c2evec(i,1:ncfblk(i)),ind(i,1:ncfblk(i)))
     totc2evecblk = 0.d0
     do k=ncfblk(i),1,-1
        totc2evecblk = totc2evecblk + c2evec(i,k)
        kk = ncfblk(i) + 1 - k
        ind2(i,kk) = ind(i,k)
        if(totc2evecblk.le.c2evecblklim) then
           checkcsf(i,ind2(i,kk)) = 1
        end if
        !write(*,*) k, kk, ind2(i,kk), c2evec(i,k), totc2evecblk, checkcsf(i,ind2(i,kk))
     end do
     !write(*,*)
  end do

  ! Open and read data from input file <state>.c
  open (7,file=file2,status='unknown',form='formatted')
  do j=1,5
     read(7,'(a)') header(j)
         !write(*,*) trim(header(j))
  end do
  i = 1
  j = 1
  do
     read(7,'(a)',end=99) string
     if(string(2:2).eq.'*') then
        !ncsf(i) = j - 1
        i = i + 1
        j = 1
        read(7,'(a)') conf(i,j)
     else
        conf(i,j) = string
     end if
     read(7,'(a)') coupling(i,j)
     read(7,'(a)') spin(i,j)
     j = j + 1
  end do
99 continue
  close(7)

  ! Open output file rcsf.out and write reduced CSF list 
  ncfblockout(:) = 0
  open (10,file='rcsf.out',status='unknown',form='formatted')
  do j=1,5
     write(10,'(a)') trim(header(j))
  end do
  do i = 1, nblock
     do j=1,ncfblk(i)
        if(checkcsf(i,j).eq.1) then
           ncfblockout(i) = ncfblockout(i) + 1
           !write(10,'(i10,i10)') j,ind2(i,j)
           write(10,'(a)') trim(conf(i,j))
           write(10,'(a)') trim(coupling(i,j))
           write(10,'(a)') trim(spin(i,j))
        end if
     end do
     if(i.lt.nblock) write(10,'(a2)') ' *'
  end do
  close(10)

  ! Print information of reduced list
  write(*,*)
  write(*,*) 'Number of CSF:s written to rcsf.out'
  write(*,*) '        block        ncf'
  do i = 1, nblock
     write(*,*) i, ncfblockout(i)
  end do
  
end program rmixreduce

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
10 continue
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
20 if(J.le.IR)then
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
