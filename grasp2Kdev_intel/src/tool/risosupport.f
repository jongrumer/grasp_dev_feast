*
************************************************************************
*
      program risosupport
* 
* written by Cedric Naz\'e and Michel Godefroid
* version  October 13, 2012
* Copyright CompAS 
*                                                                      *
************************************************************************

      implicit double precision(a-h,o-z)

      parameter (emass = 5.4857990946d-4, ry = 10973731.568539d0)
      parameter (a0 = 0.52917721092d-10, pi = 3.14159265359d0)
      doubleprecision isomass(10),el(10),eu(10),diff(10),eln(10),eun(10)
     :,rme(10),muovm(10),rydm(10),isomassau(10)
      doubleprecision wlength(10),dwlength(10),rms(10)
      doubleprecision elns(10),euns(10)
      DOUBLE PRECISION K_NMSlow, K_NMSup, K_SMSlow, K_SMSup,E_FSlow,
     & E_FSup,mostabms,mostabfs
      INTEGER trans,theone
      character*5 isolabel(10), mostab
      character*72 transition
      CHARACTER*1 ans

      open(unit=13,file='ISTrans',status='unknown')
!      open(unit=14,file='isodata',status='unknown')
      theone=0
      convmkmhz = 29.9792458000000
      convmhzau = 1.5198298460045*1.d-10
      eh2ev   = 27.21138505
      eh2Mhz  = 6.579683920729*1.d9
      eh2cmM1 = 2.194746313708*1.d5
      cmM12eh = 4.556335252760*1.d-6
!      eh2m = 
      elec = 1.602176487*1.d-19
      clum = 299792458
      ukg  = 1.660538782*1.d-27

      write(*,*) 
      write(*,*) '==================================================='
      write(*,*) ' Program for calculating transition isotope shift.'
      write(*,*) '==================================================='
      write(*,*) 
*      write(*,*) '                                           _______'
     
Cmrg  write(*,*) ' of the nucleus should be given, as separated from '
Cmrg  write(*,*) ' the equivalent uniform nuclear radius Req where   '
Cmrg  write(*,*) ' <r^2> = (3/5)*Req^2 '
CedN  write(*,*) 
Cmrg  write(*,*) ' Note that the elctron density |wfn(o)|^2 should '
Cmrg  write(*,*) ' be given, as separated from the modified density '
Cmrg  write(*,*) ' D = 4*pi*|wfn(0)|^2'
      write(*,*) '==================================================='



      write(13,*)
      write(13,*) '=============='
      write(13,*) ' Initial data'
      write(13,*) '=============='

***********************************************************************
***********************************************************************
      write(*,*) 
*      write(*,*) 
*      write(*,*) 
      write(*,*) ' Specify the transition'
      read(*,'(a72)') transition

      write(13,*)
      write(13,*) ' Transition'
      write(13,'(a72)') transition

      write(*,*) ' Nuclear charge'
      read(*,*) z

      write(13,*)
      write(13,*) ' Nuclear charge',z

CedN      write(13,*)
CedN      write(*,*) ' Energy (a.u.) infinite nucl. mass for lower state'
CedN      read(*,*) elow
CedN      write(13,5003) 'lower',elow

CedN991   write(*,*) 
CedN     1' Energy (a.u.) infinite nucl. mass for upper state (1) or '
CedN      write(*,*) 
CedN     2' energy diff. (cm^-1) between upper and lower state (2)'
CedN      read(*,*) n
CedN      if (n.eq.1) then
CedN        write(*,*) ' Energy (a.u.) infinite nucl. mass for upper state'
CedN        read(*,*) eup
CedN      elseif (n.eq.2) then
        DO
!          write(*,*) 'Please type 1 or 2 to choose units'
          write(*,*) ' What is the refernce isotope? (ex: Nd142 )'
          read (*,*) mostab
          write(*,*) ' Introduce now the transition of the',
     & ' reference isotope'
          write(*,*) ' Please type'!  1 or 2 for  units'
          write(*,*) ' 1 for transition wavenumber (cm^-1)'
          write(*,*) ' 2 for transition wavelength (Angstrom)'
          read(*,*) trans
          IF ((trans .eq.1) .OR. (trans .eq.2)) EXIT
        ENDDO
        IF (trans.eq.1) THEN
          write(*,*) ' Please introduce the value in cm^-1:'
          read(*,*) ediffcm
          ediffwave = 1.0d0/(ediffcm*1.0d2)*1.0d10
          write(*,4008) ediffcm,' cm^-1.   '
!          write(*,*) ediffwave
        ELSEIF(trans.eq.2) THEN
          write(*,*) ' Please introduce the value in Angstrom:'
          read(*,*) ediffwave
          ediffcm = 1.0d-2/(ediffwave*1.0d-10)
          write(*,4008) ediffwave,' Angstrom.'
!          write(*,*) edffcm
        ENDIF
CedN        eup = elow + 50.d0*ediff/ry
CedN      else
CedN         goto 991
CedN      endif
!      if (n.eq.1) then
!         write(13,*) ' Energy for upper state from calculation'
!      else
!         write(13,*) ' Energy for upper state from exp. energy diff.'
!      endif
CedN     write(13,5003) 'upper',eup
     

      write(13,*)
      write(*,*)' Value of nms parameter K_NMS (a.u.) for lower state'
      read(*,*) K_NMSlow 
      write(13,5000)'lower',K_NMSlow
      
      write(*,*) ' Value of nms parameter K_NMS (a.u.) for upper state'
      read(*,*) K_NMSup
      write(13,5000) 'upper',K_NMSup

      write(13,*)
      write(*,*) ' Value of sms parameter K_SMS (a.u.) for lower state'
      read(*,*) K_SMSlow
      write(13,5001) 'lower',K_SMSlow

      write(*,*) ' Value of sms parameter K_SMS (a.u.) for upper state'
      read(*,*) K_SMSup
      write(13,5001) 'upper',K_SMSup


      write(13,*)
      write(*,*)  
     1    ' Density (a_0^-3) at the nucleus for lower state'
      read(*,*) eldenslow
      write(13,5002) 'lower', eldenslow
      
      write(*,*)  
     1    ' Density (a_0^-3) at the nucleus for upper state'
      read(*,*) eldensup
      write(13,5002) 'upper', eldensup

      write(*,*) ' Type the number of isotopes considered'
      read(*,*) niso

!      write(13,*)
!      write(13,4006)
!      diffeh= (ediffcm)*cmM12eh!/eh2cmM1
!      write(13,4003) ' Trans.  ',diffeh,diffeh*eh2ev,ediffcm,
!     &   diffeh*eh2Mhz,ediffwave


      write(13,*)
      write(13,*) ' Number of isotopes',niso
      write(13,*)
      write(13,*) ' Isotope label, atomic isotope mass (u) and ',
     & 'rms nuclear charge radius (fm)'

      write(*,*) ' (If more than 1, please give them from the',
     &' lightest to the heaviest)'
      do 10 i = 1,niso
         write(*,*) ' Isotope label (character*5) (ex: 142Nd )'
         read(*,' (a5)') isolabel(i)
         IF (mostab == isolabel(i)) THEN 
            theone =i
         ENDIF
         write(*,*) ' Atomic isotope mass (u)'
         read(*,*) isomass(i)
         write(*,*)' Root mean square nuclear radius sqrt(<r**2 >) (fm)'
         read(*,*) rms(i)
         write(13,1000) isolabel(i),isomass(i),rms(i)
CcN     1   isolabel(i),'    ',isomass(i),'u','    ',rms(i),'fm'
10    continue 
       IF (niso>1) THEN
         IF (theone .eq.0) THEN
           write(*,*) ' You did not introduce the',
     &   ' reference isotope!'
          STOP
         ENDIF
         write(13,*)
         write(13,*) 'The reference isotope is ',mostab
         write(13,*)
         write(13,*)  ' The transition energy provided by the user is'
         write(13,4006)
         diffeh= (ediffcm)*cmM12eh!/eh2cmM1
         write(13,4003) '         ',diffeh,diffeh*eh2ev,ediffcm,
     &   diffeh*eh2Mhz,ediffwave
       ELSE
                  write(13,*)
         write(13,*) 'There is only one isotope: ',mostab
         write(13,*) 'The transition is ',transition
!         write(13,*)
         write(13,4006)
         diffeh= (ediffcm)*cmM12eh!/eh2cmM1
         write(13,4003) '         ',diffeh,diffeh*eh2ev,ediffcm,
     &   diffeh*eh2Mhz,ediffwave
       ENDIF

Cmrg  Z*emass is the electron mass we need to substract from the
Cmrg  atomic mass to get the nuclear mass. Z is the number of protons,
Cmrg  ie. the number of electrons for the neutral atom.
Cmrg  M_N (A,Z) = M_A (A,Z) - Z*m_e + B_e(Z)
CedN  where B_e = B=14.4381*Z^(2.39)+1.55468*10^(-6)*Z^(5.35) eV
CedN  Bm = (B*e)/c^2
CedN  Bmu = Bm/1.660538782*10^(-27)
CedN  The mass default:
      B_e =14.4381*(z**(2.39))+(1.55468*(1.0d-6))*(z**(5.35))
      Bm = (B_e*elec)/(clum**2.0d0) !The mass default in kg
      Bmu = Bm/ukg              !The mass default in u
*     Calculate the nuclear mass from atomic mass, that is subtract
*     the electron mass
      write(13,*)
      write(13,*) ' Nuclear mass(es) for the isotope(s)'
      do 20 i = 1,niso
         isomass(i) = isomass(i) - z*emass+Bmu
         isomassau(i) = isomass(i)/emass
         write(13,*) isolabel(i),'    ',isomass(i),'u'
20    continue
***********************************************************************
***********************************************************************

      write(13,*)
      write(13,*) '============================='
      write(13,*) ' Level IS on the lower state '
      write(13,*) '============================='
      write(13,*)
      write(13,*) 'The values below are the level isotope shifts',
     & ' relatively to an infinite mass point charge nucleus system'
      write(13,*)
      write(13,4000) 
      Do i = 1,niso
      write(13,4002) isolabel(i)
!      write(13,4003) ' E_0    = ',elow,elow*eh2ev,elow*eh2cmM1,
!     &   elow*eh2Mhz,1.d10/(2.d0*elow*ry)
*Normal Mass Shift
      E_NMS = K_NMSlow/isomassau(i)
      write(13,4001) ' NMS    = ',E_NMS,E_NMS*eh2ev,
     &   E_NMS*eh2cmM1,E_NMS*eh2Mhz!,(E_NMS*eh2cmM1*1.0d2)**(-1.0d00)
*Specific Mass Shift
      E_SMS = K_SMSlow/isomassau(i)
      write(13,4001) ' SMS    = ',E_SMS,E_SMS*eh2ev, 
     &   E_SMS*eh2cmM1,E_SMS*eh2Mhz!,(E_SMS*eh2cmM1*1.0d2)**(-1.0d00)

*Mass Shift
      E_MS = (K_NMSlow+K_SMSlow)/isomassau(i)
      write(13,4001) ' MS     = ',E_MS,E_MS*eh2ev, 
     &   E_MS*eh2cmM1,E_MS*eh2Mhz!,(E_SMS*eh2cmM1*1.0d2)**(-1.0d00)
      write(13,*)'+'
*Field shift
      E_FSlow=z*(2.0d0/3.0d0)*pi*eldenslow*(1.0d-15/a0)**2*rms(i)**2
      write(13,4001) ' FS     = ',E_FSlow,E_FSlow*eh2ev, 
     &E_FSlow*eh2cmM1, E_FSlow*eh2Mhz
!      Etot=elow+(K_NMSlow+K_SMSlow)/isomassau(i)+E_FSlow
      Etot=(K_NMSlow+K_SMSlow)/isomassau(i)+E_FSlow
      write(13,4005)
      write(13,4001) ' IS     = ',Etot,Etot*eh2ev, 
     &Etot*eh2cmM1, Etot*eh2Mhz
      write(13,*)
       ENDDO

      write(13,*)
      write(13,*) '============================='
      write(13,*) ' Level IS on the upper state '
      write(13,*) '============================='
      write(13,*)
      write(13,*) 'The values below are the level isotope shifts', 
     & ' relatively to an infinite mass point charge nucleus system'
      write(13,*)
      write(13,4000)
      Do i = 1,niso
      write(13,4002) isolabel(i)
!      write(13,4003) ' E_0    = ',eup,eup*eh2ev,eup*eh2cmM1,
!     &   eup*eh2Mhz,1.d10/(2.d0*eup*ry)
*Normal Mass Shift
      E_NMS = K_NMSup/isomassau(i)
      write(13,4001) ' NMS    = ',E_NMS,E_NMS*eh2ev,E_NMS*eh2cmM1,
     &   E_NMS*eh2Mhz!,(E_NMS*eh2cmM1*1.0d2)**(-1.0d00)!,Ktemp
*Specific Mass Shift
      E_SMS = K_SMSup/isomassau(i)
      write(13,4001) ' SMS    = ',E_SMS,E_SMS*eh2ev,E_SMS*eh2cmM1,
     &   E_SMS*eh2Mhz!,(E_SMS*eh2cmM1*1.0d2)**(-1.0d00)
*Mass Shift
      E_MS = (K_NMSup+K_SMSup)/isomassau(i)
      write(13,4001) ' MS     = ',E_MS,E_MS*eh2ev,
     &   E_MS*eh2cmM1,E_MS*eh2Mhz!,(E_SMS*eh2cmM1*1.0d2)**(-1.0d00)
      write(13,*)'+'
*Field shift
      E_FSup=z*(2.0d0/3.0d0)*pi*eldensup*(1.0d-15/a0)**2*rms(i)**2
      write(13,4001) ' FS     = ',E_FSup,E_FSup*eh2ev,
     &E_FSup*eh2cmM1, E_FSup*eh2Mhz

!      Etot=eup+(K_NMSup+K_SMSup)/isomassau(i)+E_FSup
      Etot=(K_NMSup+K_SMSup)/isomassau(i)+E_FSup
      write(13,4005)
      write(13,4001) ' IS     = ',Etot,Etot*eh2ev,
     &Etot*eh2cmM1, Etot*eh2Mhz
      write(13,*)
       ENDDO

      write(13,*)
      write(13,*) '==========================================='
      write(13,*) ' Transition isotope shifts for one isotope '
      write(13,*) '==========================================='
      write(13,*) 'The values below are the transition isotope shifts',
     & ' relatively to an infinite mass point charge nucleus system'

*LET US calculate all the differences
          DNMS= K_NMSup-K_NMSlow
          DSMS= K_SMSup-K_SMSlow
          DFS= eldensup-eldenslow
      write(13,*)
!      write(13,*) 'The isotope shifts in this transition are:'
!      write(13,*)
      write(13,4006)
!      write(13,4010) ' Trans.  ',diffeh,diffeh*eh2ev,ediffcm,
!     &   diffeh*eh2Mhz,ediffwave
!      write(13,*)
      DO i=1,niso
          write(13,4002)isolabel(i)
*Normal Mass Shift
          E_NMS=DNMS/isomassau(i)
          write(13,4003) ' NMS    = ',E_NMS,
     &    E_NMS*eh2ev,
     &    E_NMS*eh2cmM1,E_NMS*eh2Mhz,
     &    -((ediffwave)**2*E_NMS*eh2Mhz*(1.0d6*1.0d-10)/(clum))
*Specific Mass Shift 
          E_SMS=DSMS/isomassau(i)
          write(13,4003) ' SMS    = ',E_SMS,
     &    E_SMS*eh2ev,
     &    E_SMS*eh2cmM1,E_SMS*eh2Mhz,
     &    -((ediffwave)**2*E_SMS*eh2Mhz*(1.0d-4)/(clum))
*Mass Shift
          E_MS = (DNMS+DSMS)/isomassau(i)
          write(13,4003) ' MS     = ',E_MS,E_MS*eh2ev,
     &    E_MS*eh2cmM1,E_MS*eh2Mhz,
     &  -((ediffwave)**2*E_MS*eh2Mhz*(1.0d-4)/(clum))
          write(13,*)'+'
* Field shift
!          DE_FS  =z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*Dr2
!          DE_FSA =z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*rms(j)**2
          DE_FS  =z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*rms(i)**2
          write(13,4003) ' FS     = ',DE_FS,DE_FS*eh2ev,
     &    DE_FS*eh2cmM1, DE_FS*eh2Mhz,
     &    -((ediffwave)**2*DE_FS*eh2Mhz*(1.0d-4)/clum)
!      DEtot=diff(1)+(DNMS+DSMS)/isomassau(i)+DE_FS
          DEtot=(DNMS+DSMS)/isomassau(i)+DE_FS
          write(13,4009)
          write(13,4003) ' IS     = ',DEtot,DEtot*eh2ev,
     &  DEtot*eh2cmM1, DEtot*eh2Mhz,
     & (-(ediffwave**2*DEtot*eh2Mhz*(1.0d-4))/clum)
          write(13,*)
      ENDDO
**************************
**************************
      IF (niso .GE. 2) THEN
      write(13,*)
      write(13,*) '================================================'
!     &'==================='
      write(13,*) ' Transition energies for the different isotopes '
!     &' reference isotope '
      write(13,*) '================================================'
!     &'==================='
      write(13,*)

      write(13,*) ' The reference isotope is ',mostab
      write(13,4006)
      write(13,4003) '         ',diffeh,diffeh*eh2ev,ediffcm,
     &   diffeh*eh2Mhz,ediffwave
      write(13,*)
!      write(13,4013) transition
      mostabms=(DNMS+DSMS)/isomassau(theone)
      mostabfs=z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*rms(theone)**2
      DO i=1,niso
        IF (i .ne. theone) THEN
          DE_FS  =z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*rms(i)**2
          E_NMS=DNMS/isomassau(i)
          E_SMS=DSMS/isomassau(i)
          write(13,4002)isolabel(i)
          write(13,4006)
          DEtot=diffeh-(mostabms+mostabfs)+E_NMS+E_SMS+DE_FS
          write(13,4003) '          ',DEtot,DEtot*eh2ev,
     &  DEtot*eh2cmM1, DEtot*eh2Mhz,(-(ediffwave**2*((E_NMS+E_SMS)
     & +DE_FS-(mostabms+mostabfs))*eh2Mhz*(1.0d-4))/clum)+ediffwave
          write(13,*)
        ENDIF
      ENDDO

!      IF (niso .GE. 2) THEN
      write(13,*)
      write(13,*) '==============================================='
      write(13,*) ' Transition isotope shifts for an isotope pair '
      write(13,*) '==============================================='
      write(13,*)

!      write(13,4006)
!      write(13,4003) ' Trans.  ',diffeh,diffeh*eh2ev,ediffcm,
!     &   diffeh*eh2Mhz,ediffwave
!      write(13,*)
      write(13,4012) transition
      write(13,*)
      write(13,4006)
      DO i=1,niso
        DO j=i+1, niso
*LET US calculate all the differences
          write(13,4007) isolabel(j),isolabel(i)
!          DNMS= K_NMSup-K_NMSlow
!          DSMS= K_SMSup-K_SMSlow
!          DFS= eldensup-eldenslow
          DnuclM= (1/isomassau(j)-1/isomassau(i))**(-1)
          Dr2=rms(j)**2-rms(i)**2
          write(*,4004) isolabel(j),isolabel(i),rms(j),rms(i),Dr2
          write(*,*) 'Do you want to change this value?'
          READ(*,*) ans
          IF ((ans.eq.'y').or.(ans.eq.'Y')) THEN
            write(*,*) 'Introduce the value in fm'
            READ(*,*) Dr2             
          ENDIF
******
*Normal Mass Shift
          write(13,4003) ' NMS    = ',DNMS/DnuclM,
     &    DNMS*eh2ev/DnuclM,
     &    DNMS*eh2cmM1/DnuclM,DNMS*eh2Mhz/DnuclM,
     &    -((ediffwave)**2*DNMS*eh2Mhz*(1.0d-4)/(DnuclM*clum))
*Specific Mass Shift
          write(13,4003) ' SMS    = ',DSMS/DnuclM,
     &    DSMS*eh2ev/DnuclM,
     &    DSMS*eh2cmM1/DnuclM,DSMS*eh2Mhz/DnuclM,
     &    -((ediffwave)**2*DSMS*eh2Mhz*(1.0d-4)/(DnuclM*clum))
*Mass Shift
      E_MS = (DNMS+DSMS)/DnuclM
      write(13,4003) ' MS     = ',E_MS,E_MS*eh2ev,
     &   E_MS*eh2cmM1,E_MS*eh2Mhz,
     & -((ediffwave)**2*E_MS*eh2Mhz*(1.0d-4)/(clum))
      write(13,*)'+'

* Field shift
          DE_FS  =z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*Dr2
          DE_FSA =z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*rms(j)**2
          DE_FSAp=z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*rms(i)**2
          write(13,4003) ' FS     = ',DE_FS,DE_FS*eh2ev,
     &    DE_FS*eh2cmM1, DE_FS*eh2Mhz,
     &    -((ediffwave)**2*DE_FS*eh2Mhz*(1.0d-4)/clum)

!      DEtot=diff(1)+(DNMS+DSMS)/isomassau(i)+DE_FS
      DEtot=(DNMS+DSMS)/DnuclM+DE_FS
      write(13,4009)
      write(13,4003) ' IS     = ',DEtot,DEtot*eh2ev,
     & DEtot*eh2cmM1, DEtot*eh2Mhz,
     & (-(ediffwave**2*DEtot*eh2Mhz*(1.0d-4))/clum)
      write(13,*)
        ENDDO
      ENDDO
      ENDIF

1000  FORMAT (A5,14X, F15.10,' u',14X,ES11.4,' fm')
1001  FORMAT (ES11.4,'(a_0^-3)')
4000  FORMAT (16X,'E_h',16X,' eV',13X,' cm-1',16X,'MHz')
4001  FORMAT (A9,1X,ES17.10,F17.6,2X,ES17.10,2X,ES17.10)
4002  FORMAT ('For ',A5,':')
4003  FORMAT (A9,1X,ES17.10,2X,F17.6,2X,ES17.10,2X,ES17.10,2X,ES17.10)
4004  FORMAT ('The diff. of mean-square charge radii ',A5,'-',A5,/,
     & 'is:', ES11.4,'^2 -',ES11.4,'^2 = ',ES11.4)
4005  FORMAT (' --------------------------------------------------------
     &--------------------------')
4006  FORMAT (15X,'E_h',19X,' eV',13X,' cm-1',15X,'MHz',14X,'Angstrom')
4007  FORMAT (A5,'-',A5)
4008  FORMAT (' The transition value is ',ES17.6,A10)
4009  FORMAT (' --------------------------------------------------------
     &---------------------------------------------')
4010  FORMAT (A9,1X,ES15.7,4X,ES15.7,4X,ES15.7,4X,ES15.7,4X,ES15.7)
4011  FORMAT ('For the isotope ',A5,' the transition energy is:')
4012  FORMAT (' For the transition ',A72,/,' the isotope ',
     &'shifts between pairs of isotopes are:')
4013  FORMAT (' For the transition ',A72,/,' the line of other ',
     &'isotopes are ',/)
5000  FORMAT (' Value of nms parameter K_NMS for ', A5,' state ',
     &  ES17.10,' E_h m_e')
5001  FORMAT (' Value of sms parameter K_SMS for ', A5,' state ',
     &  ES17.10,' E_h m_e')
5002  FORMAT (' Density at the nucleus for ', A5, ' state ',
     &  ES17.10,' a_0^-3')
5003  FORMAT (' Energy infinite nucl. mass for ', A5, ' state ', 
     &  ES17.10,' E_h')
      end
