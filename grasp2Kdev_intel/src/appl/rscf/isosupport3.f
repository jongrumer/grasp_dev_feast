*
************************************************************************
*
      program isosupport
* 
* written by Per Jonsson and Michel Godefroid
* version  October 10, 2012
*
*                                                                      *
************************************************************************

      implicit double precision(a-h,o-z)

      parameter (emass = 5.4857990946d-4, ry = 10973731.568539d0)
      parameter (a0 = 0.52917721092d-10, pi = 3.14159265359d0)
      doubleprecision isomass(10),el(10),eu(10),diff(10),eln(10),eun(10)
     :,rme(10),muovm(10),rydm(10),isomassau(10)
      doubleprecision wlength(10),dwlength(10),rms(10)
      doubleprecision elns(10),euns(10)
      DOUBLE PRECISION K_NMSlow, K_NMSup, K_SMSlow, K_SMSup,E_FSlow
      character*5 isolabel(10)
      character*72 transition
      CHARACTER*1 ans

      open(unit=13,file='ISTrans',status='unknown')
      open(unit=14,file='isodata',status='unknown')

      convmkmhz = 29.9792458000000
      convmhzau = 1.5198298460045*1.d-10
      eh2ev  = 27.21138505
      eh2Mhz  = 6.579683920729*1.d9
      eh2cmM1= 2.194746313708*1.d9
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
      write(*,*) 
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
      write(*,*) 
      write(*,*) 
      write(*,*) ' Specify the transition'
      read(*,'(a72)') transition

      write(13,*)
      write(13,*) ' Transition'
      write(13,'(a72)') transition

      write(*,*) ' Nuclear charge'
      read(*,*) z

      write(13,*)
      write(13,*) ' Nuclear charge',z

      write(13,*)
      write(*,*) ' Energy (a.u.) infinite nucl. mass for lower state'
      read(*,*) elow
      write(13,5003) 'lower',elow

991   write(*,*) 
     1' Energy (a.u.) infinite nucl. mass for upper state (1) or '
      write(*,*) 
     2' energy diff. (cm^-1) between upper and lower state (2)'
      read(*,*) n
      if (n.eq.1) then
        write(*,*) ' Energy (a.u.) infinite nucl. mass for upper state'
        read(*,*) eup
      elseif (n.eq.2) then
        write(*,*) ' Energy diff. (cm^-1) between upper and lower state'
        read(*,*) ediff
        eup = elow + 50.d0*ediff/ry
      else
         goto 991
      endif
      if (n.eq.1) then
         write(13,*) ' Energy for upper state from calculation'
      else
         write(13,*) ' Energy for upper state from exp. energy diff.'
      endif
      write(13,5003) 'upper',eup

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
     1    ' Mod. density (a_0^-3) at the nucleus for lower state'
      read(*,*) eldenslow
      write(13,5002) 'lower', eldenslow
      
      write(*,*)  
     1    ' Mod. density (a_0^-3) at the nucleus for upper state'
      read(*,*) eldensup
      write(13,5002) 'upper', eldensup

      write(*,*) ' Number of isotopes'
      write(*,*) ' If more than 1, please give them from the',
     &' lightest to the heaviest'
      read(*,*) niso 
      write(13,*)
      write(13,*) ' Number of isotopes',niso
      
CedN      if (niso. gt. 2) then
CedN        write(*,*)
CedN     1     'No more than one pair for a given f(Z)_MM'
CedN        stop
CedN      end if

      write(13,*)
      write(13,*) ' Isotope label, atomic isotope mass (u) and root mean square
     1 nuclear radius (fm)'
      do 10 i = 1,niso
         write(*,*) ' Isotope label (character*5)'
         read(*,'(a5)') isolabel(i)
         write(*,*) ' Atomic isotope mass (u)'
         read(*,*) isomass(i)
         write(*,*) ' Root meansquare nuclear radius sqrt(<r**2 >) (fm)'
         read(*,*) rms(i)
         write(13,1000) isolabel(i),isomass(i),rms(i)
CcN     1   isolabel(i),'    ',isomass(i),'u','    ',rms(i),'fm'
10    continue 

Cmrg  write(*,*) ' Number of electrons' 
Cmrg  read(*,*) ne
Cmrg  write(13,*)
Cmrg  write(13,*) ' Number of electrons',ne 
Cmrg  Z*emass is the electron mass we need to substract from the
Cmrg  atomic mass to get the nuclear mass. Z is the number of protons,
Cmrg  ie. the number of electrons for the neutral atom.
Cmrg  M_N (A,Z) = M_A (A,Z) - Z*m_e + B_e(Z)
CedN  where B_e = B=14.4381*Z^(2.39)+1.55468*10^(-6)*Z^(5.35) eV
CedN  Bm = (B*e)/c^2
CedN  Bmu = Bm/1.660538782*10^(-27)
CedN      ne = ifix(z)
CedN  The mass default:
      B_e =14.4381*(z**(2.39))+(1.55468*10**(-6))*(z**(5.35))
      Bm = (B_e*elec)/(clum**2) !The mass default in kg
      Bmu = Bm/ukg              !The mass default in u

*     Calculate the nuclear mass from atomic mass, that is subtract
*     the electron mass
      write(13,*)
      write(13,*) ' Nuclear mass for the isotopes'
      do 20 i = 1,niso
         isomass(i) = isomass(i) - z*emass+Bmu
         isomassau(i) = isomass(i)/emass
         write(13,*) isolabel(i),'    ',isomass(i),'u'
*         muovm(i) = isomass(i)/(isomass(i) + emass)
*         write(13,*) '     mu/m ratio = ',muovm(i)
*         rydm(i) = ry * muovm(i)
*         write(13,*) '     Rydberg_M  = ',rydm(i)*0.01D0,' cm-1 '
*         write(13,*) ' (2*)Rydberg_M  = ',
*     :                                rydm(i)*0.02D0,' cm-1 '
*         write(13,*)
*         rme(i) = emass * muovm(i)
20    continue

CCC      write(13,*)
CCC      write(13,*) ' Nuclear mass for the isotopes'
CCC      do 20 i = 1,niso
CCC         isomass(i) = isomass(i) - ne*emass
CCC         write(13,*) isolabel(i),'    ',isomass(i),'u'
CCC         muovm(i) = isomass(i)/(isomass(i) + emass)
CCC         write(13,*) '     mu/m ratio = ',muovm(i)
CCC         rydm(i) = ry * muovm(i)
CCC         write(13,*) '     Rydberg_M  = ',rydm(i)*0.01D0,' cm-1 '
CCC         write(13,*) ' (2*)Rydberg_M  = ',
CCC     :                                rydm(i)*0.02D0,' cm-1 '
CCC         write(13,*)
CCC         rme(i) = emass * muovm(i)
CCC20    continue

      write(13,*)
      write(13,*) '============='
      write(13,*) ' Lower state'
      write(13,*) '============='
      write(13,4000) 
*2321  FORMAT (1X,A12,'= ',ES17.10,' au',F17.5,' eV',2X,ES22.15,' cm-1')
      Do i = 1,niso
      write(13,4002) isolabel(i)
      write(13,4001) ' E_0    = ',elow,elow*eh2ev,elow*eh2cmM1,
     &   elow*eh2Mhz
      write(13,4001) ' NMS    = ',K_NMSlow/isomassau(i),
     &   K_NMSlow*eh2ev/isomassau(i),
     &   K_NMSlow*eh2cmM1/isomassau(i),K_NMSlow*eh2Mhz/isomassau(i)
      write(13,4001) ' SMS    = ',K_SMSlow/isomassau(i),
     &   K_SMSlow*eh2ev/isomassau(i), 
     &   K_SMSlow*eh2cmM1/isomassau(i),K_SMSlow*eh2Mhz/isomassau(i)
      E_FSlow=z*(2.0d0/3.0d0)*pi*eldenslow*(1.0d-15/a0)**2*rms(i)**2
      write(13,4001) ' FS     = ',E_FSlow,E_FSlow*eh2ev, 
     &E_FSlow*eh2cmM1, E_FSlow*eh2Mhz
      Etot=elow+(K_NMSlow+K_SMSlow)/isomassau(i)+E_FSlow
      write(13,4005)
      write(13,4001) ' SUM    = ',Etot,Etot*eh2ev, 
     &Etot*eh2cmM1, Etot*eh2Mhz
      write(13,*)
       ENDDO

CedN      write(13,*)
CedN      write(13,*) ' Energy for infinite nuclear mass'

CedN      write(13,*) '        ',elow,'E_h'

*     Lower energy for the different isotopes

*      write(13,*)
*      write(13,*) ' Energy corrected for normal mass effect'

*      do 25 i = 1,niso
Cmrg     el(i) = elow*isomass(i)/(isomass(i) + emass) 
*         el(i) = elow*muovm(i)
*         eln(i) = el(i)
*         write(13,*) isolabel(i),'    ', el(i),'E_h'
*25    continue

*      write(13,*)
CedN      write(13,*)' Energy corrected for normal and specific mass effect'

CedN      do 30 i = 1,niso
Cmrg     el(i) = el(i) + emass*gradlow/isomass(i)
CedN         el(i) = rme(i)*muovm(i)*(K_NMSlow+K_SMSlow)/isomass(i)
CedN         elns(i) = el(i)
CedN         write(13,*) isolabel(i),'    ', el(i),'E_h'
CedN30    continue        

CedN      write(13,*)
CedN      write(13,*)' Energy corrected for normal and specific mass effect'
CedN     1          ,' and fieldshift'

CedN      do 31 i = 1,niso
CedN        el(i) = el(i) + eldenslow*fz*rms(i)*rms(i)*convmhzau/(4.d0*z)
Cmrg     rms(i) = rms(i)*1.d-15/a0
Cmrg     el(i) = el(i) + 2.d0*pi*z*eldenslow*rms(i)*rms(i)/3.d0
CedN         write(13,*) isolabel(i),'    ', el(i),'E_h'
CedN31    continue        

      write(13,*)
      write(13,*) '============='
      write(13,*) ' Upper state'
      write(13,*) '============='
      write(13,4000)
      Do i = 1,niso
      write(13,4002) isolabel(i)
      write(13,4001) ' E_0    = ',eup,eup*eh2ev,eup*eh2cmM1,
     &   eup*eh2Mhz
      write(13,4001) ' NMS    = ',K_NMSup/isomassau(i),
     &   K_NMSup*eh2ev/isomassau(i),
     &   K_NMSup*eh2cmM1/isomassau(i),K_NMSup*eh2Mhz/isomassau(i)
      write(13,4001) ' SMS    = ',K_SMSup/isomassau(i),
     &   K_SMSup*eh2ev/isomassau(i),
     &   K_SMSup*eh2cmM1/isomassau(i),K_SMSup*eh2Mhz/isomassau(i)
      E_FSup=z*(2.0d0/3.0d0)*pi*eldensup*(1.0d-15/a0)**2*rms(i)**2
      write(13,4001) ' FS     = ',E_FSup,E_FSup*eh2ev,
     &E_FSup*eh2cmM1, E_FSup*eh2Mhz
      Etot=eup+(K_NMSup+K_SMSup)/isomassau(i)+E_FSup
      write(13,4005)
      write(13,4001) ' SUM    = ',Etot,Etot*eh2ev,
     &Etot*eh2cmM1, Etot*eh2Mhz
      write(13,*)
       ENDDO

CedN      write(13,*)
CedN      write(13,*) ' Energy for infinite nuclear mass'

CedN      write(13,*) '        ',eup,'E_h'

*     Upper energy for the different isotopes

CedN      write(13,*)
CedN      write(13,*) ' Energy corrected for normal mass effect'

CedN      do 37 i = 1,niso
Cmrg     eu(i) = eup*isomass(i)/(isomass(i) + emass) 
CedN         eu(i) = eup*muovm(i)
CedN         eun(i) = eu(i)
CedN         write(13,*) isolabel(i),'    ', eu(i),'E_h'
CedN37    continue        

CedN      write(13,*)
CedN      write(13,*)' Energy corrected for normal and specific mass effect'

CedN      do 40 i = 1,niso
Cmrg     eu(i) = eu(i) + emass*gradup/isomass(i)
CedN         eu(i) = eu(i) + rme(i)*muovm(i)*gradup/isomass(i)
CedN         euns(i) = eu(i)
CedN         write(13,*) isolabel(i),'    ', eu(i),'E_h'
CedN40    continue        

CedN      write(13,*)
CedN      write(13,*)' Energy corrected for normal and specific mass effect'
CedN     1          ,' and fieldshift'

CedN      do 41 i = 1,niso
CedN         eu(i) = eu(i) + eldensup*fz*rms(i)*rms(i)*convmhzau/(4.d0*z)
Cmrg     eu(i) = eu(i) + 2.d0*pi*z*eldensup*rms(i)*rms(i)/3.d0 
CedN         write(13,*) isolabel(i),'    ', eu(i),'E_h'
CedN41    continue        

      write(13,*)
      write(13,*) '===================='
      write(13,*) ' Energy differences'
      write(13,*) '===================='

      write(13,*)
      write(13,*) ' Energy difference between upper and lower state'
      write(13,*)
CedN      write(13,*) ' Difference for infinite nuclear mass:'
CedN      write(13,*) 
     
*     Energy difference

      diff(1) = eup - elow
      write(13,4006)
CedN      write(13,*) '        ',diff(1),'E_h'
CedN      write(13,*) '        ',diff(1)*ry*2.d0*0.01d0,'cm-1'

CedN      write(13,*) ' Transition wavelength (Angstrom)'
CedN      wlength(1) = 1.d10/(2.d0*diff(1)*ry)
CedN      write(13,*) '        ', wlength(1),'A'

      write(13,*)
*      write(13,*) ' Differences with normal mass correction'
      DO i=1,niso
        DO j=i+1, niso
*LET US calculate all the differences
          write(13,4007) isolabel(j),isolabel(i)
          write(13,4003) ' DiffE_0= ',diff(1),diff(1)*eh2ev, 
     &      diff(1)*eh2cmM1,diff(1)*eh2Mhz,1.d10/(2.d0*diff(1)*ry)
          DNMS= K_NMSup-K_NMSlow
          DSMS= K_SMSup-K_SMSlow
          DFS= eldensup-eldenslow
          DnuclM= (1/isomassau(j)-1/isomassau(i))**-1
          Dr2=rms(j)**2-rms(i)**2
          write(*,4004) isolabel(j),isolabel(i),rms(j)**2,rms(i)**2,Dr2
          write(*,*) 'Do you want to change it?'
          READ(*,*) ans
          IF ((ans.eq.'y').or.(ans.eq.'Y')) THEN
            write(*,*) 'Introduce the value in fm'
          ENDIF
******
!      write(13,4001) ' E_0    = ',eup,eup*eh2ev,eup*eh2cmM1,
!     &   eup*eh2Mhz
          write(13,4001) ' NMS    = ',DNMS/DnuclM,
     &    DNMS*eh2ev/DnuclM,
     &    DNMS*eh2cmM1/DnuclM,DNMS*eh2Mhz/DnuclM
CN
          write(13,4001) ' SMS    = ',DSMS/DnuclM,
     &    DSMS*eh2ev/DnuclM,
     &    DSMS*eh2cmM1/DnuclM,DSMS*eh2Mhz/DnuclM
CN
          DE_FS=z*(2.0d0/3.0d0)*pi*DFS*(1.0d-15/a0)**2*Dr2
          write(13,4001) ' FS     = ',DE_FS,DE_FS*eh2ev,
     &    DE_FS*eh2cmM1, DE_FS*eh2Mhz
      DEtot=diff(1)+(DNMS+DSMS)/isomassau(i)+DE_FS
      write(13,4005)
      write(13,4001) ' SUM    = ',DEtot,DEtot*eh2ev,
     &DEtot*eh2cmM1, DEtot*eh2Mhz
      write(13,*)
        ENDDO
      ENDDO
*     Energy difference

*      do 50 i = 1,niso
*         diff(i) = eun(i) - eln(i)
*         write(13,*) isolabel(i),'    ',diff(i),'E_h'
*50    continue
*
*      write(13,*)
*      do 55 i = 1,niso
*         write(13,*) isolabel(i),'    ',diff(i)*ry*2.d0*0.01d0,'cm-1'
*55    continue
*
*      write(13,*)
*      write(13,*) isolabel(1),'-',isolabel(niso),'    ',
*     1   diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0,'cm-1'
*      dif1 = diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0
      

*      write(13,*)
*      write(13,*) ' Transition wavelength (Angstrom)'

*      do 60 i = 1,niso
*         wlength(i) = 1.d10/(2.d0*diff(i)*ry)
*         write(13,*) isolabel(i),'    ', wlength(i),'A'
*60    continue

*      write(13,*)
*      write(13,*) ' Wavelength difference between the isotopes'

*      do 70 i = 2,niso
*         dwlength(i) = dabs(wlength(1) - wlength(i))
*         write(13,*) isolabel(1),isolabel(i),dwlength(i),'   A'
*70    continue

*      write(13,*)
*      write(13,*) 
*     :' Differences with normal and specific mass correction'
*
*     Energy difference

*      do 150 i = 1,niso
*         diff(i) = euns(i) - elns(i)
*         write(13,*) isolabel(i),'    ',diff(i),'E_h'
*150    continue
*
*      write(13,*)
*      do 155 i = 1,niso
*         write(13,*) isolabel(i),'    ',diff(i)*ry*2.d0*0.01d0,'cm-1'
*155    continue
*
*      write(13,*)
*      write(13,*) isolabel(1),'-',isolabel(niso),'    ',
*     1   diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0,'cm-1'
*      dif2 = diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0

*      dif21 = dif2 - dif1
*      write(13,*)
*      write(13,*) ' Transition SMS (cm-1) = ', dif21
*      dif21 = dif21 * 2.99792458D+04
*      write(13,*) ' Transition SMS (MHz)  = ', dif21
*
*
*      write(13,*)
*      write(13,*) ' Transition wavelength (Angstrom)'
*
*      do 160 i = 1,niso
*         wlength(i) = 1.d10/(2.d0*diff(i)*ry)
*         write(13,*) isolabel(i),'    ', wlength(i),'A'
*160    continue

*      write(13,*)
*      write(13,*) ' Wavelength difference between the isotopes'
*
*      do 170 i = 2,niso
*         dwlength(i) = dabs(wlength(1) - wlength(i))
*         write(13,*) isolabel(1),isolabel(i),dwlength(i),'   A'
*170    continue
*
*      write(13,*)
*      write(13,*)
*      write(13,*) ' Differences with normal and specific mass',
*     1            ' and field shift correction'
*
**     Energy difference
*
*      do 180 i = 1,niso
*         diff(i) = eu(i) - el(i)
*         write(13,*) isolabel(i),'    ',diff(i),'E_h'
*180    continue
*
*      write(13,*)
*      do 185 i = 1,niso
*         write(13,*) isolabel(i),'    ',diff(i)*ry*2.d0*0.01d0,'cm-1'
*185    continue
*
*      write(13,*)
*      write(13,*) isolabel(1),'-',isolabel(niso),'    ',
*     1   diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0,'cm-1'
*      write(13,*)
*      dif3 = diff(1)*ry*2.d0*0.01d0-diff(niso)*ry*2.d0*0.01d0
*
*      dif32 = dif3 - dif2
*      write(13,*)
*      write(13,*) ' Transition FS  (cm-1) = ', dif32
*      dif32 = dif32 * 2.99792458D+04
*      write(13,*) ' Transition FS  (MHz)  = ', dif32
*
*      write(13,*) ' Transition wavelength (Angstrom)'
*
*      do 190 i = 1,niso
*         wlength(i) = 1.d10/(2.d0*diff(i)*ry)
*         write(13,*) isolabel(i),'    ', wlength(i),'A'
*190    continue
*
*      write(13,*)
*      write(13,*) ' Wavelength difference between the isotopes'
*
*      do 200 i = 2,niso
*         dwlength(i) = dabs(wlength(1) - wlength(i))
*         write(13,*) isolabel(1),isolabel(i),dwlength(i),'   A'
*200    continue

1000  FORMAT (A5,4X, F15.10,'u',4X,ES11.4,' Fm')
1001  FORMAT (ES11.4,'(a_0^-3)')
4000  FORMAT (15X,'E_h',17X,' eV',12X,' cm-1',15X,'MHz')
4001  FORMAT (A9,1X,ES17.10,F17.5,2X,ES17.10,2X,ES17.10)
4002  FORMAT ('For ',A5,':')
4003  FORMAT (A9,1X,ES17.10,F17.5,2X,ES17.10,2X,ES17.10,ES17.10)
4004  FORMAT ('The diff. of mean-square charge radii ',A5,'-',A5,/,
     & 'is:', ES11.4,'^2-',ES11.4,'^2= ',ES11.4)
4005  FORMAT ('---------------------------------------------------------
     &---------------------------')
4006  FORMAT (15X,'E_h',17X,' eV',12X,' cm-1',17X,'MHz',12X,'Angstrom')
4007  FORMAT (A5,'-',A5)
5000  FORMAT (' Value of nms parameter K_NMS for ', A5,' state ',
     &  ES17.10,' E_h m_e')
5001  FORMAT (' Value of sms parameter K_SMS for ', A5,' state ',
     &  ES17.10,' E_h m_e')
5002  FORMAT (' Mod. density at the nucleus for ', A5, ' state ',
     &  ES17.10,' a_0^-3')
5003  FORMAT (' Energy infinite nucl. mass for ', A5, ' state ', 
     &  ES17.10,' E_h')
      end
