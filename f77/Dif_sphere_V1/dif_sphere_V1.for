
      PROGRAM Diff_sphere
c KWH this version edited May 2010 to include a more robust Bessel function evaluation
c
c

c this program adapted to model diffusion in a sphere
c with additional delta function
c KWH 3/12/2007
c this program and its associated subroutines
c have been extensively modified by KWH 
c 2/26/96 for fitting multiple spectra
c at their associated Qs simultaneously
c also modified 4/22/97 to include a dector
c normalization routine
c modified 9/27/99 to include more arrays for new QENS

      implicit double precision(a-h,o-z)

      EXTERNAL SQWIN3,FUNCT3
      COMMON/TITLES/PNAMES,TX(4),TY(4)
      COMMON/PARTIT/NPARAS
      COMMON/CALDAT/CALPAR(100),YCALC(5000),WT(5000)
      COMMON/PARCTL/NumL,IBG,IWT,INST
      COMMON/FILES/RFILE,LENRF,SFILE,LENSF,SPNUM
      COMMON/PLTINF/SPEED,RTEMP,STEMP
      common/tofunc/MomTra, MolRad
c declaration of new arrays 9/21/99
      integer totspc,fitspc,ispec(25),npspec,nsetsp(25),nlimsp(25,6),
     +        nfitsp(25),kxmin(25),kxmax(25),npltsp(25),lxminsp(25),
     +        lxmaxsp(25),nepsp(25),nrptot(25)
      real*8  xlimsp(25,6),xmin(25),xmax(25),ymin(25),ymax(25),
     +        xscale(25,2),yscale(25,2),sysp(25,5000),sxsp(25,5000),
     +        syersp(25,5000),wtsp(25,5000),fiwtsp(25,5000),
     +        ycalcsp(25,5000),
     +        parsp(500),stepsp(500),esumsp(25),rysp(25,5000),
     +        ryersp(25,5000),qspec(25),
     +        rxminsp(25),rxmaxsp(25),wxsp(25,5000),risp(25,5000)
c new variables
	real*8 xnl(100)
	integer nind(100),lind(100)
c the global common for these arrays has been established as below:
      COMMON/TOTDAT/ totspc,fitspc,ispec,npspec,nsetsp,xlimsp,nlimsp,
     +               nfitsp,
     +               xmin,xmax,ymin,ymax,xscale,yscale,kxmin,kxmax,
     +               npltsp,
     +               sysp,sxsp,syersp,wtsp,fiwtsp,ycalcsp,parsp,stepsp,
     +               lxminsp,lxmaxsp,nepsp,esumsp,rysp,ryersp,rxminsp,
     +               rxmaxsp,nrptot,wxsp,risp,qspec

	common/xcoef/ xnl,nind,lind

      CHARACTER*5 TX,TY,ANS*1
      CHARACTER*15 PNAMES(40)
      CHARACTER*20 RFILE,SFILE
      real*8 MomTra, MolRad
      INTEGER SPNUM,SPEED


C Define Fortran TTY channel number and open.
      IIN=5
      IO = 6
20    WRITE(IO,40)
40    FORMAT(/,1X,'FIT_TRL, Fits QENS data using model for',
     +       'isotropic diffusion in a sphere',
     +       /1x,'plus an extra delta function')
      Write (io, 60)
60    format(/,1x,'Last update March 15, 2007')
C Define titles for x and y.
      TX(1) = 'Energ'
      TX(2) = 'y meV'
      TX(3) = '  (me'
      TX(4) = 'V)   '
      TY(1) = 'S(Q,w'
      TY(2) = ')    '
c open graphics window
	call pgbegin(0, '/gw', 1, 1)
	call pgpap(0.0,1.25)
80	continue
c open the data file for the xnl and read in the coefficients
	open(unit=1,file='xnl_file.txt',status='old')
	do 35 i=1,98
		read(1,*)lind(i),nind(i),xnl(i)
35	continue
	close(1)
140   continue


	ibg=2
220   WRITE(IO,240)
240   FORMAT(1X,'Weight using esd''s ? (<CR> = yes)',$)
      READ(IIN,180) ANS
180   format(a)
      IWT = 0
      IF (ANS.NE.'N'.AND.ANS.NE.'n') IWT = 1
C     Use unit weight if IWT is 0 or use esd if IW is 1.
C     Set array WT to 1.0 - this is useful if IWT is 0.
      DO 260 I=1,5000
260   WT(I) = 1.0e-5
C     Now calculate number of parameters......
      write(6,261)
261   format(1x,'enter number of spectra to be fit:',$)
      read(5,*)fitspc
	write(io,*) "DEBUG: fitspc=",fitspc
      npspec=7+ibg
      NPARAS = npspec*fitspc

C Variables for the fit are an overall scale factor, and the
C translational and rotational diffusion constants Diffusion constant
c is in 10^-5 cm^2/sec if energy is meV
C An extra Delta function intensity is needed for background
C contributions etc
      PNAMES(1)  = 'Scale factor   '
      PNAMES(2)  = '<u^2>^0.5 sph  '
      PNAMES(3)  = 'elastic pos    '
      PNAMES(4)  = 'Extra delta Int'
      PNAMES(5)  = '<u^2>^0.5 elast'
      PNAMES(6)  = 'Dsubt (10^-5)  '
      PNAMES(7)  = 'Sphere Radius  '
      pnames(8)  = 'Background Cons'
      PNAMES(9)  = 'Background Line'

C     Call main fitting routine.
360   CALL FITFUN(SQWIN3,FUNCT3)
380   WRITE(IO,400)
400   FORMAT(1X,'Continue in FITCON ? (<CR> = no)',$)
      READ(IIN,180) ANS
      IF (ANS.EQ.' '.OR.ANS.EQ.'N'.OR.ANS.EQ.'n') GOTO 420
C     Clear screen.
      CALL FTERAS
      GOTO 80
C     Reset terminal to text mode.
420   CALL FTRSET
c close graphics window
	call pgend
      END



c SQWIN3 has now been modified to read in all the spectra and stuff them into
c associated arrays, not very fancy or sophisticated

      SUBROUTINE SQWIN3
      implicit double precision(a-h,o-z)
C **********************************************
C Data must be output from the QENS QUELL program
c ***********************************************
      COMMON/SPLDAT/NSPNT,SX(5000),SY(5000),SYER(5000)
      COMMON/RESDAT/RY(5000),RYER(5000),WX(5000),RI(5000)
      COMMON/RESLIM/RTOF,NEP,LXMIN,LXMAX,ERRSUM
      COMMON/PARCTL/NumL,IBG,IWT,INST
      COMMON/FILES/RFILE,LENRF,SFILE,LENSF,NEXT
      COMMON/PLTINF/SPEED,RTEMP,STEMP
      common/tofunc/MomTra,MolRad
c declaration of new arrays 9/21/99
      integer totspc,fitspc,ispec(25),npspec,nsetsp(25),nlimsp(25,6),
     +        nfitsp(25),kxmin(25),kxmax(25),npltsp(25),lxminsp(25),
     +        lxmaxsp(25),nepsp(25),nrptot(25)
      real*8  xlimsp(25,6),xmin(25),xmax(25),ymin(25),ymax(25),
     +        xscale(25,2),yscale(25,2),sysp(25,5000),sxsp(25,5000),
     +        syersp(25,5000),wtsp(25,5000),fiwtsp(25,5000),
     +        ycalcsp(25,5000),
     +        parsp(500),stepsp(500),esumsp(25),rysp(25,5000),
     +        ryersp(25,5000),qspec(25),
     +        rxminsp(25),rxmaxsp(25),wxsp(25,5000),risp(25,5000)
c the global common for these arrays has been established as below:
      COMMON/TOTDAT/ totspc,fitspc,ispec,npspec,nsetsp,xlimsp,nlimsp,
     +               nfitsp,
     +               xmin,xmax,ymin,ymax,xscale,yscale,kxmin,kxmax,
     +               npltsp,
     +               sysp,sxsp,syersp,wtsp,fiwtsp,ycalcsp,parsp,stepsp,
     +               lxminsp,lxmaxsp,nepsp,esumsp,rysp,ryersp,rxminsp,
     +               rxmaxsp,nrptot,wxsp,risp,qspec

      CHARACTER*1 ANS,LINE*80
      CHARACTER*20 RFILE,SFILE,FILE
	character*80 cdum
      REAL*8 TEMPR,TEMPS
      real*8 twoth(25)
      real*8 MomTra, MolRad
      real temp,xtmp(5000),ytmp(5000),dytmp(5000),efobin,thobin,rx(5000)
      INTEGER SPEED

c*******************************************************************************
c
c KWH 1/2010 - this new modification reads in a grouped ASCII file of DAVE format
c
c*******************************************************************************

C Initialize variables.
      DATA NRF,RFILE(1:4)/0,'exit'/
      DATA LENFIL,NCHMAX/-1,5000/
      DATA RXMIN/0.0/,RXMAX/0.0/
C Define LPT channel number.
      LPT = 8
C Define I/O channel numbers.
C TTY channel.
      IIN=5
      IO = 6
c no fancy jump to different parts of the subroutine
c simply read in the resolution file first
	Do 15 i=1,20
		file(i:i)=' '
15	continue
20    WRITE(IO,40) 
40    FORMAT(/1X,'Resolution filename ? ',$)
      READ(IIN,60) FILE
60    FORMAT(A)
	i=0
17	i=i+1
	if (file(i:i) .ne. ' ') goto 17
	lenfil=i-1
      IF (LENFIL.EQ.0) GOTO 1000
80    NRF = 1
      RFILE = FILE
      LENRF=LENFIL
C Read in the number of spectra and runs etc
      OPEN(UNIT=22,FILE=RFILE,STATUS='OLD',err=800)
c no error checking in this version you are on your own

	read(22,109)cdum
109	format(a)
	read(22,*)ncr
	read(22,109)cdum
	read(22,*)totspc
	read(22,109)cdum
	do 91 i=1,ncr
	   read(22,*)rx(i)
91	continue
	read(22,109)cdum
	do 92 i=1,totspc
	   read(22,*)qspec(i)
92	continue
	do 94 i=1,totspc
	   read(22,109)cdum
  	   do 93 j=1,ncr
	      read(22,*)ytmp(j),dytmp(j)
	      rysp(i,j)=ytmp(j)
		if (dytmp(j) .le. 0.0) then
		   ryersp(i,j)=100000.0
		else
	           ryersp(i,j)=dytmp(j)
		end if
93	   continue
94	continue


c open sample file

120   WRITE(IO,140)
140   FORMAT(/1X,'Sample filename ? ',$)
      READ(IIN,60) FILE
	i=0
147	i=i+1
	if (file(i:i) .ne. ' ') goto 147
	lenfil=i-1
      IF (LENFIL.EQ.0) GOTO 1000
      SFILE=FILE
      LENSF=LENFIL
      OPEN(UNIT=23,FILE=SFILE,STATUS='OLD',ERR=840)
c no error checking in this version you are on your own

	read(23,109)cdum
	read(23,*)ncs
	read(23,109)cdum
	read(23,*)totspc
	read(23,109)cdum
	do 191 i=1,ncs
	   read(23,*)xtmp(i)
191	continue
	read(23,109)cdum
	do 192 i=1,totspc
	   read(23,*)qspec(i)
192	continue
	do 194 i=1,totspc
	   read(23,109)cdum
  	   do 193 j=1,ncs
	      read(23,*)ytmp(j),dytmp(j)
	      sxsp(i,j)=xtmp(j)
	      sysp(i,j)=ytmp(j)
		if (dytmp(j) .le. 0.0) then
		   syersp(i,j)=100000.0
		else
	           syersp(i,j)=dytmp(j)
		end if
193	   continue
194	continue

c kwh 1/2010 end of changes

c a few more checks
C Check that the energy arrays of the two spectra match by looking
C at first and last elements
c      IF (RX(1).NE.xtmp(1).OR.RX(NCR).NE.xtmp(NCS)) THEN
c      WRITE (IO,220) RX(1),RX(NCR),xtmp(1),xtmp(NCS)
c220   FORMAT(1X,'The energy arrays do not match',
c     +     /,1X,'Resolution min and max = ',2(1X,F5.2),
c     +     /,1X,'Sample min and max     = ',2(1X,F5.2))
c      ENDIF
c ok now put thobin into twoth
c      twoth(i)=thobin
c calculate the momentum transfer
c now calculate twoth incase used later from the Q
c pretend final wavelength is 4.5 angstroms always
	do 81 i=1,totspc
		thobin = (asin(qspec(i)/(12.566*4.5)))/0.0087266
		twoth(i)=thobin
		write(6,*)qspec(i),twoth(i)
81	continue
c close the files
      close(22)
      close(23)
C Search for elastic peak in each spectra
      do 241 isp=1,totspc
      rxminsp(isp) = -1.
      rxmaxsp(isp) = 1.
      rymax=0.0
      DO 240 I=1,NCR
      IF (RYsp(isp,I).GT.RYMAX) THEN
      IF (ABS(SXsp(isp,I)).LT.2) THEN
      RYMAX = RYsp(isp,I)
      NEPsp(isp) = I
      ENDIF
      ENDIF
240   CONTINUE
241   continue

C Get limits for resolution function.
      isp=totspc
260   WRITE(IO,280) RXMINsp(isp),RXMAXsp(isp)
280   FORMAT(/,1X,'Current limits for use of resolution function are ',
     +  F7.2,' to ',F7.2,' meV (mmeV)'
     +  /1X,'Type Y to change or <CR> to retain old ',$)
      READ(IIN,300) ANS
300   FORMAT(A1)
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
      WRITE(IO,320)
320   FORMAT(1X,'What are the new values :')
      READ(IIN,*,ERR=260) R1,R2
c modification/correction made by KWH 12/16/2002
	do 327 isp=1,totspc
		RXMINsp(isp)=R1
		RXMAXsp(isp)=R2
327	continue
        isp=totspc
      ENDIF
C     Check if unsuitable - if so repeat question !
340   IF ((RXMAXsp(isp) - RXMINsp(isp)).GT.0.01) GOTO 400
360   WRITE(IO,380)
380   FORMAT(1X,'Values not suitable !')
      GOTO 260
C     Find channel numbers for use with resolution data.
400   CONTINUE

c again must find these as a function of each spectra
c now loop over all spectra and establish the corresponding
c channels for the resolution functions

      do 765 isp=1,totspc
      if (twoth(isp).eq.0) goto 764
      DO 420 I=1,NCR
      NXMIN = I
C     When first limit found jump out of loop
      IF (RXMINsp(isp).LE.SXsp(isp,I)) GOTO 440
420   CONTINUE
440   CONTINUE
      DO 460 I=NXMIN,NCR
      NXMAX = I
C     When second limit found jump out of loop
      IF (RXMAXsp(isp).LT.SXsp(isp,I)) GOTO 480
460   CONTINUE
C Check enough points present for a resolution function - say 10
480   NRPTOT(isp) = NXMAX - NXMIN
      NXMAX = NXMAX - 1
      IF (NRPTOT(isp).GT.9) GOTO 520
      WRITE(IO,500)
500   FORMAT(1X,'Insufficient data to define resolution function')
      GOTO 260
C Check elastic peak is requested !
520   IF (NXMIN.GT.(NEPsp(isp)-4).OR.NXMAX.LT.(NEPsp(isp)+4)) GOTO 360
      LXMINsp(isp) = NEPsp(isp) - NXMIN
      LXMAXsp(isp) = NXMAX - NEPsp(isp)
C Take away the background defined as the lowest value at nxmin or nxmax
      IF(RYsp(isp,NXMAX).LT.RYsp(isp,NXMIN)) GOTO 540
      RYBKG=RYsp(isp,NXMIN)
      GOTO560
540   RYBKG=RYsp(isp,NXMAX)
560   CONTINUE
      IF (RYBKG . LT. 0.) RYBKG=0.0E0
      DO 580 ICHN=NXMIN,NXMAX
      RYsp(isp,ICHN)=RYsp(isp,ICHN)-RYBKG
580   CONTINUE

C Set constant channel width
      DELT = (SXsp(isp,NCR) - SXsp(isp,1)) / FLOAT(NCR - 1)
      DO 640 I=1,NCR
640   WXsp(isp,I) = DELT

660   FACTOR = 0.0
      DO 680 I=NXMIN,NXMAX
      RIsp(isp,I)= WXsp(isp,I) * RYsp(isp,I)
680   FACTOR = FACTOR + WXsp(isp,I) * RYsp(isp,I)
C     Normalize resolution function now.
      DO 700 I=NXMIN,NXMAX
      RIsp(isp,I) = RIsp(isp,I) / FACTOR
      RYsp(isp,I) = RYsp(isp,I) / FACTOR
C and change errors !
700   RYERsp(isp,I) = RYERsp(isp,I) / FACTOR
C Set remainder of resolution function to zero.
      DO 720 I=1,(NXMIN-1)
      risp(isp,i)=0.0
720   RYsp(isp,I) = 0.0
      DO 740 I=(NXMAX+1),NCR
      risp(isp,i)=0.0
740   RYsp(isp,I) = 0.0
C Calculate total error in resolution function
      ERRSUM = 0.0
      DO 760 I=NXMIN,NXMAX
760   ERRSUM = ERRSUM + RYER(I)
      esumsp(isp)=errsum
764   continue
765   continue

c this ends the final loop over all spectra
c

C Write out title and number of points read in for sample.
      WRITE(IO,780) NCS,RXMINsp(totspc),RXMAXsp(totspc),
     +              NrpTOT(totspc)
780   FORMAT(/,1X,I4,' points read in',
     +  /,1X,'Resolution data from ',F6.2,' to ',F5.2,' meV (mmeV)',
     +    1X,' contains ',I5,' points'/)

C     Set number of points to number of channels.
      NRPNT = NCR
      NSPNT = NCS
C     Return to FITFUN.
      RETURN

C************************************************
C      Error handling messages
C************************************************

C If error in opening resolution file...
800   WRITE(IO,820) FILE(1:LENFIL)
820   FORMAT(1X,'Error opening file ',A)
      GOTO 20

C If error in opening resolution file...
840   WRITE(IO,860) FILE(1:LENFIL)
860   FORMAT(1X,'Error opening file ',A)
      GOTO 120

C Error during initial read of resolution file
880   WRITE (IO,900)
900   FORMAT(1X,'Error during initial read of resolution file')
      GOTO 120

C Error during initial read of sample file
920   WRITE (IO,940)
940   FORMAT(1X,'Error during initial read of sample file')
      GOTO 120

C Error while opening print file
c960   WRITE (IO,980)
c980   FORMAT(1X,'Error while opening output file (enough space ?)')
c      PAUSE

C     Return here to FITFUN - escape useful when debugging !
1000  NSPNT = 0
      RETURN
C     *******
      END



      SUBROUTINE FUNCT3(NX1,NX2)
      implicit double precision(a-h,o-z)
C *******************************************************************
C Subroutine which calculates Y values and weights for array elements
C     NX1 to NX2 where 1 <= NX1 < NX2 <= NSPNT.
C********************************************************************
      COMMON/TITLES/PNAMES,TX(4),TY(4)
      COMMON/PARTIT/NPARAS
      COMMON/SPLDAT/NSPNT,SX(5000),SY(5000),SYER(5000)
      COMMON/TEXT1/TEXT,STEXT
      COMMON/CALDAT/CALPAR(100),YCALC(5000),WT(5000)
      COMMON/PARCTL/NumL,IBG,IWT,INST
      COMMON/RESDAT/RY(5000),RYER(5000),WX(5000),RI(5000)
      COMMON/RESLIM/RTOF,Epos,LXMIN,LXMAX,ERRSUM
      common/tofunc/MomTra,MolRad
c parameter list
C Variables for the fit are an overall scale factor, and the
C translational and rotational diffusion constants Diffusion constant
c is in 10^-5 cm^2/sec if energy is meV
C An extra Delta function intensity is needed for background
C contributions etc
c      PNAMES(1)  = 'Scale factor   '
c      PNAMES(2)  = '<u^2>^0.5 sph  '
c      PNAMES(3)  = 'elastic pos    '
c      PNAMES(4)  = 'Extra delta Int'
c      PNAMES(5)  = '<u^2>^0.5 elast'
c      PNAMES(6)  = 'Dsubt (10^-5)  '
c      PNAMES(7)  = 'Sphere Radius  '
c      pnames(8)  = 'Background Cons'
c      PNAMES(9)  = 'Background Line'

      CHARACTER*5 TX,TY,TEXT*40,STEXT*20
      CHARACTER*15 PNAMES(40)
      real*8 TraLor(5000), SI(5000), Ta(5000)
      real*8 MomTra, MolRad
	real*8 eisf,floatl,xnl2,anl
	real*8 D,width
	real*8 xnl(100)
	real*4 a,qa,bj(0:25)
	integer nind(100),lind(100),l
      integer Epos
	common/xcoef/ xnl,nind,lind


      PI = 3.141592654
      a=CALPAR(7)
      d = calpar(6)
      qa=MomTra*a
      call sbess(bj,qa)
      elpos=CalPar(3)
      do 138 i=1,nspnt
      si(i) = 0.0
      ta(i) = 0.0
138   continue

      DO 20 I=NX1,NX2
      YCALC(I) = 0.0
20    continue

      NSFT = int(10000 + CALPAR(3)/WX(1) - 0.5) - 10000
      DELCH = CALPAR(3)/WX(1) - NSFT - 0.5

	eisf=(3.*bj(1)/qa)**2
	const=eisf
c add in delta function
        Si(Epos+NSFT)=
     +       Si(Epos+NSFT)+const*(1-Delch)/wx(1)
        Si(Epos+NSFT+1)=
     +       Si(Epos+NSFT+1)+const*Delch/wx(1)

c now loop through the remaining 98 terms

	do 130 j=1,98
	   l=lind(j)
	   floatl=float(l)
	   xnl2=xnl(j)**2
c probably a problem if xnl(j) gets too close to qa as well
	   if (xnl(j) .eq. qa) then
		anl=1.5*(bj(l)**2)*(xnl2-floatl*(floatl+1.))/xnl2
	   else
		anl=6.*xnl2/(xnl2-floatl*(floatl+1.))
		anl=anl*((qa*bj(l+1)-floatl*bj(l))/(qa**2-xnl2))**2
	   endif
	if (anl .gt. 1000) then
	write(6,*) '***',anl
	do 1333 iijjll=0,25
	write(6,*)iijjll,bj(iijjll)
1333    continue
	end if
	   width=xnl2*D*0.065821/(a*a)
c now have anl let's sum up in the Lorentzian
c if widths are too small then must include as a 
c delta function contribution
	const=anl*(2.*floatl+1.)
      if (width .le. wx(1)/2.) then
        Si(Epos+NSFT)=
     +       Si(Epos+NSFT)+const*(1-Delch)/wx(1)
        Si(Epos+NSFT+1)=
     +       Si(Epos+NSFT+1)+const*Delch/wx(1)
      else
c otherwise add it in as a Lorentzian
        const=const/pi
        DO 127 I=1,NSPNT
        Pos = Sx(I) - ElPos
127     Si(i) = Si(i) +  Const * ALRNTZ(Pos, Const, Width)
      endif
130	continue

c KWH change, multiply by scale factor here, so that when
c add on additional Lorentzian and delta-functions, their 
c intensities can be compared with the other fitting
c programs
      STMP=Calpar(1)*
     +    exp(-1.*MomTra*MomTra*calpar(2)*calpar(2)/3.0)

      Do 114 i=nx1,nx2
        ta(i) = si(i) * STMP
114   Continue

C Add on the extra delta function.
c same problem as above, in order for the fitted
c parameter to be an intensity, you should divide by
c the channel width
      xdum=Calpar(4)*
     +     exp(-1.*MomTra*MomTra*calpar(5)*calpar(5)/3.0)

      Ta(Epos + NSFT) = Ta(Epos + NSFT) + xdum*(1 - Delch)/WX(1)
      Ta(Epos + NSFT+1) = Ta(EPos + NSFT+1) + xdum * Delch/WX(1)

C Add in the extra Lorentzian
c KEEP THIS CODE if NEEDED FOR FUTURE MODIFICATIONS
c      ElPos = CALPAR(5)
c      Width = CalPar(11)
c same problem as above the Lorentzian should not
c be multiplied by the width
c      Const = Wx(1) * CalPar(6) / Pi 
c      Const = CalPar(10) 
c      if (width .le. wx(1)/2.) then
c        Ta(Epos+NSFT)=
c     +       Ta(Epos+NSFT)+const*(1-Delch)/wx(1)
c        Ta(Epos+NSFT+1)=
c     +       Ta(Epos+NSFT+1)+const*Delch/wx(1)
c      else
c otherwise add it in as a Lorentzian
c      const=const/pi
C Get Lorentzian function.
c      DO 142 I=nx1,nx2
c      Pos = Sx(I) - ElPos
c142   Ta(i) = Ta(i) +  Const * ALRNTZ(Pos, Const, Width)
c      endif

c Convolute in the resolution function
      DO 120 I=NX1,NX2
      CONV = 0.0
      DO 100 J=(I-LXMAX),(I+LXMIN)
      IF (J.LT.1.OR.J.GT.NSPNT) GOTO 100
      CONV = CONV + RI(I-J+Epos) * Ta(J)
100   CONTINUE
      YCALC(I) = YCALC(I) + CONV 
120   CONTINUE

      

C     Now add background function.
160   GOTO (200) IBG + 1
      DO 180 I=NX1,NX2
      YCALC(I) = YCALC(I) + CALPAR(8)
c	write(6,*)i,ycalc(i),calpar(8)
      GOTO (180) IBG
      YCALC(I) = YCALC(I) + CALPAR(9) * SX(I)

180   CONTINUE

C     Add error in observed data points to those in the calculated
C     function (which are due to the use of an experimental resolution
C     function).
200   GOTO (240) IWT + 1
C     Use only sample weights for now.
      DO 220 I=NX1,NX2
      WT(I) = SYER(I)
C     THEN CHECK THAT NONE OF THE ERRORS ARE ZERO..
      IF (WT(I).EQ.0.0) WT(I) = 1.0E6
220   CONTINUE
240   RETURN
      END

      double precision FUNCTION ALRNTZ(X,H,W)
      implicit double precision(a-h,o-z)
      
C Function to calculate Lorentzian curve.
      alrntz = w / (w*w + x*x)
      END

C===============================================================================
C
      SUBROUTINE READQL(FN,J,XD,SUMDD,SUMVD,NCHMX,TEMP,THOBIN,EFOBIN,
     *                  NCH,IERR)
C
C     ON ENTRY:  IERR = -1  => PRINT ERROR MESSAGES
C                     =  0  => DO NOT PRINT ERROR MESSAGES
C
C     ON EXIT:   IERR =  0  => NO ERRORS DETECTED
C                     =  1  => FILE NOT FOUND
C                     =  2  => SPECTRUM NUMBER J NOT IN FILE
C                     =  3  => ARRAY IS TOO SMALL TO HOLD ALL OF SPECTRUM
C
C===============================================================================

      implicit double precision(a-h,o-z)
      REAL*8          XD(NCHMX), SUMDD(NCHMX), SUMVD(NCHMX)
      REAL            X(5000), SUMD(5000), SUMV(5000)
      real v1,v2,ttemp
      CHARACTER*(*)   FN

      CHARACTER*1     BLOT
c this subroutine no longer used
	return
c      IERR2 = IERR
c     IERR = 0

c      OPEN(UNIT=22,FILE=FN,FORM='UNFORMATTED',STATUS='OLD',
c    +        recordtype='segmented',convert='vaxd',ERR=80)
c      READ (22,ERR=60,END=60) NBINS, tTEMP, NBN
c      temp=ttemp
c      IF (J.LT.1 .OR. J.GT.NBINS) THEN
c        IERR = 2
c20      FORMAT(' SPECTRUM NUMBER ', I2, ' NOT HELD IN FILE')
c        IF (IERR2.EQ.-1) WRITE (6,20) J
c      ELSE
c        DO I = 1, 4 * (J - 1)
c          READ(22,ERR=60,END=60) 
c        ENDDO
c        READ(22,ERR=60,END=60) NCHO, v1, v2
c        EFOBIN = v1
c        THOBIN = v2
c        NCH = MIN(NCHO, NCHMX)
c        READ(22,ERR=60,END=60) (X(I), I = 1, NCH)
c        READ(22,ERR=60,END=60) (SUMD(I), I = 1, NCH)
c        READ(22,ERR=60,END=60) (SUMV(I), I = 1, NCH)c
c
C Convert to double precision
c        DO 25 I=1,NCH
c        XD(I)=X(I)
c        SUMDD(I)=SUMD(I)
c        SUMVD(I)=SUMV(I)
c25      CONTINUE
c
c        IF (NCH.LT.NCHO) THEN
c          IERR = 3
c40        FORMAT(' ONLY ', I4, ' OF THE ', I4, ' STORED CHANNELS COULD
c     * BE READ FOR BIN ', I2, '.')
c          IF (IERR2.EQ.-1) WRITE(6,40) NCH, NCHX, J
c        ENDIF
c60      CLOSE(UNIT=22)
c     ENDIF
c      RETURN
c
c80    CONTINUE
c      IERR = 1
c100   FORMAT(' FILE NOT FOUND: ', A64)
c      IF (IERR2.EQ.-1) WRITE (6,100) FN
c      RETURN
      END

      subroutine sbess(bj,qa)
      implicit real*4(a-h,o-z)
      real*4 bj(0:25)
      real*4 qa
      INTEGER n
      REAL sj,sjp,sy,syp,x

	x=qa
	do 200 ij=1,26
	   n=ij-1
	   call sphbes(n,x,sj,sy,sjp,syp)
	   bj(n)=sj
200	continue
	return     
	end

C Numerical Recipe Subroutines and functions (F77 versions)

      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      INTEGER n
      REAL sj,sjp,sy,syp,x
CU    USES bessjy
      REAL factor,order,rj,rjp,ry,ryp,RTPIO2
      PARAMETER (RTPIO2=1.2533141)
      if(n.lt.0.or.x.le.0.)pause 'bad arguments in sphbes'
      order=n+0.5
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.*x)
      syp=factor*ryp-sy/(2.*x)
      return
      END

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      INTEGER MAXIT
      REAL rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-10,FPMIN=1.e-30,MAXIT=10000,XMIN=2.,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     *f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,
     *r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
     *temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessjy; try asymptotic expansion'
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        pause 'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=5,NUSE2=5)
CU    USES chebev
      REAL xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371168d0,6.5165112670737d-3,3.087090173086d-4,
     *-3.4706269649d-6,6.9437664d-9,3.67795d-11,-1.356d-13/
      DATA c2/1.843740587300905d0,-7.68528408447867d-2,
     *1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8,2.423096d-10,
     *-1.702d-13,-1.49d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.,1.,c1,NUSE1,xx)
      gam2=chebev(-1.,1.,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END

      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL chebev,a,b,x,c(m)
      INTEGER j
      REAL d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END
