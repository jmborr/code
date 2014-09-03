c 9/21/99 KWH modified for larger arrays to be used with the 
c new QENS upgrade
c this is the version with substantial modifications by K. W. Herwig
c the goal here is to provide a routine which can fit multiple spectra
c (different Q's) simultaneously. In this fashion some parameters of the
c fit may be tied or constrained across the various spectra in order to 
c more readily fit the data. The constraints are handled in the same 
c fashion as individual spectra fits. Namely, each spectra's list of fitting
c parameters is identical to all others, and parameters may be tied across the
c various spectra. Recommendations are that you tie the various parameters to
c those in the earliest spectra fit, that is the lowest number in the overall
c parameter list at which the parameter first appears. You do, however, have 
c the ability to tie parameters between subsets of the total number of 
c fitted parameters. You simply have to pay attention to where in the
c overall parameter list you are.  Note the order of parameters can
c change if you change the number of fitted spectra, BEWARE!!
c 
c new modification 4/21/97
c there is a new command called change class 'd' which changes all of the
c parameters in a given class of parameter throughout the entire list
c of spectra
c new modification 5/1/97
c there is a new command which acts like the 'g' command except
c it affects an entire class of parameters throughout the list
c of spectra
c.............................................................................      
c modifications begun 11/14/95
      SUBROUTINE FITFUN(READIN,CALSUB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL READIN,CALSUB
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
      COMMON/TITLES/PNAMES(40),TX(4),TY(4)
      COMMON/PARTIT/NPARAS
      COMMON/CONFIT/DMAX,ACC,H,MAXNO
      COMMON/TIN/IO,VARIN(10)
      COMMON/FITLIM/XLIM(6),NSET,NLIM(6),NFIT,NPLT
      COMMON/FITPAR/PAR(500),STEP(500),OLD(500),INDEX(500)
      COMMON/SPLDAT/NSPNT,SX(5000),SY(5000),SYER(5000)
      COMMON/RESDAT/RY(5000),RYER(5000),WX(5000),RI(5000)
      COMMON/RESLIM/RTOF,NEP,LXMIN,LXMAX,ERRSUM
      COMMON/TEXT1/TEXT,STEXT
      COMMON/CALDAT/CALPAR(100),YCALC(5000),WT(5000)
      COMMON/PARFIL/PFILE,LENPF
      COMMON/PARFL1/IN
      COMMON/RESDL/REQUED
      COMMON/TOSPLT/FIWT,R1,R2,CURSOR
      common/tofunc/MomTra, MolRad
      COMMON/FILES/RFILE,LENRF,SFILE,LENSF,SPNUM


      LOGICAL REQUED,CURSOR,SAVLOG,IFU
      CHARACTER ANS
      CHARACTER*5 TX,TY,TEXT*40,STEXT*20
      CHARACTER PFILE*20,Filnam*20
      CHARACTER*20 RFILE,SFILE
      INTEGER SPNUM
      character*15 pnames

      real plims(4),pxdata(5000),pydata(5000),presid(5000),pycalc(5000)
      real pyerr(5000)

      REAL*8 FIWT(5000),fi(5000),TEMSAV(500),ESTERR(500),DELCHI(6)
      REAL*8 YOLD(5000)
      real*8 MomTra,MolRad

      INTEGER R1,R2

c      DATA VERSHN/'V2.0'/
      DATA DELCHI/1.00,2.30,3.53,4.72,5.89,7.04/

	write(io,*) "DEBUG Entering FITFUN"
	
      REQUED=.TRUE.
      IFU=.FALSE.

      IN = 10
      IIN=5
      IO = 6
      R1=1
      R2=R1
      WRITE(IO,20)
20    FORMAT(' Intermediate residuals required (<CR>=Yes)?',/' -> ',$)
      READ(IIN,'(A1)') ANS
      IF(ANS.EQ.'N') REQUED=.FALSE.
C     Set hardcopy flag to zero.
      IHC = 0
C     Set "Q" command flag to zero.
      IFQ = 0
C     Set number of points to zero.
      NSPNT = 0
      NFIT = 0
C     Initialize parameters before reading old ones in.

c now must loop over all spectra
      do 65 isp=1,totspc
      NSETsp(isp) = 1
      DO 40 I=2,6,2
      XLIMsp(isp,I-1) = 0.0
40    XLIMsp(isp,I) = 0.0
      DO 60 I=1,2
      XSCALE(isp,I) = 0.0
60    YSCALE(isp,I) = 0.0
65    continue

      DO 80 I=1,NPARAS
      PAR(I) = 0.0
80    STEP(I) = 1.0
C     Set constants required by VAO5A.
      DMAX = 100.0
      ACC = 0.01
      H = 0.1
      MAXNO = 100
C     Get filename for use in dumping parameters.
      WRITE(IO,100)
100   FORMAT(' Name of parameter file (<CR>=CON) ?',/' -> ',$)
c need parsing routine for pfile
	do 101 iijj=1,20
		pfile(iijj:iijj)=' '
101	continue

      READ(IIN,120) PFILE
120   FORMAT(A)
	iijj=20
111	iijj=iijj-1
	if (pfile(iijj:iijj) .eq. ' ') goto 111
	lenpf=iijj
	

      IF (LENPF.EQ.0) THEN
      LENPF=7
      PFILE(1:LENPF)='CON.PAR'
      ELSE
      PFILE(LENPF+1:LENPF+4)='.PAR'
      LENPF=LENPF+4
      ENDIF
C     Open old parameter file and use its parameters if it exists.
      OPEN(UNIT=IN,ACCESS='SEQUENTIAL',FILE=PFILE(1:LENPF),ERR=200,
     &	 STATUS='UNKNOWN')
C     Read in old "Only" limits.
      read(in,151,err=200,end=200)totspc,newspc,npspec
151   format(1x,3i4)
      isptot=0
      iftpar=0


c now must loop over all spectra
c store parameters etc into a master parameter list
c when get to fitting routine is where you change to list
c required by fitting algorithms
c
      ipar=0
      do 185 isp=1,totspc
c check integer flag if 1 then this spectra is in the fitting list
      read(in,*)iflag
	
	write(io,*) "DEBUG iflag=",iflag

      if (iflag .eq. 1) then
         isptot=isptot+1
         ispec(isptot)=isp
      endif
      read(in,*)nsetsp(isp)
      DO 140 I=2,6,2
140   READ(IN,161,ERR=200,END=200) XLIMsp(isp,I-1),XLIMsp(isp,I)
      nfitsp(isp)=0
      do 141 i=1,nsetsp(isp)
141   NFITsp(isp) = NFITsp(isp) + NLIMsp(isp,i+1) - NLIMsp(isp,i)
160   FORMAT(1X,4G12.4)
161   FORMAT(1X,2G12.4)
C     Read in old plotting limits.
      READ(IN,161,END=200,ERR=200) XSCALE(isp,1),XSCALE(isp,2)
      READ(IN,161,END=200,ERR=200) YSCALE(isp,1),YSCALE(isp,2)
C     Read in old parameter values - all of them !
      DO 180 I=1,npspec
      ipar=ipar+1
c read in values into a global par and step array
      READ(IN,161,ERR=200,END=200) PARsp(ipar),STEPsp(ipar)
c also set up par and step arrays
      if (iflag .eq. 1) then
      iftpar=iftpar+1
      par(iftpar)=parsp(ipar)
      step(iftpar)=stepsp(ipar)
      endif
180   continue
185   continue
c reset fitspc and nparas
      fitspc=isptot

	write(io,*) "DEBUG fitspc=",fitspc

      nparas=fitspc*npspec

c check to make sure read in everything correctly
      ntemp=npspec*totspc
      if (ntemp .ne. ipar) then
        write(6,191)
191    format(1x,'number of total input parameters in PAR file does')
       write(6,192)
192    format(1x,'not agree with current number of parameters, 
     + proceed with caution')
       write(6,187)ipar,ntemp
187    format(1x,' counted number:',i5,' calc. number:',i5)
      endif
      if (nparas .ne. iftpar) then
       write(6,197)
197    format(1x,'number of total fitted parameters in PAR file does')
       write(6,192)
      endif
c check to make sure current fitspc is the same as the one selected
      if (newspc .ne. fitspc) then
       write(6,193)
193    format(1x,'number of fitted spectra in PAR file disagrees 
     +with')
       write(6,194)
194    format(1x,'desired number of fitted spectra, need to use 
     + C command')
      endif

C     Close old parameter file.
200   CLOSE(UNIT=IN)
C     Set old parameter values to zero.
      DO 220 I=1,NPARAS
220   OLD(I) = 0.0

C Clear screen and leave in graphics mode.
      CALL FTERAS

C Give command list first time round.
C     CALL FTHELP

C Get command.
240   CALL FTCMD(IVAL)
C Commands are H,A,X,Y,O,V,R,P,F,B,Q,L,S,C,D,E,G,W
260   IF (IVAL.LE.0) GOTO 1080
      GOTO (280,300,320,340,360,400,420,440,540,860,900,980,1000,1020,
     &	 1040,1080,1060,1100),IVAL
      GOTO 240

C If "H" command.
280   CALL FTERAS
      CALL FTHELP
C Get next command.
      GOTO 240

C If "A" command. 
300   CALL FTERAS
c FTVALS rewritten 2/13/96
      CALL FTVALS
      GOTO 240

C If "X" command.
c FTXSET rewritten 2/13/96
320   CALL FTXSET
C Reset "Q" command flag.
      IFQ = 0
      GOTO 240

C If "Y" command.
c modified FTYSET 2/13/96
340   CALL FTYSET
C Reset "Q" command flag.
      IFQ = 0
      GOTO 240

C If "O" command.
360   continue
      isp=Int(varin(1))
      DO 380 I=1,6
380   XLIMsp(isp,I) = VARIN(I+1)
C Determine ranges and number of points.
c FTLIMT modified 2/14/96
      CALL FTLIMT(isp)
      GOTO 240

C If "V" command.
c FTNEWP modified 2/14/96
400   CALL FTNEWP(ESTERR)
      GOTO 240

C If "R" command.
420   CALL FTERAS
      CALL READIN
      do 425 isp=1,totspc
      CALL FTDATN(isp)
c ftlimt modified 2/14/96
      CALL FTLIMT(isp)
      varin(1)=isp
      VARIN(2) = XSCALE(isp,1)
      VARIN(3) = XSCALE(isp,2)
      CALL FTXSET
425   continue
      GOTO 240

C If "P" command.
C If no points - no plot !
c first of all ask which of the fitted spectra you wish to plot
440   write(6,441)
441   format(1x,'enter spectra to plot:',$)
      read(5,*)isp
c note that isp may not be in the fitted list and should
c be treated differently in the two cases
c search fitted list to see
      ifit=0
      do 442 i=1,fitspc
      if (isp .eq. ispec(i)) ifit = 1
442   continue
c440   GOTO (460),NSPNT + 1
c      GOTO (500) NPLT + 1
      CALL FTERAS
C-----------------------------------------------------------------------
C Put in a special call to splt here to take care of C language version
C-----------------------------------------------------------------------
c need to determine r2 and r1 here in this routine
      R1=NLIMsp(isp,1)
      R2=NLIMsp(isp,NSETsp(isp)*2)
c also create the original YCALC and FIWT arrays
      do 118 i=1,nspnt
        ycalc(i)=ycalcsp(isp,i)
        fiwt(i)=fiwtsp(isp,i)
c select out the current data set
        sx(i)=sxsp(isp,i)
        sy(i)=sysp(isp,i)
c set up WT array only used in plotting program and
c in Weighting the residuals in the fitting part of this routine
        wt(i)=syersp(isp,i)
118   continue
c
c	write(6,*)' You got this far 1'
      numplt=r2-r1+1
      j=0
      do 123 i=r1,r2
      j=j+1
      pxdata(j)=sxsp(isp,i)
      pydata(j)=sysp(isp,i)
      pyerr(j)=syersp(isp,i)
c these are different if no fit
      if (ifit .eq. 1) then
        presid(j)=fiwtsp(isp,i)
        pycalc(j)=ycalcsp(isp,i)
      endif
123   continue
      plims(1)=xscale(isp,1)
      plims(2)=xscale(isp,2)
      plims(3)=yscale(isp,1)
      plims(4)=yscale(isp,2)
C     call splt(plims,pxdata,pydata,pyerr,presid,pycalc,numplt,1)
c provide the ID number to the plotting routine
      spnum=isp

      CALL SPLT(XSCALE(isp,1),XSCALE(isp,2),YSCALE(isp,1),
     +          YSCALE(isp,2),SX,SY,NSPNT)
C      CALL ENDGR(0)

      CALL FTCMD(IVAL)
      CALL FTERAS
      GOTO 260
C Give error message when no points.
460   WRITE(IO,480)
480   FORMAT(1X,'No data !')
      GOTO 240

C Give error message if no points in current X range.
c500   WRITE(IO,520)
c520   FORMAT(1X,'No data for plotting !')
c      GOTO 240

C If "F" command.
C If no points - no fit !
540   GOTO (460) NSPNT + 1
C Check number of parameters !
      IPAR = 0
c
c initialize nfit
      nfit=0
	write(IO,*) "debug: fitspc",fitspc
      do 541 i=1,fitspc
      isp=ispec(i)
	write(IO,*) "debug: isp",isp," nfitsp(isp)",nfitsp(isp)
	nfit=nfit+nfitsp(isp)
541   continue
c
      DO 560 I=1,NPARAS
560   IF (STEP(I).GT.0.0) IPAR = IPAR + 1
      IF (NFIT.GT.IPAR) GOTO 600

      WRITE(IO,580) NFIT,IPAR
580   FORMAT(1X,'Fitting ',I4,' points with ',I3,' parameters',
     +  1X,'is strictly forbidden')
      GOTO 240

C Reset the flag for the error estimation
600   CONTINUE
      DO 610 I=1,NPARAS
      ESTERR(I)=0.0
610   CONTINUE
      MAXFUN = INT(VARIN(1))
      LOOP=INT(VARIN(2))
C If no iterations required simply calculate residuals and YCALC
C values.
C Otherwise make current parameter values into old parameter values.
      IF (MAXFUN.LE.0) GOTO 680
      IF (LOOP.LE.0) LOOP=1
      DO 660 ILOOP=1,LOOP
      DO 620 I=1,NPARAS
620   OLD(I)=PAR(I)
      CALL FTERAS
      WRITE(IO,640) ILOOP,LOOP
640   FORMAT(' FIT NO. ',I3,' STARTED OUT OF ',I3)
      CALL FTFFFF(CALSUB,MAXFUN,ESTERR)
660   CONTINUE
c no plot anyway only allowing for plots now within the
c "P" command as inquiry must be made about which spectra
c to plot, note must calculate the residuals etc. for all
c spectra so that they are available to the "Plot" command
C If no points - no plot !
c680   GOTO (500) NPLT + 1
680   continue

      RESID = 0.0
      ACF = 0.0
      NACF=0

c loop over all fitted spectra
      do 9200 iij=1,fitspc
c select the fitted spectra
      isp=ispec(iij)
      DO 700 I=1,NPspec
c after fit transfer a copy of current fit values into global arrays
c so that if number of spectra is changed after this
c results are not lost
      parsp((isp-1)*npspec+i)=PAR(I+(iij-1)*npspec)
      stepsp((isp-1)*npspec+i)=step(I+(iij-1)*npspec)
c initialize paramters for CALSUB for this spectra
700   CALPAR(I) = PAR(I+(iij-1)*npspec)

c ok now we need to create all the arrays and elements 
c necessary for calculating the fit etc.

      lxmin=lxminsp(isp)
      lxmax=lxmaxsp(isp)
      nep=nepsp(isp)
      errsum=esumsp(isp)
c also need MOMTRA
      Momtra=qspec(isp)

c create arrays needed 

      do 710 i=1,nspnt
        sx(i)=sxsp(isp,i)
        sy(i)=sysp(isp,i)
        syer(i)=syersp(isp,i)
        ry(i)=rysp(isp,i)
        ryer(i)=ryersp(isp,i)
        wx(i)=wxsp(isp,i)
        ri(i)=risp(isp,i)
710   continue

c ok now call CALSUB just like normal

      CALL CALSUB(KXMIN(isp),KXMAX(isp))

c store ycalc in ycalcsp so it is available for plotting 
      do 715 i=kxmin(isp),kxmax(isp)
      ycalcsp(isp,i)=ycalc(i)
715   continue

c calculate the incremental to the residuals
c as well as the array used in SPLT

      DO 720 I=1,5000
      FIWTsp(isp,I)=0
720   CONTINUE
      R1=NLIMsp(isp,1)
      R2=NLIMsp(isp,NSETsp(isp)*2)
      DO 760 J=1,NSETsp(isp)
      J2 = J * 2 - 1
      DO 740 I=NLIMsp(isp,J2),NLIMsp(isp,J2+1)
      FI(I) = (SY(I) - YCALC(I))
      FIWTsp(isp,I) = FI(I) / syersp(isp,i)
740   RESID = RESID + FIWTsp(isp,I) * FIWTsp(isp,I)
760   CONTINUE

C     Calculate auto-correlation function.

      NACF = nacf+NFITsp(isp)
      DO 800 J=1,NSETsp(isp)
      J2 = J * 2 - 1
C If only one point in range - jump.
      GOTO (800) NLIMsp(isp,J2) - NLIMsp(isp,J2+1)
      NACF = NACF - 1
      DO 780 I=NLIMsp(isp,J2),(NLIMsp(isp,J2+1)-1)
      TEMP=ABS(SYER(I)-SYER(I+1))
      IF (TEMP.GT.0) THEN
      ACF = ACF + ABS(FI(I) - FI(I+1)) / TEMP
      ENDIF
780   CONTINUE
800   CONTINUE

c end loop over all spectra
9200  CONTINUE

C     Divide residuals by N - P + C
      RESID = RESID / (NFIT - IPAR)

      ACF = ACF / NACF

c no plotting here, must call plot subroutine to plot
      goto 820
C Erase screen and plot latest !
C      CALL FTERAS
C-----------------------------------------------------------------------
C Put in a special call to splt here to take care of C language version
C-----------------------------------------------------------------------
c      write (6,*) r1,r2
c      numplt=r2-r1+1
c      j=0
c      do 124 i=r1,r2
c      j=j+1
c      pxdata(j)=sx(i)
c      pydata(j)=sy(i)
c      pyerr(j)=syer(i)
c      presid(j)=fiwt(i)
c      pycalc(j)=ycalc(i)
c124   continue
c      plims(1)=xscale(1)
c      plims(2)=xscale(2)
c      plims(3)=yscale(1)
c      plims(4)=yscale(2)
C     call splt(plims,pxdata,pydata,pyerr,presid,pycalc,numplt,1)
c
c      CALL SPLT(XSCALE(1),XSCALE(2),YSCALE(1),YSCALE(2),SX,SY,NSPNT)
820   CONTINUE
      WRITE(IO,840) RESID,ACF
840   FORMAT(1X,'Value of the residuals is ',G12.4
     +  /1X,'Goodness of fit is ',G12.4)
C Get next command.
      CALL FTCMD(IVAL)
      CALL FTERAS
      GOTO 260

C If "B" command.
C If no points - no plot !
c "B" is a null command not active in this version
c return to sender
860   Write(6,861)
861   format(1x,'NULL command return to menu')
      goto 240

c860   GOTO (460) NSPNT + 1
c      GOTO (500) NPLT + 1
C Transfer parameters to array used in CALSUB.
c      DO 880 I=1,NPARAS
c880   CALPAR(I) = PAR(I)
c      CALL CALSUB(KXMIN,KXMAX)
C-----------------------------------------------------------------------
C Put in a special call to splt here to take care of C language version
C-----------------------------------------------------------------------
c      numplt=r2-r1+1
c      j=0
c      do 125 i=r1,r2
c      j=j+1
c      pxdata(j)=sx(i)
c      pydata(j)=sy(i)
c      pyerr(j)=syer(i)
c      presid(j)=fiwt(i)
c      pycalc(j)=ycalc(i)
c125   continue
c      plims(1)=xscale(1)
c      plims(2)=xscale(2)
c      plims(3)=yscale(1)
c      plims(4)=yscale(2)
c     call splt(plims,pxdata,pydata,pyerr,presid,pycalc,numplt,2)
c
c      CALL HPLT(XSCALE(1),XSCALE(2),YSCALE(1),YSCALE(2),SX,SY,NSPNT)
C Set hardcopy flag to 1 so that ENDPLT is called.
c      IHC = 1
C Set flag for use with "Q" command.
c      IFQ = 1
C Get next command.
c      GOTO 240

C If "Q" command.
c this is a change step variable for class of parameters command
c this command is now rewritten as a change
c class of variable command
900   continue
      Write(6,1041)
      write(6,1042)
      read(5,*)ip
      if ((ip .le. 0) .or. (ip .gt. npspec)) then
c return
         write(6,1043)
         goto 240
      endif
c otherwise a valid command so change all parameters even in spectra
c not currently being fit, change master list
      write(6,944)
944   format(1x,'enter new value for step:',$)
      write(6,1047)
      read(5,*)istep
      do 948 isp=1,totspc
         ipar=(isp-1)*npspec+ip
         step(ipar)=istep
948   continue
      if (istep .lt. 0) then
      write(6,1049)
      endif
      goto 240

C If "L" command.
c this routine must be modified
980   CONTINUE
      CALL FTLIST(CALSUB,ESTERR)
      GOTO 240

C If "S" command.
c this command must be modified
1000  CALL FTSAVE
      GOTO 240

C If "C" command.
c "C" command is now the command for choosing spectra to fit
c via a new subroutine lspec
1020  call lspec
c must now resequence parameter info from master list 
      do 1021 iij=1,fitspc
c select the fitted spectra
      isp=ispec(iij)
      DO 1019 I=1,NPspec
        PAR(I+(iij-1)*npspec)=parsp((isp-1)*npspec+i)
        step(I+(iij-1)*npspec)=stepsp((isp-1)*npspec+i)
1019  continue
1021  continue
      nparas=npspec*fitspc
c initialize paramters for CALSUB for this spectra
c1020  CURSOR=.TRUE.
      GOTO 240

C If "D" command.
c kkk
c this command is now rewritten as a change
c class of variable command
1040  continue
      Write(6,1041)
1041  format(1x,'this command varies an entire class of variables')
      write(6,1042)
1042  format(1x,'enter paramter number as appears in first spectra:',$)
      read(5,*)ip
      if ((ip .le. 0) .or. (ip .gt. npspec)) then
c return
         write(6,1043)
1043     format(1x,'return to main menu')
         goto 240
      endif
c otherwise a valid command so change all parameters even in spectra
c not currently being fit, change master list
      write(6,1044)
1044  format(1x,'enter new value for the parameters:',$)
      read(5,*)xval
      write(6,1046)
1046  format(1x,'enter value for step:')
      write(6,1047)
1047  format(1x,'0 to fix, 1 to vary, (-1*par#) to tie:',$)
      read(5,*)istep
      do 1048 isp=1,totspc
         ipar=(isp-1)*npspec+ip
         par(ipar)=xval
         step(ipar)=istep
1048  continue
      if (istep .lt. 0) then
      write(6,1049)
1049  format(1x,'check the value of STEP for first parameter')
      endif
      goto 240
c     CURSOR=.FALSE.
c      GOTO 240

C If "G" command
c this command works exactly as it always did
1060  STEP(INT(VARIN(1)))=VARIN(2)
      GOTO240

C If "W" command
c this command slightly modified 2/21/96
c will only write out particular identified spectra
1100  continue
      write (6,1160)
1160  format(1x,'Filename : ',$)
      read (5,1200) Filnam
      write(6,1161)
1161  format(1x,'enter spectra to write out:',$)
      read(5,*)isp
1200  format(a)
      open(11,file=filnam,status='unknown')

c initialize r1,r2
      R1=NLIMsp(isp,1)
      R2=NLIMsp(isp,NSETsp(isp)*2)      
      do 1220 i=r1,r2
      write (11,1240) sxsp(isp,i),sysp(isp,i),syersp(isp,i)
     +                ,ycalcsp(isp,i)
1240  format(4(1x,e10.3e2))
1220  continue
      close(11)
      goto 240

C If "E" command.
c will need to modify FTSAVE
1080  WRITE(IO,1120)
1120  FORMAT(1X,'Do you wish to save the current parameter values ?',
     +  1X,' (<CR> = no)'/1X,'->  ',$)
      READ(IIN,1140) ANS
1140  FORMAT(A1)
      IF (ANS.NE.' ') CALL FTSAVE
C     If hardcopies requested then call ENDPLT.
C      CALL DONEPL
C      IF (IHC.EQ.1) CALL ENDPLT
C     Now return user to main program !
      RETURN
C     *******
      END

c FTNEWP modified for multi-spectra fitting 2/17/96
      SUBROUTINE FTNEWP(ESTERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/TITLES/PNAMES(40),TX(4),TY(4)
      COMMON/PARTIT/NPARAS
      COMMON/CONFIT/DMAX,ACC,H,MAXNO
      COMMON/FITPAR/PAR(500),STEP(500),OLD(500),INDEX(500)
      COMMON/TIN/IO,VARIN(10)
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


      CHARACTER*5 TX,TY
      CHARACTER*15 PNAMES
      REAL*8 ESTERR(500)
C     Get parameter number.
      IP = INT(VARIN(1))
C     Give values if IP negative or zero.
      IF (IP.LE.0) GOTO 140
C     If special parameter then not stored in PAR !
      IF (IP.GT.NPARAS) GOTO 20
C     If constraint present check its possible !
      IF (VARIN(3).LT.0.0.AND.INT(ABS(VARIN(3))).GT.NPARAS) GOTO 100
      PAR(IP) = VARIN(2)
      STEP(IP) = VARIN(3)
C     Parameters set so return.
      RETURN
20    ISP = IP / 100
      GOTO (40,60,80) ISP
C     If not a special parameter number then give message.
      GOTO 220
40    DMAX = VARIN(2)
      RETURN
60    ACC = VARIN(2)
      RETURN
80    H = VARIN(2)
      RETURN
C     Error message if large negative step size.
100   WRITE(IO,120) VARIN(3)
120   FORMAT(1X,G12.4,' is not a valid constraint !')
      RETURN
C     Give values.
140   IF (IP.LT.0) GOTO 180
c this is the print paramter list section
c supposed to list all fitted paramters
c do so in the same manner as in FTVALS
      CALL FTERAS
      ipar=0
      do 845 iij=1,fitspc
      isp=ispec(iij)
      write(io,839) isp
839   format(1x,'parameters for spectra id number: ',i3)
      write(io,840)
840   format(1x,'Number',5x,'Name',6x,'Value',5x,'( Old value  )',
     +   3x,'Step',9x,'Esd')
      do 843 j=1,npspec
         ipar=ipar+1
         write(io,842)ipar,PNAMES(j),par(ipar),old(ipar),step(ipar),
     +                esterr(ipar)
842      format(1x,i4,3x,a15,2x,g12.4,'(',g12.4,')',g12.4,2x,g12.4)
843   continue
845   continue
      RETURN
C     Write out values requested.
180   IP = -IP
      IF (IP.GT.NPARAS) GOTO 220
c only tricky thing here is to get the correct name listed
      iij=mod(ip,npspec)
      WRITE(IO,200) IP,PNAMES(Iij),PAR(IP),OLD(IP),STEP(IP)
200   FORMAT(1X,I4,3X,A15,2X,G12.4,'(',G12.4,')',G12.4/)
      RETURN
C     Give error if IP too large.
220   WRITE(IO,240) IP
240   FORMAT(1X,I5,' is not a valid parameter number !')
      RETURN
C     *******
      END
c
c
c FTSAVE modified by KWH 2/1/96
      SUBROUTINE FTSAVE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/PARFIL/PFILE,LENPF
      COMMON/PARFL1/IN
      COMMON/FITPAR/PAR(500),STEP(500),OLD(500),INDEX(500)
      COMMON/FITLIM/XLIM(6),NSET,NLIM(6),NFIT,NPLT
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

      CHARACTER*10 PFILE*20
	character*10 QFILE*20
      character*15 PNAMES

C     This writes to file PFILEthe current values of the
C     "O", "X", "Y", and "V" commands.
      IW = IN
	IW = 9
	QFILE(1:LENPF)=PFILE(1:LENPF)
	QFILE(LENPF+1:LENPF+4)='.PAR'
C     Open file for output.
      OPEN(UNIT=IW,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &	 FILE=QFILE(1:LENPF+4))

      write(iw,151)totspc,fitspc,npspec
151   format(1x,3i4)
      isptot=0
      iftpar=0

c now must loop over all spectra
c store parameters etc into a master parameter list
c when get to fitting routine is where you change to list
c required by fitting algorithms
c
      do 185 isp=1,totspc
c check if this spectra is in the ispec list
c if so modify the flag in the xxx.PAR file
      iflag=0
      do 35 i=1,fitspc
        if (ispec(i) .eq. isp) iflag=1
35    continue
      write(iw,*)iflag
      write(iw,*)nsetsp(isp)
      DO 140 I=2,6,2
140   write(iw,161) XLIMsp(isp,I-1),XLIMsp(isp,I)
160   FORMAT(1X,4G12.4)
161   FORMAT(1X,2G12.4)
C     write out old plotting limits.
      write(Iw,161) XSCALE(isp,1),XSCALE(isp,2)
      write(Iw,161) YSCALE(isp,1),YSCALE(isp,2)
C     write out old parameter values - all of them !
      DO 180 I=1,npspec
      ipar=(isp-1)*npspec+i
c write in values into a global par and step array
      write(iw,161) PARsp(ipar),STEPsp(ipar)
180   continue
c end loop over all spectra
185   continue
C     Close output file.
      CLOSE(UNIT=IW)
      RETURN
C     *******
      END
c FTDATN modified 2/14/96
      SUBROUTINE FTDATN(isp)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/SPLDAT/NSPNT,SX(5000),SY(5000),SYER(5000)
      COMMON/TEXT1/TEXT,STEXT
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

      CHARACTER*40 TEXT,STEXT*20
C     If no points readin don't change limits.
      IF (NSPNT.LT.2) GOTO 120
      XMIN(isp) = SXsp(isp,1)
      XMAX(isp) = SXsp(isp,NSPNT)
      YMIN(isp) = SYsp(isp,1)
      YMAX(isp) = YMIN(isp)
      DO 20 I=2,NSPNT
C     Check that X values are already sorted from -ve to + ve.
C     If not set error flag to one.
      IF (SXsp(isp,I).LE.SXsp(isp,I-1)) GOTO 80
C     Get limits for Y values.
      TMP = SYsp(isp,I)
      IF (TMP.LT.YMIN(isp)) YMIN(isp)=TMP
20    IF (TMP.GT.YMAX(isp)) YMAX(isp)=TMP
C     If no plot limits set then might as well take defaults from
C     data limits.
      IF (XSCALE(isp,1).NE.XSCALE(isp,2)) GOTO 40
      XSCALE(isp,1) = XMIN(isp)
      XSCALE(isp,2) = XMAX(isp)
40    IF (YSCALE(isp,1).NE.YSCALE(isp,2)) GOTO 60
      YSCALE(isp,1) = YMIN(isp)
      YSCALE(isp,2) = YMAX(isp)
60    RETURN
C     If x values not in ascending order, ignore data (for now).
80    WRITE(6,100)
100   FORMAT(1X,'Data set not suitable - x values not in ascending',
     +  1X,'order')
C     Set number of points back to zero.
120   NSPNT = 0
      RETURN
C     *******
      END
c modified FTLIMT for multi-spectra fitting 2/14/96
      SUBROUTINE FTLIMT(isp)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/SPLDAT/NSPNT,SX(5000),SY(5000),SYER(5000)
      COMMON/TEXT1/TEXT,STEXT
      COMMON/TIN/IO,VARIN(10)
      COMMON/FITLIM/XLIM(6),NSET,NLIM(6),NFIT,NPLT

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


      CHARACTER*40 TEXT,STEXT*20
C     Determine number of ranges set by "O" command.
      NSETsp(isp) = 0
      DO 20 J=2,6,2
C     If no range set or stupid range set take default.
      IF (XLIMsp(isp,J).LT.XLIMsp(isp,J-1)) GOTO 40
      IF (XLIMsp(isp,J).EQ.XLIMsp(isp,J-1).AND.NSETsp(isp).EQ.0) 
     +   GOTO 40
20    IF (XLIMsp(isp,J).GT.XLIMsp(isp,J-1)) 
     +    NSETsp(isp) = NSETsp(isp) + 1
      GOTO 60
C     If no ranges set, take full range.
40    NSETsp(isp) = 1
      XLIMsp(isp,1) = XMIN(isp)
      XLIMsp(isp,2) = XMAX(isp)
      NFITsp(isp) = NSPNT
      NLIMsp(isp,1) = 1
      NLIMsp(isp,2) = NSPNT
      GOTO 220
C     Check that ranges do not overlap.
60    GOTO (100) NSETsp(isp)
      DO 80 J=2,NSETsp(isp)
      J2 = J * 2 - 1
C     If ranges overlap, take full range.
      IF (XLIMsp(isp,J2).LT.XLIMsp(isp,J2-1)) GOTO 40
80    CONTINUE
C     Get total number of points.
100   NFITsp(isp) = 0
C     If no data readin the set return.
      GOTO (220) NSPNT + 1
      IBEG = 1
      DO 200 J=1,NSETsp(isp)
      J2 = J * 2 - 1
      DO 120 I=IBEG,NSPNT
      NLIMsp(isp,J2) = I
C     When lower limit found jump out of loop.
      IF (XLIMsp(isp,J2).LE.SXsp(isp,I)) GOTO 140
120   CONTINUE
140   DO 160 I=NLIMsp(isp,J2),NSPNT
      NLIMsp(isp,J2+1) = I
C     When upper limit found jump out of loop.
      IF (XLIMsp(isp,J2+1).LT.SXsp(isp,I)) GOTO 180
160   CONTINUE
      NLIMsp(isp,J2+1) = NSPNT + 1
180   NFITsp(isp) = NFITsp(isp) + NLIMsp(isp,J2+1) - NLIMsp(isp,J2)
      IBEG = NLIMsp(isp,J2+1)
200   NLIMsp(isp,J2+1) = NLIMsp(isp,J2+1) - 1
C     Write out number of points.
220   WRITE(IO,240) NFITsp(isp),isp
240   FORMAT(1X,'There are ',I4,' points in current fitting range for
     + spectra',i3/)
      RETURN
C     *******
      END
      SUBROUTINE FTVALS
c modified 2/13/96
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/TITLES/PNAMES(40),TX(4),TY(4)
      COMMON/PARTIT/NPARAS
      COMMON/FITPAR/PAR(500),STEP(500),OLD(500),INDEX(500)
      COMMON/CONFIT/DMAX,ACC,H,MAXNO
      COMMON/FITLIM/XLIM(6),NSET,NLIM(6),NFIT,NPLT
      COMMON/TIN/IO,VARIN(10)
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

      CHARACTER*5 TX,TY
      CHARACTER*15 PNAMES
      WRITE(IO,20) TY,TX
20    FORMAT(1X,'Parameters in fitting Y : ',4A5/16X,'Versus X : ',4A5)
      WRITE(IO,60) DMAX,ACC,H
60    FORMAT(1X,' 100',3X,'Maximum step',G12.4
     +  /1X,' 200',3X,'Accuracy    ',G12.4
     +  /1X,' 300',3X,'Step width  ',G12.4)
c loop over all fitted spectra
      ipar=0
      do 45 iij=1,fitspc
      isp=ispec(iij)
      write(io,39) isp
39    format(1x,'parameters for spectra id number: ',i3)
      write(io,40)
40    format(1x,'Number',5x,'Name',6x,'Value',5x,'( Old value  )',
     +   3x,'Step')
      do 43 j=1,npspec
         ipar=ipar+1
         write(io,42)ipar,PNAMES(j),par(ipar),old(ipar),step(ipar)
42       format(1x,i4,3x,a15,2x,g12.4,'(',g12.4,')',g12.4)
43    continue
      WRITE(IO,80) XMIN(isp),XMAX(isp),YMIN(isp),YMAX(isp),
     +   XSCALE(isp,1),xscale(isp,2),YSCALE(isp,1),yscale(isp,2)
80    FORMAT(/1X,'Current data limits : Xmin = ',G12.4,' Xmax = ',
     +  G12.4/23X,'Ymin = ',G12.4,' Ymax = ',G12.4
     +  /1X,'Plot limits',9X,': Xmin = ',G12.4,' Xmax = ',G12.4
     +  /23X,'Ymin = ',G12.4,' Ymax = ',G12.4)
      WRITE(IO,100) NFITsp(isp)
100   FORMAT(/1x,'There are ',I5,' points in the current',
     +  1X,'range set by the "O" command')
      I2 = NSETsp(isp) * 2
      WRITE(IO,120) (XLIMsp(isp,I-1),XLIMsp(isp,I),I=2,I2,2)
120   FORMAT(1X,'Current limits are :',
     +  1X,G12.4,'  to  ',G12.4,2(/22X,G12.4,'  to  ',G12.4)/)
c concludes loop over all fitted spectra
45    continue
      RETURN
C     *******
      END
      SUBROUTINE FTHELP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/TIN/IO,VARIN(10)
      CHARACTER*5 DD(2),TT
C     Get date and time.
c      CALL DATE(DD)
c      CALL TIME(TT)
C     Type out options available.
c      WRITE(IO,20)DD,TT
	write(io,20)
20    FORMAT(10X,8X,/
     +  /1X,'Basic fitting and control commands',
     +  /1X,'H                type this list',
     +  /1X,'A                verify all current values',
     +  /1X,'X min max        change current x-scale to min<x<max',
     +  /1X,'Y min max        change current y-scale to min<y<max',
     +  /1X,'O I a b c d e f  only fit spectrum I between range a<x<b',
     +  /1X,'                 and c<x<d etc. (ignored if both are zero)',
     +  /1X,'V n v s          change parameter number n to v with',
     +  /1X,'                 associated step s (about 1.0) if this is',
     +  /1X,'                 to be fitted.',
     +  /1X,'                 if s is negative parameter n is tied',
     +  /1X,'                 to preceding parameter #-s',
     +  /1X,'R                read in fresh data')
      WRITE(IO,40)
40    FORMAT(1X,'P                plot data',
     +  /1X,'F                plot data and calculated function',
     +  /1X,'F n m            fit for n cycles repeat m times and plot',
     +  /1X,'B                NULL command in this version',
     +  /1X,'                                             ',
     +  /1X,'Q                change class of parameters step value',
     +  /1X,'                 previous hard copy plot (be careful)',
     +  /1X,'L                list summary of current fit to LPT',
     +  /1X,'S                save current parameters',
     +  /1X,'C                choose the subset of fitted spectra',
     +  /1X,'D                change class of parameters',
     +  /1X,'E                EXIT to main program.',
     +  /1X,'G m v            Change parameter m''s step value to v'/)
      RETURN
C     *******
      END

      SUBROUTine FTXSET
c FTXSET rewritten 2/13/96

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/FITLIM/XLIM(6),NSET,NLIM(6),NFIT,NPLT
      COMMON/SPLDAT/NSPNT,SX(5000),SY(5000),SYER(5000)
      COMMON/TEXT1/TEXT,STEXT
      COMMON/TIN/IO,VARIN(10)

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

      CHARACTER*40 TEXT,STEXT*20

c first element in VARIN is spectra id of interest
      isp=INT(Varin(1))
      IF (VARIN(2).LT.VARIN(3)) GOTO 60
      IF (VARIN(2).GT.VARIN(3)) WRITE(IO,20)
20    FORMAT(1X,'Illegal limits (Min => Max)'/)
      WRITE(IO,40) XMIN(isp),XMAX(isp),XSCALE(isp,1),xscale(isp,2)
40    FORMAT(1X,'Current data limits : Xmin = ',G12.4,' Xmax = ',G12.4
     +  /1X,'Plot limits',9X,': Xmin = ',G12.4,' Xmax = ',G12.4/)
      RETURN
60    XSCALE(isp,1)=VARIN(2)
      XSCALE(isp,2)=VARIN(3)
C     Get array limits for use with plotting.
C     Get total number of points.
      NPLT = 0
C     If no data readin jump this.
      GOTO (160) NSPNT + 1
      DO 80 I=1,NSPNT
      KXMIN(isp) = I
C     When lower limit found jump out of loop.
      IF (XSCALE(isp,1).LE.SXsp(isp,I)) GOTO 100
80    CONTINUE
100   DO 120 I=KXMIN(isp),NSPNT
      KXMAX(isp) = I
C     When upper limit found jump out of loop.
      IF (XSCALE(isp,2).LT.SXsp(isp,I)) GOTO 140
120   CONTINUE
      KXMAX(isp) = NSPNT + 1
140   NPLTsp(isp) = KXMAX(isp) - KXMIN(isp)
      KXMAX(isp) = KXMAX(isp) - 1
C     Write out number of points.
160   WRITE(IO,180) NPLTsp(isp)
180   FORMAT(1X,'There are ',I4,' points in current plotting range'/)
      RETURN
C     *******
      END


      SUBROUTINE FTYSET
c rewritten 2/13/96
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/TIN/IO,VARIN(10)
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

      isp=INT(varin(1))
      IF (VARIN(2).LT.VARIN(3)) GOTO 60
      IF (VARIN(2).GT.VARIN(3)) WRITE(IO,20)
20    FORMAT(1X,'Illegal limits (Min => Max)'/)
      WRITE(IO,40) YMIN(isp),YMAX(isp),YSCALE(isp,1),yscale(isp,2)
40    FORMAT(1X,'Current data limits : Ymin = ',G12.4,' Ymax = ',G12.4
     +  /1X,'Plot limits',9X,': Ymin = ',G12.4,' Ymax = ',G12.4/)
      RETURN
60    YSCALE(isp,1)=VARIN(2)
      YSCALE(isp,2)=VARIN(3)
      RETURN
C     *******
      END

c this subroutine modified by KWH 2/21/96
      SUBROUTINE FTLIST(CALSUB,ESTERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      EXTERNAL CALSUB
      COMMON/FITLIM/XLIM(6),NSET,NLIM(6),NFIT,NPLT
      COMMON/SPLDAT/NSPNT,SX(5000),SY(5000),SYER(5000)
      COMMON/RESDAT/RY(5000),RYER(5000),WX(5000),RI(5000)
      COMMON/RESLIM/RTOF,NEP,LXMIN,LXMAX,ERRSUM
      COMMON/TEXT1/TEXT,STEXT
      COMMON/TITLES/PNAMES(40),TX(4),TY(4)
      COMMON/PARTIT/NPARAS
      COMMON/CONFIT/DMAX,ACC,H,MAXNO
      COMMON/FITPAR/PAR(500),STEP(500),OLD(500),INDEX(500)
      COMMON/CALDAT/CALPAR(100),YCALC(5000),WT(5000)
      COMMON/FILES/RFILE,LENRF,SFILE,LENSF,SPNUM
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

      real*8 fi(5000)

      CHARACTER*5 DD*9,TEXT*40,STEXT*20,TX,TY,TT*8
      CHARACTER*20 RFILE,SFILE,PNAMES*15,FILNAM*24
      INTEGER SPNUM
      LOGICAL NEWPRT
      REAL*8 ESTERR(500)
      real*8 Momtra,MolRad

C Subroutine writes out current parameters and value of the residual
      LPT = 8
C Get date and time.
c      CALL DATE(DD)
c      CALL TIME(TT)

c always open a new file for this output
c no longer sequentially updating individual spectra

      LENFIL=LENSF+4
      FILNAM(1:LENFIL)=SFILE(1:LENSF)//'.PRT'
      OPEN(LPT,FILE=FILNAM(1:LENFIL),STATUS='NEW')

C Write out program name, date, time, and sample info.
      WRITE(LPT,20) RFILE(1:LENRF),SFILE(1:LENSF)
20    FORMAT('1',5X,'Output by program MULSPEC ',
     +   //,6X,'Resolution filename =  ',A,
     +    /,6X,'Sample filename = ',A)
c
c need to loop through individual spectra
c slight resequencing of results from earlier versions
c go ahead and write out only spectra that are being fit
      do 999 ij=1,fitspc
      isp=ispec(ij)
      WRITE (LPT,40) isp
40    FORMAT(6X,'Spectrum number = ',I2)
      WRITE(LPT,60) TY,TX

60    FORMAT(/,6X,'Parameters in fitting Y : ',4A5/16X,
     +            'Versus X : ',4A5)
c ok loop over this spectrum's parameters
      write(LPT,79)
79    FORMAT(6X,'Number',5X,'Name',6X,'Value',6X,'( Old value  )',
     +  3X,'Step',9X,'Esd')
      do 85 ip=1,npspec
        j=(ij-1)*npspec+ip
        WRITE(LPT,80) J,PNAMES(ip),PAR(J),OLD(J),
     +    STEP(J),ESTERR(J)
80      FORMAT(6X,I4,3X,A15,2X,G12.4,' (',G12.4,')',
     +     G12.4,2X,G12.4)
85    continue
      WRITE(LPT,120) XMIN(isp),XMAX(isp),YMIN(isp),YMAX(isp)
     +  ,XSCALE(isp,1),xscale(isp,2),YSCALE(isp,1),yscale(isp,2)
120   FORMAT(/,6X,'Current data limits : Xmin = ',G12.4,' Xmax = ',
     +  G12.4,/,29X,'Ymin = ',G12.4,' Ymax = ',G12.4
     +  /,6X,'Plot limits',9X,': Xmin = ',G12.4,' Xmax = ',G12.4
     +  /29X,'Ymin = ',G12.4,' Ymax = ',G12.4)

	WRITE(LPT,140) NFITsp(isp)
140   FORMAT(/6X,'There are ',I5,' points in the current',
     +  1X,'range set by the "O" command')
      I2 = NSETsp(isp) * 2
      WRITE(LPT,160) (XLIMsp(isp,I-1),XLIMsp(isp,I),I=2,I2,2)
160   FORMAT(6X,'Current limits are :',
     +  1X,G12.4,'  to  ',G12.4,2(/22X,G12.4,'  to  ',G12.4)/)
      WRITE (LPT,180) NRPTOT(isp),RXMINsp(isp),RXMAXsp(isp)
180   FORMAT(/,6X,'Resolution function has ',I5,' points',
     1 1X,'ranging from ',F7.2,' to ',F6.2,' meV')
c end loop over fitted spectra
999   continue
      WRITE(LPT,100) DMAX,ACC,H
100   FORMAT(6X,' 100',3X,'Maximum step',G12.4
     +  /,6X,' 200',3X,'Accuracy    ',G12.4
     +  /,6X,' 300',3X,'Step width  ',G12.4)
C If insufficient points to define fit then return.
c must calculate resids over entire number of fitted spectra
      RESID = 0.0
c loop over all spectra
      do 9200 iij=1,fitspc
c select the fitted spectra
      isp=ispec(iij)
      DO 700 I=1,NPspec
c after fit transfer a copy of current fit values into global arrays
c so that if number of spectra is changed after this
c results are not lost
      parsp((isp-1)*npspec+i)=PAR(I+(iij-1)*npspec)
      stepsp((isp-1)*npspec+i)=step(I+(iij-1)*npspec)
c initialize paramters for CALSUB for this spectra
700   CALPAR(I) = PAR(I+(iij-1)*npspec)

c ok now we need to create all the arrays and elements 
c necessary for calculating the fit etc.
c also need to create MomTra
      Momtra=qspec(isp)
      lxmin=lxminsp(isp)
      lxmax=lxmaxsp(isp)
      nep=nepsp(isp)
      errsum=esumsp(isp)

c create arrays needed 

      do 710 i=1,nspnt
        sx(i)=sxsp(isp,i)
        sy(i)=sysp(isp,i)
        syer(i)=syersp(isp,i)
        ry(i)=rysp(isp,i)
        ryer(i)=ryersp(isp,i)
        wx(i)=wxsp(isp,i)
        ri(i)=risp(isp,i)
710   continue

c ok now call CALSUB just like normal

      CALL CALSUB(KXMIN(isp),KXMAX(isp))

c calculate the incremental to the residuals
c as well as the array used in SPLT

      DO 720 I=1,5000
      FIWTsp(isp,I)=0
720   CONTINUE
      R1=NLIMsp(isp,1)
      R2=NLIMsp(isp,NSETsp(isp)*2)
      DO 760 J=1,NSETsp(isp)
      J2 = J * 2 - 1
      DO 740 I=NLIMsp(isp,J2),NLIMsp(isp,J2+1)
      FI(I) = (SY(I) - YCALC(I))
      FIWTsp(isp,I) = FI(I) / WT(I)
740   RESID = RESID + FIWTsp(isp,I) * FIWTsp(isp,I)
760   CONTINUE
c end loop over all spectra
9200  CONTINUE
      IPAR = 0
c
c initialize nfit
      nfit=0
      do 541 i=1,fitspc
      isp=ispec(i)
      nfit=nfit+nfitsp(isp)
541   continue
c
      DO 560 I=1,NPARAS
560   IF (STEP(I).GT.0.0) IPAR = IPAR + 1

C     Divide residuals by N - P + C
      RESID = RESID / (NFIT - IPAR)

C Write out to LPT file.
      WRITE(LPT,280) RESID
280   FORMAT(/,6X,'Value of the residuals is ',G12.4)
c close lpt file
      close(unit=lpt,status='save')
300   RETURN
C     *******
      END


      SUBROUTINE FTERAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******


c this routine is unnecessary in the version running on windows


C     Erases screen and leaves terminal in ALPHA mode with cursor
C     in TLH corner.
c     IO = 6
c     ICOM = 33
C      WRITE(IO,1) (ICOM,I=1,4)
C    1 FORMAT(1X,A1,'WOR 33 H'/1X,A1,'GRA 1,35'/1X,A1,'SHR B'/1X,A1,'UP 35')
c      WRITE (6,*) CHAR(27),'[2J'
      RETURN
C     *******
      END
      SUBROUTINE FTRSET
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******

c this routine is unneceesary in windows mode
C     Returns to text mode.
c      IO = 6
c      ICOM = 33
C      WRITE(IO,1) (ICOM,I=1,2)
C    1 FORMAT(1X,A1,'WOR 0'/1X,A1,'MON H K')
C      WRITE (6,*) CHAR(2)
      RETURN
C     *******
      END
      SUBROUTINE FTCMD(IVAL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *******
      COMMON/TIN/IO,VARIN(10)
      DIMENSION ICOM(20)
      CHARACTER*75 IBUF
	character*1 cchar

C     Command options in ICOM.
      DATA ICOM/'H','A','X','Y','O','V','R','P','F','B','Q','L',
     +  'S','C','D','E','G','W',2*' '/
C     Define number of commands.
      IIN=5
      NCOM = 18
C     Give list of options.
20    WRITE(IO,40) (ICOM(I),I=1,NCOM)
40    FORMAT(1X,'Type help or option : ',A1,18(',',A1))
      WRITE(IO,60)
60    FORMAT(1X,'->  ',$)
      DO 80 I=1,10
80    VARIN(I)=0.0
C     Read in next command string.
c rewrite so as to decode assuming spaces are the delimeter
	do 81 i=1,75
		ibuf(i:i)=' '
81	continue
	read(iin,100,err=20,end=20)ibuf
100	format(a)
c first character in the string is the command itself
	cchar(1:1) = ibuf(1:1)
	if ((cchar(1:1) .eq. 'H') .or. (cchar(1:1) .eq. 'h')) then
		ival=1
	else if ((cchar(1:1) .eq. 'A') .or. (cchar(1:1) .eq. 'a')) then
		ival=2
	else if ((cchar(1:1) .eq. 'X') .or. (cchar(1:1) .eq. 'x')) then
		ival=3
	else if ((cchar(1:1) .eq. 'Y') .or. (cchar(1:1) .eq. 'y')) then
		ival=4
	else if ((cchar(1:1) .eq. 'O') .or. (cchar(1:1) .eq. 'o')) then
		ival=5
	else if ((cchar(1:1) .eq. 'V') .or. (cchar(1:1) .eq. 'v')) then
		ival=6
	else if ((cchar(1:1) .eq. 'R') .or. (cchar(1:1) .eq. 'r')) then
		ival=7
	else if ((cchar(1:1) .eq. 'P') .or. (cchar(1:1) .eq. 'p')) then
		ival=8
	else if ((cchar(1:1) .eq. 'F') .or. (cchar(1:1) .eq. 'f')) then
		ival=9
	else if ((cchar(1:1) .eq. 'B') .or. (cchar(1:1) .eq. 'b')) then
		ival=10
	else if ((cchar(1:1) .eq. 'Q') .or. (cchar(1:1) .eq. 'q')) then
		ival=11
	else if ((cchar(1:1) .eq. 'L') .or. (cchar(1:1) .eq. 'l')) then
		ival=12
	else if ((cchar(1:1) .eq. 'S') .or. (cchar(1:1) .eq. 's')) then
		ival=13
	else if ((cchar(1:1) .eq. 'C') .or. (cchar(1:1) .eq. 'c')) then
		ival=14
	else if ((cchar(1:1) .eq. 'D') .or. (cchar(1:1) .eq. 'd')) then
		ival=15
	else if ((cchar(1:1) .eq. 'E') .or. (cchar(1:1) .eq. 'e')) then
		ival=16
	else if ((cchar(1:1) .eq. 'G') .or. (cchar(1:1) .eq. 'g')) then
		ival=17
	else if ((cchar(1:1) .eq. 'W') .or. (cchar(1:1) .eq. 'w')) then
		ival=18
	else
		write(io,*)'invalid entry - try again'
		goto 20
	end if

c first character in string is the command itself
c how long is ibuf?
	i=75
110	i=i-1
	if (ibuf(i:i) .eq. ' ') goto 110
	lenbuf=i
c now parse for the number of entries
	ip=2
	i=0
120	i=i+1
c if no commands
140	if (ip .gt. lenbuf) goto 160
	if (ibuf(ip:ip) .eq. ' ') then
		ip=ip+1
		goto 140
	else
c ok at least one command read it
		read(ibuf(ip:lenbuf),*)varin(i)
c find the next blank and start over
150		ip=ip+1
		if (ip .gt. lenbuf) goto 160
		if (ibuf (ip:ip) .ne. ' ') goto 150
		goto 120
	end if
160	continue	
      RETURN
C     *******
      END

c this routine modified by KWH for multi-spectra analysis

      subroutine ftffff(calsub,maxfun,esterr)
      implicit double precision(a-h,o-z)

      external lsfun1

      common/fitpar/par(500),step(500),old(500),index(500)
      common/fitlim/xlim(6),nset,nlim(6),nfit,nplt
      common/titles/pnames(40),tx(4),ty(4)
      common/partit/nparas
      common/confit/dmax,acc,h,maxno
      common/spldat/nspnt,sx(5000),sy(5000),syer(5000)
      common/text1/text,stext
      common/caldat/calpar(100),ycalc(5000),wt(5000)
      common/tin/io,varin(10)

      character*5 tx,ty,text*40,stext*20
      character*15 pnames
      real*8 espar(100),esterr(500)

      real*8 hesian(100,100),cj(100),resids(25000)
      real*8 ftol,xtol,gtol,epsfcn,factor
      real*8 diag(100),fjac(25000,100),qtf(100),
     +       wa1(100),wa2(100),wa3(100),wa4(25000)
      integer ipvt(100),maxfev,mode,nprint,info,nfev,ldfjac
      real*8 rcond,det
      integer kpvt(100),job,inert

      character infstr(9)*120

C-----------------------------------------------------------------------
C Put the messages for the lmdif info parameter into infstr
C-----------------------------------------------------------------------
      infstr(1)= 'info = 0  improper input parameters.'
      infstr(2)= 'info = 1  both actual and predicted relative'//
     +        ' reductions in the sum of squares are at most ftol.'
      infstr(3)= 'info = 2  relative error between two consecutive'//
     +        ' iterates is at most xtol.'
      infstr(4)= 'info = 3  conditions for info = 1 and info = 2'//
     +        ' both hold.'
      infstr(5)= 'info = 4  the cosine of the angle between fvec and'//
     +        ' any column of the jacobian is at most gtol in'//
     +        ' absolute value.'
      infstr(6)= 'info = 5  number of calls to fcn has reached or'//
     +        ' exceeded maxfev.'
      infstr(7)= 'info = 6  ftol is too small. no further reduction'//
     +        ' in the sum of squares is possible.'
      infstr(8)= 'info = 7  xtol is too small. no further improvement'//
     +        ' in the approximate solution x is possible.'
      infstr(9)= 'info = 8  gtol is too small. fvec is orthogonal to'//
     +        ' the columns of the jacobian to machine precision.'



c this routine looks just like the original version
C-----------------------------------------------------------------------
C Determine number of floating parameters.
C-----------------------------------------------------------------------
      nvar = 0
      do 20 i=1,nparas
C-----------------------------------------------------------------------
C Keep an INDEX array for use with constraints.
C The method chosen allows any parameter to be tied to any parameter
C whether or not it is floated.
C-----------------------------------------------------------------------
      index(i) = 0
      if (abs(step(i)).lt.1e-4) goto 20
      if (step(i).lt.0.0) then
        index(i) = index(abs(step(i)))
      else
        nvar = nvar + 1
        index(i) = nvar
C-----------------------------------------------------------------------
C Use step size as initial estimate.
C-----------------------------------------------------------------------
        espar(nvar) = 1.0 / step(i)
      endif
20    continue

C-----------------------------------------------------------------------
C go to least squares routine lmdif from minpack-1
C-----------------------------------------------------------------------
      ftol=1e-8
c       ftol=1e-6
      xtol=1e-6
      gtol=1e-4
      maxfev=1500
      epsfcn=1e-8
      mode=1
      factor=100.0
      nprint=1
      ldfjac=25000
      call lmdif(lsfun1,nfit,nvar,espar,resids,ftol,xtol,gtol,maxfev,
     +   epsfcn,diag,mode,factor,nprint,info,nfev,fjac,ldfjac,ipvt,
     +   qtf,wa1,wa2,wa3,wa4)
      write (6,695) infstr(info+1),nfev
695   format(1x,A,
     +     /,1x,'Number of function calls = ',i3)

C-----------------------------------------------------------------------
C Put best parameters back into PAR.
C-----------------------------------------------------------------------
      do 120 i=1,nparas
      stp = step(i)
      if (stp.lt.0.0) stp = step(int(abs(stp)))
      if (stp.lt.1e-4) goto 120
      par(i) = par(i) * espar(index(i)) * stp
120   continue


C-----------------------------------------------------------------------
c get jacobian matrix
C-----------------------------------------------------------------------
      call fdjac2(lsfun1,nfit,nvar,espar,resids,fjac,ldfjac,
     +                  iflag,epsfcn,wa4)
C-----------------------------------------------------------------------
c Use this to create Hessian
C-----------------------------------------------------------------------
      do 1800 i=1,nvar
      do 1810 j=i,nvar                                        
      hesian(i,j)=0.0
      do 1820 k=1,nfit
      hesian(i,j)=hesian(i,j)+fjac(k,i)*fjac(k,j)
1820  continue
1810  continue
1800  continue
C-----------------------------------------------------------------------
C now take the inverse to get the covariance matrix (using Linpack)
C-----------------------------------------------------------------------
      call dsico(hesian,100,nvar,kpvt,rcond,wa1)
      if ((rcond+1.0).eq.1.0) then
      write (6,1825)
1825  format(1x,'--------Hessian matrix is singular-------')
      goto 100
      else
      job=1
      call dsidi(hesian,100,nvar,kpvt,det,inert,wa1,job)
C-----------------------------------------------------------------------
C As the matrix is symmetric, need to copy elements to lower half
C-----------------------------------------------------------------------
      do 1827 i=1,nvar
      do 1829 j=i+1,nvar
      hesian(j,i)=hesian(i,j)
1829  continue
1827  continue
      endif
C-----------------------------------------------------------------------
C Caculate Chi-squared(v) for v (nvar) variables from the incomplete
C Gamma function for 95% confidence limit
C-----------------------------------------------------------------------
      chi2v=gchi2v(nvar,0.68D0)
C-----------------------------------------------------------------------
C Calculate the estimated standard deviations (1 sigma = 68 %)
C-----------------------------------------------------------------------
      do 1831 i=1,nvar
      cj(i)=chi2v*hesian(i,i)
1831  continue

      do 1901 i=1,nparas
      stp = step(i)
      if (stp.lt.1e-4) then
      esterr(i)=0.0
      goto 1901
      endif
      if (stp.lt.0.0) then
        stp = step(int(abs(stp)))
      endif
      esterr(i) = par(i) * sqrt( cj(index(i)) )
1901  continue
100   continue
      end

c this subroutine has seen some changes by KWH
c it is the routine which returns the vector of residuals 
c back to the fitting routines
c this routine must keep track of which spectra are being fit
c the parameter order in PAR is in order of 
c fitted spectra
c modified 2/21/96 

      subroutine lsfun1(m,n,xc,fvecc,iflag)
      implicit double precision(a-h,o-z)
      common/spldat/nspnt,sx(5000),sy(5000),syer(5000)
      COMMON/RESDAT/RY(5000),RYER(5000),WX(5000),RI(5000)
      COMMON/RESLIM/RTOF,NEP,LXMIN,LXMAX,ERRSUM
      common/caldat/calpar(100),ycalc(5000),wt(5000)
      common/fitpar/par(500),step(500),old(500),index(500)
      common/fitlim/xlim(6),nset,nlim(6),nfit,nplt
      common/partit/nparas
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

      integer iflag
      real*8 fvecc(m),xc(n)
      real*8 MOMTRA,MolRad

C-----------------------------------------------------------------------
C Transfer parameters to array used in funct3
C-----------------------------------------------------------------------
c transfer parameters from par which contains the varied parameters
c to each varied spectra's parameter list individually, and individually 
c calculate the residuals.

c these initialization statements had to be moved outside the
c loop over fitted spectra
      k=0
      resid=0.0
      
      do 1999 jspec=1,fitspc
         isp=ispec(jspec)
C-----------------------------------------------------------------------
C Transfer parameters to array used in funct3
C-----------------------------------------------------------------------
      do 10 iij=1,npspec
c this little modification ensures that you are looking at the correct
c element in the PAR array for spectra isp
      i=(jspec-1)*npspec+iij
      stp = step(i)
      if (stp.lt.0.0) stp = step(int(abs(stp)))
C-----------------------------------------------------------------------
C copy over parameters without xc when stp=0.0
C-----------------------------------------------------------------------
      if (stp.gt.1e-4) then
        calpar(iij) = par(i) * xc(index(i)) * stp
        else
        calpar(iij) = par(i)
      endif
10    continue


c ok each entry into calpar is now correct for spectra ISP
c now initialize the sy,syer, and sx arrays for spectra isp
c inintialize all arrays required for fitting this spectra
c all spectra have the same number of points
c ok now we need to create all the arrays and elements 
c necessary for calculating the fit etc.

      lxmin=lxminsp(isp)
      lxmax=lxmaxsp(isp)
      nep=nepsp(isp)
      errsum=esumsp(isp)
c also need to initialize MOMTRA
      MOMTRA=qspec(isp)

c create arrays needed 

      do 710 i=1,nspnt
        sx(i)=sxsp(isp,i)
        sy(i)=sysp(isp,i)
        syer(i)=syersp(isp,i)
        ry(i)=rysp(isp,i)
        ryer(i)=ryersp(isp,i)
        wx(i)=wxsp(isp,i)
        ri(i)=risp(isp,i)
710   continue

c modifications for particular spectra being examined
      do 100 j=1,nsetsp(isp)
      j2 = j * 2 - 1
      nx1 = nlimsp(isp,j2)
      nx2 = nlimsp(isp,j2+1)
      call funct3(nx1,nx2)
      do 80 i=nx1,nx2
      k=k+1
c note that wt(i) is created in funct3 and
c is nothing more than syer(i)
      fvecc(k) = (sy(i) - ycalc(i)) / wt(i)
      resid=resid+fvecc(k)*fvecc(k)
c	if (ycalc(i) .gt. 1000) then
c	write(6,*) isp,resid,sy(i),ycalc(i),wt(i)
c	end if
80    continue
100   continue

c here is the termination for loop over fitted spectra
1999  continue

      if (iflag.eq.0) then
      write (6,120) resid/m
120   format(1x,'Residuals = ',e12.4e2)

c write out the varied parameter list to try and figure out why this blows up for 6th call at
c error estimation
c	do 2999 jspec=1,fitspc
c	do 999 iij=1,npspec
c           i=(jspec-1)*npspec+iij
c	   write(6,*) jspec,par(i)
c999	continue
c2999    continue

      endif
      end


      double precision function erfc(x)
C-----------------------------------------------------------------------
C Calculates the complementary error function
C-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter(p=0.3275911,a1=0.254829592,a2=-0.284496736,
     +       a3=1.421413741,a4=-1.453152027,a5=1.061405429)

      real*8 t,x

      t= 1.0 / ( 1.0 + p*abs(x) )
      erfc=((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-x*x)

C-----------------------------------------------------------------------
C Make the correction for negative x
C-----------------------------------------------------------------------
      if (x.lt.0) erfc=2-erfc

      end

      double precision function gchi2v(nvar,conlim)
      implicit real*8(a-h,o-z)
      real*8 chicon,conlim,a,x1,x2,dx
      integer nvar

      a=nvar/2.0
C start search at delta chi2=2.0 (x=1) and then double if not in interval
      x1=0.0
      x2=1.0
10    chicon=gammaq(a,x2)
      if (chicon.lt.conlim) then
      x1=x2
      x2=2*x2
      goto 10
      endif
      do 20 i=1,100
      dx=0.5*(x2-x1)
      x1=x1+dx
      chicon=gammaq(a,x1)
      if (abs(chicon-conlim).lt.1.0e-5) goto 30
      if (chicon.gt.conlim) then
      x2=x2-dx
      x1=x1-dx
      endif
20    continue
      write (6,40) conlim,chicon-conlim
40    format(1x,'Root differs from requested value ',f6.3,' by ',f7.4)
30    continue
      gchi2v=2.0*x1
      end

      double precision function gammaq(a,x)
      implicit real*8(a-h,o-z)
      real*8 a,x,eps,ap,del,sum,gln
      real*8 gold,g,fac,b1,b0,anf,ana,an,a1,a0

      gold=0.0
      fac=1.0
      b1=1.0
      b0=0.0
      a0=1.0
      eps=3e-7

      if (x.lt.(a+1.0)) then
        gln=gammln(a)
        ap=a
        del=1.0/a
        sum=del
        do 10 i=1,100
        ap=ap+1.0
        del=del*x/ap
        sum=sum+del
         if (abs(del).lt.(abs(sum)*eps)) then
         gammaq=sum*exp(-x+a*log(x)-gln)
         goto 100
         endif
10      continue
        write (6,20)
20      format(1x,'Not enough iterations in gammaq for convergence')
      else
        gln=gammln(a)
        a1=x
        do 30 i=1,100
        an=i
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if (a1.ne.0) then
          fac=1.0/a1
          g=b1*fac
          if (abs((g-gold)/g).lt.eps) then
            gammaq=1.0-exp(-x+a*log(x)-gln)*g
          goto 100
          endif
          gold=g
        endif
30      continue
        write (6,20)
      endif
100   continue
      end

      double precision function gammln(xx)
      implicit real*8(a-h,o-z)
      real*8 xx,x,tmp,ser,cof(6)
      integer j

      data cof/76.18009173d0,-86.50532033d0,24.01409822d0,
     +      -1.231739516d0,0.120858003d-2,-0.536382d-5/

      x=xx-1.0
      tmp=x+5.5
      tmp=tmp-(x+0.5)*log(tmp)
      ser=1.0
      do 10 j=1,6
      x=x+1.0
      ser=ser+cof(j)/x
10    continue
      gammln=-tmp+log(2.50662827465d0*ser)
      end
