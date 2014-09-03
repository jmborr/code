C this subroutine written by KWH 11/95
c its purpose is to list the spectra available to be fit
c and also asks the user how many and which spectra should
c be fit. it modifies a common which contains the number
c of spectra to be fit and an integer array containing the spectra
c that are to be fit

      SUBROUTINE LSPEC
      
      character*1 answ

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

c the routine also has access via a common statement to the number of
c available spectra and their Q-settings

c list number of available spectra
      write(6,100)totspc
100   format(1x,'the number of spectra available to be fit is:',i3)
      write(6,110)
110   format(1x,'more complete information on spectra (Y/N)?',$)
      read(5,121)answ
121   format(a)
      if ((answ .eq. 'Y').or.(answ .eq. 'y')) then
        write(6,115)
115     format(1x,' spectra id     Q (A**-1)')
        do 150 i=1,totspc
           write(6,*)i,qspec(i)
120        format(4x,i2,8x,f6.3)
150     continue
      endif

      write(6,160)fitspc
160   format(1x,'currently fitting ',i2,' spectra')
      write(6,161)
161   format(1x,'change the list of spectra being fit?',$)
      read(5,121)answ
      if ((answ .eq. 'y').or.(answ .eq. 'Y')) then
165      write(6,170)
170      format(1x,'enter number of spectra to be fit:',$)
         read(5,*)fitspc
         if ((fitspc .gt. totspc) .or. (fitspc .le. 0)) goto 165
         do 200 i=1,fitspc
            write(6,175)i
175         format(1x,'enter id number of',i2,'th spectra to be fit:',$)
            read(5,*)ispec(i)
200      continue
      endif

      return
      end
