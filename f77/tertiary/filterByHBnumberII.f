*     pgf77 -g -Mbounds -Mchkfpstk -Mextend  -Wl,-static -o filterByHBnumberII.x filterByHBnumberII.f
*     To debug, type in fec02: pgdbg filterByHBnumber.x
*     pgf77 -Mextend -O -s -fast -Wl,-static -o filterByHBnumberII.x filterByHBnumberII.f
*
*     This program filters a TASSER trajectory by number of hydrogen bonds

      program filterByHBnumber
      implicit none
      integer nmaxres
      parameter(nmaxres=2000) !maximum number of residues
      logical readcommand

      character*512 train     !input trajectory file
      character*512 traout     !output trajectory file
      character*512 secf      !seq.dat file
      character*512 header
      character*512 snapshot(nmaxres)
      integer avhbn           !number of hydrogen bonds averaged over the trajectory
      real rav                !running average of average number of HB's
      real newrav
      integer nacc            !number of accepted snapshots
      integer maxacc
      real xyz(3,nmaxres)     !coordinates array
      integer lch               !sequence length
      integer sec(nmaxres)    !secondary structure from seq.dat
      integer trainp
      integer traoutp
      real bondvec(3,nmaxres),bonduni(3,nmaxres)
      real bisector(3,nmaxres), hbvec(3,nmaxres)
      real hbmap(nmaxres,nmaxres)
      integer nb,nexcess
      integer getnhb
      real getehb
      real ehb !energy of the conformation
      real eexcess !energy of the conformation allowing more than two bonds per residue
c     manage input command arguments
      if(readcommand(train,traout,secf,avhbn).eq..true.) call exit(1)
c     read input secondary structure assignments
      call readSeqDatFile(secf,lch,sec)
c     cycle over the trajectory file
      rav=-1.0
      nacc=0
      trainp=1 !file unit for input trajectory file
      traoutp=2
      open(unit=trainp,file=train,status='old',err=666)
      open(unit=traoutp,file=traout,status='unknown') !will overwrite
      maxacc=100000
 10   call readsnapshot(trainp,xyz,header,snapshot,lch)
      if(lch.eq.0) goto 20        !we finished reading the trajectory file
c     prepare bond, bisector, and hydrogen-bond vectors
      call prepare_vectors(xyz,bondvec,bonduni,bisector,hbvec,lch)
c     calculate hydrogen bond contact map
      call hydrogenbondmap(xyz,bondvec,bonduni,bisector,hbvec,lch,
     $    sec,hbmap,nexcess,eexcess)
      nb=getnhb(hbmap,lch) !calculate number of hydrogen bonds
      ehb=getehb(hbmap,lch)
      newrav=(rav*nacc+nb)/(nacc+1)
      if( (newrav.gt.rav) .or. (newrav.gt.avHBn) )then
         call savesnapshot(traoutp,header,snapshot,lch)
         rav=newrav  !update running average
         nacc=nacc+1 !update number of acceptede snapshots
      endif
      if(nacc.gt.maxacc) goto 20
      write(6,'2(f7.1x)')eexcess,ehb
      goto 10 !read next snapshot
 20   close(2)
      close(1)
c      write(6,*)nacc,rav
      return
 666  call exit(1)
      end !end program filterByHBnumber
c     ################################################


c     ################################################
      function readcommand(train,traout,secf,avHBn)
c     read command line arguments
c     input variables
      logical readcommand
      character*512 train  !input trajectory file
      character*512 traout !output trajectory file
      character*512 secf   !seq.dat file
      integer avhbn     !average number of hydrogen bonds
c     local variables
      character*512 fnam
      integer ireq,nreq !nreq: number of Required arguments
      integer narg      !number of words in the command line
      integer i
      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h'
     $     .or.fnam.eq.'-help')then
 100     call system('clear')
         write(*,*)'Usage: filterByHBnumber.x [options]:'
         write(*,*)' Required arguments:'
         write(*,*)' -train input trajectory file'
         write(*,*)' -traout output trajectory file'
         write(*,*)' Optional arguments:'
         write(*,*)' -avHBn average number of hydrogen bonds in',
     $        'trajectory' 
         readcommand=.true. !signal error
         return
      endif

      write(secf,*) 'seq.dat'  !fix this later

      i=0
      nreq=2
      ireq=0
      narg=iargc()
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-train')then
         i=i+1
         call getarg(i,fnam)
         write(train,*) fnam(1:lnblnk(fnam))
         ireq=ireq+1
      elseif(fnam.eq.'-traout')then
         i=i+1
         call getarg(i,fnam)
         write(traout,*) fnam(1:lnblnk(fnam))
         ireq=ireq+1
      elseif(fnam.eq.'-secf')then
         i=i+1
         call getarg(i,fnam)
         write(secf,*) fnam(1:lnblnk(fnam))
      elseif(fnam.eq.'-avHBn')then
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)avhbn
c         write(avhbn,*) fnam(1:lnblnk(fnam))
      endif
      if(i.lt.narg) goto 115
      if(ireq.lt.nreq) goto 100 !print welcome message and output program

      end !end of subroutine readcommand
c     ################################################


c     ###################################################
      subroutine readSeqDatFile(secf,lch,sec)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
c     input variables
      character*512 secf
      integer lch
      integer sec(nmaxres)
c     local variables
      integer i,j
      character*3 res
      open(unit=10,file=secf,status='old',err=666)
      lch=1
 10   read(10,*,end=20)i,res,sec(lch),j
      lch=lch+1
      goto 10
 20   lch=lch-1
      return
 666  call exit(1)
      end
c     ################################################


c     ################################################
      subroutine readsnapshot(trainp,xyz,header,snapshot,lch)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
c     input variables
      integer trainp
      real xyz(3,nmaxres)       !coordinates array
      character*512 header
      character*512 snapshot(nmaxres)
      integer lch
c     local variables
      character*512 line
      integer i,ilch
      read(trainp,'a512',err=666,end=666)header(1:512)
      write(lch,*)header !first number of header is sequence length
      do ilch=1,lch
         read(trainp,'a512')line(1:512)
         write(snapshot(ilch),'a512')line
         read(snapshot(ilch),*)xyz(1,ilch),xyz(2,ilch),xyz(3,ilch) 
      enddo
      do i=1,3
         xyz(i,lch)=xyz(i,lch)/0.87 !from Angstroms to lattice units
      enddo
      return
 666  lch=0 !signal error
      end !end of subroutine readsnapshot
c     ################################################


c     ###################################################
      subroutine prepare_vectors(xyz,bondvec,bonduni,
     $     bisector,hbvec,lch)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
c     input variables
      integer lch
      real xyz(3,nmaxres)
      real bondvec(3,nmaxres),bisector(3,nmaxres), hbvec(3,nmaxres)
      real bonduni(3,nmaxres)
c     local variables
      integer i,j,ilch
      real dd

c     calculate bond vectors and unit bond vectors
      do ilch=2,lch
         dd=0.0
         do j=1,3
            bondvec(j,ilch)=xyz(j,ilch)-xyz(j,ilch-1)
            dd=dd+bondvec(j,ilch)*bondvec(j,ilch)
         enddo
         do j=1,3
            bonduni(j,ilch)=bondvec(j,ilch)/sqrt(dd)
         enddo
      enddo

c     calculate bisector vectors
      do ilch=2,lch-1
         dd=0.0
         do j=1,3
            bisector(j,ilch)=bondvec(j,ilch+1)-bondvec(j,ilch)
            dd=dd+bisector(j,ilch)*bisector(j,ilch)
         enddo
         do j=1,3
            bisector(j,ilch)=bisector(j,ilch)/sqrt(dd) !normalize
         enddo
      enddo

c     calculate hydrogen bond vectors
      do ilch=2,lch-1
         dd=0.0
         hbvec(1,ilch)=bonduni(2,ilch)*bonduni(3,ilch+1)-
     $        bonduni(3,ilch)*bonduni(2,ilch+1)
         dd=dd+hbvec(1,ilch)*hbvec(1,ilch)
         hbvec(2,ilch)=bonduni(3,ilch)*bonduni(1,ilch+1)-
     $        bonduni(1,ilch)*bonduni(3,ilch+1)
         dd=dd+hbvec(2,ilch)*hbvec(2,ilch)
         hbvec(3,ilch)=bonduni(1,ilch)*bonduni(2,ilch+1)-
     $        bonduni(2,ilch)*bonduni(1,ilch+1)
         dd=dd+hbvec(3,ilch)*hbvec(3,ilch)
         do j=1,3
            hbvec(j,ilch)=hbvec(j,ilch)/sqrt(dd)
         enddo
      enddo

      end !end of subroutine prepare_vectors
c     ################################################


c     #################################################
      subroutine hydrogenbondmap(xyz,bondvec,bonduni,
     $     bisector,hbvec,lch,sec,hbmap,nexcess,eexcess)
c     not optimized, coded for easyness of reading.
c     input parameters
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
      integer longrange
      parameter(longrange=20)
      integer lch
      real xyz(3,nmaxres)
      real bondvec(3,nmaxres),bisector(3,nmaxres), hbvec(3,nmaxres)
      real bonduni(3,nmaxres)
      integer sec(nmaxres)
      real hbmap(nmaxres,nmaxres)
c     local variables
      integer i,j,ll
      real ehbij(nmaxres,nmaxres) !reinforcement for bond between similar secondary assignments
c     init hydrogen bond parameters
      real cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      real cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      character *512 la,lb
      integer ires,jres,kres,lres,idist
      real relvec(3) !relative vector
      real cr2,bb,cc,av11,av12,av21,av22,dotproduct,dotproductii
      real fact,energyHBa,energyHBb
      integer haccept(nmaxres),hdonate(nmaxres)
      real ahacc(nmaxres),ahdon(nmaxres)
      integer nacceptor,nexcess
      real nehb,eexcess
      nexcess=1
      eexcess=0.0
      write(la,*)'48 0.45 0.1 0.00'
      read(la,*)cr2a,acut_bb,acut_cc,acut_vv
      write(lb,*)'50 0.25 0.4 0.35'
      read(lb,*)cr2b,bcut_bb,bcut_cc,bcut_vv
      call set_ehb(sec,ehbij,lch)
      do ires=1,lch !reset the hydrogen bond network
         haccept(ires)=0
         hdonate(ires)=0
         do jres=1,lch
            hbmap(ires,jres)=0.0
         enddo
      enddo
      do ires=2,lch-2
         do jres=ires+1,lch-1
            nacceptor=0
            fact=-1.0
            cr2=0.0
            do i=1,3 !initialize relative vector
               relvec(i)=xyz(i,jres)-xyz(i,ires)
               cr2=cr2+relvec(i)*relvec(i)
            enddo
            idist=jres-ires
            bb=dotproduct(hbvec,ires,hbvec,jres)
            cc=dotproduct(bisector,ires,bisector,jres)
            av11=dotproduct(bondvec,ires,bondvec,jres)
            av12=dotproduct(bondvec,ires,bondvec,jres+1)
            av21=dotproduct(bondvec,ires+1,bondvec,jres)
            av22=dotproduct(bondvec,ires+1,bondvec,jres+1)
c     alpha-helix 
            if(idist.eq.3 .and.
     $           sec(ires).eq.2.and.sec(jres).eq.2 .and. !both not assigned as helixes
     $           cr2.lt.cr2a .and. !|r_i-r_k|^2 < Cr2a=36Angstroms^2
     $           cc.gt.acut_cc .and. !c_i*c_k < acut_cc=0.1, c_i:bisector vector
     $           bb.gt.acut_bb .and. !acut_bb=0.45
     $           av11.gt.acut_vv .and. !acut_vv=0.0
     $           av22.gt.acut_vv
     $           )then
               fact=(1-abs(cc-0.4))*(1-abs(bb-0.815)) !deviation from optimal
               nehb=energyHBa(ires,jres,fact,relvec,hbvec,ehbij)
               if(dotproductii(hbvec,ires,relvec).gt.0)then
                  nacceptor=ires
               else
                  nacceptor=jres
               endif

c     antiparallel-sheet, short range
            elseif(idist.gt.4.and.idist.lt.longrange .and.
     $              sec(ires).eq.4.and.sec(jres).eq.4 .and.
     $              cr2.lt.cr2b .and. !cr2b=38Angstroms^2
     $              cc.gt.bcut_cc .and. !bcut_cc=0.4
     $              bb.lt.-bcut_bb .and. !antiparallel, bcut_bb=0.25
     $              av12.lt.-bcut_vv .and. !bcut_vv=0.35
     $              av21.lt.-bcut_vv
     $              )then
               fact=abs(bb)*cc  !bb->-1,cc->1
               nehb=energyHBb(ires,jres,fact,relvec,hbvec,ehbij)
               if(ires.lt.jres)then
                  nacceptor=ires
               else
                  nacceptor=jres
               endif
c     (anti)parallel sheet, long range
            elseif( (idist.ge.longrange .and.
     $              sec(ires).eq.4.and.sec(jres).eq.4 .and.
     $              cr2.lt.cr2b.and.cc.gt.bcut_cc))then
               if(bb.lt.-bcut_bb .and. !antiparallel-sheet
     $              av12.lt.-bcut_vv .and.
     $              av21.lt.-bcut_vv)then
                  if(ires.lt.jres)then
                     nacceptor=ires
                  else
                     nacceptor=jres
                  endif
               elseif(bb.gt.bcut_bb .and. !parallel-sheet
     $                 av11.gt.bcut_vv .and.
     $                 av22.gt.bcut_vv)then
                  if(dotproductii(hbvec,ires,relvec).gt.0)then
                     nacceptor=ires
                  else
                     nacceptor=jres
                  endif
               endif
               fact=abs(bb)*cc       !bb->1,cc->1
               nehb=energyHBb(ires,jres,fact,relvec,hbvec,ehbij)
            endif

            if(nacceptor.gt.0)then !update new HB network, add energy
               nexcess=nexcess+1
               eexcess=eexcess+nehb
               if(nacceptor.eq.ires)then
                  kres=ires !kres allways acceptor
                  lres=jres
               else
                  kres=jres
                  lres=ires
               endif
               if(haccept(kres).gt.0)then !previous bond
                  if(ahacc(kres).lt.nehb) goto 10 !previous bond is more stable
                  ll=haccept(kres)
                  hdonate(ll)=0
               endif
               if(hdonate(lres).gt.0)then
                  if(ahdon(lres).lt.nehb) goto 10
                  ll=hdonate(lres)
                  haccept(ll)=0
               endif
               haccept(kres)=lres
               hdonate(lres)=kres
               ahacc(kres)=nehb
               ahdon(jres)=nehb
            endif
         enddo !jres-->jres+1
 10   enddo    !ires-->ires+1
c     fill hmbap
      do ires=1,lch
         if(haccept(ires).gt.0)then
            jres=haccept(ires)
            !write(6,*)ires,jres
            hbmap(ires,jres)=ahacc(ires)
            hbmap(jres,ires)=ahdon(jres)
            !write(6,*)hbmap(ires,jres),hbmap(jres,ires)
         endif
      enddo
      end !end of subroutine hbmap
c     ################################################


c     ###########################################################
      function dotproduct(vectorsi,i,vectorsj,j)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
c     input variables
      real vectorsi(3,nmaxres),vectorsj(3,nmaxres)
      integer i,j,k
c     local variables
      real dotproduct
      dotproduct=0.0
      do k=1,3
         dotproduct=dotproduct+vectorsi(k,i)*vectorsj(k,j)
      enddo
      end
c     ################################################

c     ###########################################################
      function dotproductii(vectorsi,i,vectorsj)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
c     input variables
      real vectorsi(3,nmaxres),vectorsj(3)
      integer i,j,k
c     local variables
      real dotproductii
      dotproductii=0.0
      do k=1,3
         dotproductii=dotproductii+vectorsi(k,i)*vectorsj(k)
      enddo
      end
c     ################################################


c     ##########################################################
      subroutine set_ehb(sec,ehbij,lch)
c     initialize ehbij matrix
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
      integer longrange
      parameter(longrange=20)
c     input variables
      integer sec(nmaxres)
      real ehbij(nmaxres,nmaxres)
      integer lch
c     local variables 
      integer ires,jres,is,js,idist
      do ires=1,lch
         is=sec(ires)
         do jres=1,lch
            js=sec(jres)
            ehbij(ires,jres)=0
            if(iabs(ires-jres).eq.3.and.is.eq.2.and.js.eq.2)then
               ehbij(ires,jres)=1.5 !hydrogen bond for helix enhanced
            elseif(is.eq.4.and.js.eq.4) then !both beta
               idist=iabs(ires-jres)
               if(idist.gt.longrange)then !longrange
                  ehbij(ires,jres)=3.0
               elseif(idist.gt.4)then
                  ehbij(ires,jres)=1.5
               endif
            endif
         enddo
      enddo
      end
c     ################################################


c     ##########################################################
      function energyHBa(ires,jres,fact,relvec,hbvec,ehbij)
      implicit none
      integer nmaxres      
      parameter(nmaxres=2000)
c     input variables
      integer ires,jres
      real fact, relvec(3),hbvec(3,nmaxres),ehbij(nmaxres,nmaxres)
c     local variables
      integer i,kres,lres,p
      real acut_hh,dif,ibr,jbr,f
      real energyHBa
      energyHBa=0.0
      acut_hh=5.0
      p=1
      if(ires.gt.jres) p=-1
      ibr=0.0
      jbr=0.0
      do i=1,3
         dif=5.7*hbvec(i,ires)-p*relvec(i) 
         ibr=ibr+dif*dif
         dif=5.7*hbvec(i,jres)-p*relvec(i)
         jbr=jbr+dif*dif
      enddo
      f=ehbij(ires,jres)*fact
      if(ibr.lt.acut_hh) energyHBa=energyHBa-f/(1+abs(ibr-3.2))
      if(jbr.lt.acut_hh) energyHBa=energyHBa-f/(1+abs(jbr-3.2))
      end
c     ################################################


c     ##########################################################
      function energyHBb(ires,jres,fact,relvec,hbvec,ehbij)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
c     input/output variables
      real energyHBa
      integer ires,jres
      real fact, relvec(3),hbvec(3,nmaxres),ehbij(nmaxres,nmaxres)
c     local variables
      integer i,kres,lres,p
      real bcut_hh,dif,ibrp,ibrm,jbrp,jbrm,f
      real energyHBb
      bcut_hh=20.0
      energyHBb=0.0
      ibrp=0.0
      jbrp=0.0
      ibrm=0.0
      jbrm=0.0
      do i=1,3
         dif=5.7*hbvec(i,ires)-relvec(i) 
         ibrm=ibrm+dif*dif
         dif=5.7*hbvec(i,ires)+relvec(i) 
         ibrp=ibrp+dif*dif
         dif=5.7*hbvec(i,jres)-relvec(i)
         jbrm=jbrm+dif*dif
         dif=5.7*hbvec(i,jres)+relvec(i) 
         jbrp=jbrp+dif*dif
      enddo
      f=ehbij(ires,jres)*fact
      if(ibrm.lt.bcut_hh) energyHBb=energyHBb-f/(2.+ibrm)
      if(ibrp.lt.bcut_hh) energyHBb=energyHBb-f/(2.+ibrp)
      if(jbrm.lt.bcut_hh) energyHBb=energyHBb-f/(2.+jbrm)
      if(jbrp.lt.bcut_hh) energyHBb=energyHBb-f/(2.+jbrp)
      end
c     ################################################


c     ################################################
      function getnhb(hbmap,lch)
      implicit none
      integer nmaxres
      real zero
      parameter(nmaxres=2000) !maximum number of residues
      parameter(zero=-0.00001)
c     input variables
      real hbmap(nmaxres,nmaxres)
      integer lch
c     local variables
      integer getnhb
      integer i,j
      getnhb=0
      do i=1,lch-1
         do j=i+1,lch
            if(hbmap(i,j).lt.zero) getnhb=getnhb+1
         enddo
      enddo
      end !end of function getnhb
c     ################################################


c     ################################################
      function getehb(hbmap,lch)
      implicit none
      integer nmaxres
      real zero
      parameter(nmaxres=2000) !maximum number of residues
      parameter(zero=-0.00001)
c     input variables
      real hbmap(nmaxres,nmaxres)
      integer lch
c     local variables
      real getehb
      integer i,j
      getehb=0.0
      do i=1,lch-1
         do j=i+1,lch
            if(hbmap(i,j).lt.zero) getehb=getehb+hbmap(i,j)
         enddo
      enddo      
      end !end of function getehb
c     ################################################


c     ################################################
      subroutine savesnapshot(traoutp,header,snapshot,lch)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000) !maximum number of residues
c     input variables
      integer traoutp
      character *512 header
      character *512 snapshot(nmaxres)
      integer lch
c     local variables
      integer lnblnk
      integer ilch
      character*512 line
      write(traoutp,'a')header(1:lnblnk(header))
      do ilch=1,lch
         write(line,'a') snapshot(ilch)
         write(traoutp,'a')line(1:lnblnk(line))
      enddo
      end !end of subroutine savesnapshot
c     ################################################
