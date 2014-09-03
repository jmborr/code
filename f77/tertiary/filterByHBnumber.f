*     pgf77 -g -Mbounds -Mchkfpstk -Mextend  -Wl,-static -o filterByHBnumber.x filterByHBnumber.f
*     To debug, type in fec02: pgdbg filterByHBnumber.x
*     pgf77 -Mextend -O -s -fast -Wl,-static -o filterByHBnumber.x filterByHBnumber.f
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
      integer L               !sequence length
      integer sec(nmaxres)    !secondary structure from seq.dat
      integer trainp
      integer traoutp
      real bondvec(3,nmaxres),bonduni(3,nmaxres)
      real bisector(3,nmaxres), hbvec(3,nmaxres)
      real hbmap(nmaxres,nmaxres)
      integer nb
      integer getnhb
c     manage input command arguments
      if(readcommand(train,traout,secf,avhbn).eq..true.) call exit(1)
c     read input secondary structure assignments
      call readSeqDatFile(secf,L,sec)
c     cycle over the trajectory file
      rav=-1.0
      nacc=0
      trainp=1 !file unit for input trajectory file
      traoutp=2
      open(unit=trainp,file=train,status='old',err=666)
      open(unit=traoutp,file=traout,status='old') !will overwrite
      maxacc=100000
 10   call readsnapshot(trainp,xyz,header,snapshot,L)
      if(L.eq.0) goto 20        !we finished reading the trajectory file
c     prepare bond, bisector, and hydrogen-bond vectors
      call prepare_vectors(xyz,bondvec,bonduni,bisector,hbvec,L)
c     calculate hydrogen bond contact map
      call hydrogenbondmap(xyz,bondvec,bonduni,bisector,hbvec,L,
     $    sec,hbmap)
      nb=getnhb(hbmap,L) !calculate number of hydrogen bonds
      newrav=(rav*nacc+nb)/(nacc+1)
      if( (newrav.gt.rav) .or. (newrav.gt.avHBn) )then
         call savesnapshot(traoutp,header,snapshot,L)
         rav=newrav  !update running average
         nacc=nacc+1 !update number of acceptede snapshots
      endif
      if(nacc.gt.maxacc) goto 20
      goto 10 !read next snapshot
 20   close(2)
      close(1)
      print *,'nacc=',nacc,rav
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
      subroutine readSeqDatFile(secf,L,sec)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
c     input variables
      character*512 secf
      integer L
      integer sec(nmaxres)
c     local variables
      integer i,j
      character*3 res
      open(unit=10,file=secf,status='old',err=666)
      L=1
 10   read(10,*,end=20)i,res,sec(L),j
      L=L+1
      goto 10
 20   L=L-1
      return
 666  call exit(1)
      end
c     ################################################


c     ################################################
      subroutine readsnapshot(trainp,xyz,header,snapshot,L)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
c     input variables
      integer trainp
      real xyz(3,nmaxres)       !coordinates array
      character*512 header
      character*512 snapshot(nmaxres)
      integer L
c     local variables
      character*512 line
      integer i,iL
      read(trainp,'a512',err=666,end=666)header(1:512)
      write(L,*)header !first number of header is sequence length
      do iL=1,L
         read(trainp,'a512')line(1:512)
         write(snapshot(iL),'a512')line
         read(snapshot(iL),*)xyz(1,iL),xyz(2,iL),xyz(3,iL) 
      enddo
      do i=1,3
         xyz(i,L)=xyz(i,L)/0.87 !from Angstroms to lattice units
      enddo
      return
 666  L=0 !signal error
      end !end of subroutine readsnapshot
c     ################################################


c     ###################################################
      subroutine prepare_vectors(xyz,bondvec,bonduni,
     $     bisector,hbvec,L)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
c     input variables
      integer L
      real xyz(3,nmaxres)
      real bondvec(3,nmaxres),bisector(3,nmaxres), hbvec(3,nmaxres)
      real bonduni(3,nmaxres)
c     local variables
      integer i,j,iL
      real dd

c     calculate bond vectors and unit bond vectors
      do iL=2,L
         dd=0.0
         do j=1,3
            bondvec(j,iL)=xyz(j,iL)-xyz(j,iL-1)
            dd=dd+bondvec(j,iL)*bondvec(j,iL)
         enddo
         do j=1,3
            bonduni(j,iL)=bondvec(j,iL)/sqrt(dd)
         enddo
      enddo

c     calculate bisector vectors
      do iL=2,L-1
         dd=0.0
         do j=1,3
            bisector(j,iL)=bondvec(j,iL+1)-bondvec(j,iL)
            dd=dd+bisector(j,iL)*bisector(j,iL)
         enddo
         do j=1,3
            bisector(j,iL)=bisector(j,iL)/sqrt(dd) !normalize
         enddo
      enddo

c     calculate hydrogen bond vectors
      do iL=2,L-1
         dd=0.0
         hbvec(1,iL)=bonduni(2,iL)*bonduni(3,iL+1)-
     $        bonduni(3,iL)*bonduni(2,iL+1)
         dd=dd+hbvec(1,iL)*hbvec(1,iL)
         hbvec(2,iL)=bonduni(3,iL)*bonduni(1,iL+1)-
     $        bonduni(1,iL)*bonduni(3,iL+1)
         dd=dd+hbvec(2,iL)*hbvec(2,iL)
         hbvec(3,iL)=bonduni(1,iL)*bonduni(2,iL+1)-
     $        bonduni(2,iL)*bonduni(1,iL+1)
         dd=dd+hbvec(3,iL)*hbvec(3,iL)
         do j=1,3
            hbvec(j,iL)=hbvec(j,iL)/sqrt(dd)
         enddo
      enddo

      end !end of subroutine prepare_vectors
c     ################################################


c     #################################################
      subroutine hydrogenbondmap(xyz,bondvec,bonduni,
     $     bisector,hbvec,L,sec,hbmap)
c     not optimized, coded for easyness of reading.
c     input parameters
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
      integer L
      real xyz(3,nmaxres)
      real bondvec(3,nmaxres),bisector(3,nmaxres), hbvec(3,nmaxres)
      real bonduni(3,nmaxres)
      integer sec(nmaxres)
      real hbmap(nmaxres,nmaxres)
c     local variables
      integer i,j
      real ehbij(nmaxres,nmaxres) !reinforcement for bond between similar secondary assignments
c     init hydrogen bond parameters
      real cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      real cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      character *512 la,lb
      integer ires,jres,idist
      real relvec(3) !relative vector
      real cr2,bb,cc,av11,av12,av21,av22,dotproduct
      real fact,energyHBa,energyHBb
      write(la,*)'48 0.45 0.1 0.00'
      read(la,*)cr2a,acut_bb,acut_cc,acut_vv
      write(lb,*)'50 0.25 0.4 0.35'
      read(lb,*)cr2b,bcut_bb,bcut_cc,bcut_vv
      call set_ehb(sec,ehbij,L)    
      do ires=2,L-2
         do jres=ires+1,L-1
            hbmap(ires,jres)=0.0
            fact=-1.0
            cr2=0.0
            do i=1,3 !initialize relative vector
               relvec(i)=xyz(i,jres)-xyz(i,ires)
               cr2=cr2+relvec(i)*relvec(i)
            enddo
            idist=jres-ires
            bb=dotproduct(hbvec,ires,jres)
            cc=dotproduct(bisector,ires,jres)
            av11=dotproduct(bondvec,ires,jres)
            av12=dotproduct(bondvec,ires,jres+1)
            av21=dotproduct(bondvec,ires+1,jres)
            av22=dotproduct(bondvec,ires+1,jres+1)
c     alpha-helix 
            if(idist.eq.3 .and.
     $           sec(ires).ne.4.and.sec(jres).ne.4 .and. !both not assigned as helixes
     $           cr2.lt.cr2a .and. !|r_i-r_k|^2 < Cr2a=36Angstroms^2
     $           cc.gt.acut_cc .and. !c_i*c_k < acut_cc=0.1, c_i:bisector vector
     $           bb.gt.acut_bb .and. !acut_bb=0.45
     $           av11.gt.acut_vv .and. !acut_vv=0.0
     $           av22.gt.acut_vv
     $           )then
               fact=(1-abs(cc-0.4))*(1-abs(bb-0.815)) !deviation from optimal
               hbmap(ires,jres)=energyHBa(ires,jres,fact,relvec,hbvec,
     $              ehbij)
c     antiparallel-sheet, short range
            elseif(idist.gt.4.and.idist.lt.20 .and.
     $              sec(ires).ne.2.and.sec(jres).ne.2 .and.
     $              cr2.lt.cr2b .and. !cr2b=38Angstroms^2
     $              cc.gt.bcut_cc .and. !bcut_cc=0.4
     $              bb.lt.-bcut_bb .and. !antiparallel, bcut_bb=0.25
     $              av12.lt.-bcut_vv .and. !bcut_vv=0.35
     $              av21.lt.-bcut_vv
     $              )then
               fact=abs(bb)*cc  !bb->-1,cc->1
               hbmap(ires,jres)=energyHBb(ires,jres,fact,relvec,hbvec,
     $              ehbij)
c     (anti)parallel sheet, long range
            elseif( (idist.ge.20 .and.
     $              sec(ires).ne.2.and.sec(jres).ne.2 .and.
     $              cr2.lt.cr2b.and.cc.gt.bcut_cc) .and.
     $              ((bb.lt.-bcut_bb .and. !antiparallel-sheet
     $              av12.lt.-bcut_vv .and.
     $              av21.lt.-bcut_vv) .or.
     $              (bb.gt.bcut_bb .and. !parallel-sheet
     $              av11.gt.bcut_vv .and.
     $              av22.gt.bcut_vv))
     $              )then
               fact=bb*cc       !bb->1,cc->1
               hbmap(ires,jres)=energyHBb(ires,jres,fact,relvec,hbvec,
     $              ehbij)
            endif
         enddo
      enddo
      end !end of subroutine hbmap
c     ################################################


c     ###########################################################
      function dotproduct(vectors,i,j)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
c     input variables
      real vectors(3,nmaxres)
      integer i,j,k
c     local variables
      real dotproduct
      dotproduct=0.0
      do k=1,3
         dotproduct=dotproduct+vectors(k,i)*vectors(k,j)
      enddo
      end
c     ################################################


c     ##########################################################
      subroutine set_ehb(sec,ehbij,L)
c     initialize ehbij matrix
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
c     input variables
      integer sec(nmaxres)
      real ehbij(nmaxres,nmaxres)
      integer L
c     local variables
      integer ires,jres,is,js
      do ires=1,L
         is=sec(ires)
         do jres=1,L
            js=sec(jres)
            ehbij(ires,jres)=1
            if(iabs(ires-jres).eq.3.and.is.eq.2.and.js.eq.2)then
               ehbij(ires,jres)=ehbij(ires,jres)+0.5 !hydrogen bond for helix enhanced
            endif
            if( (is.eq.4.or.js.eq.4).and. !at least one strand 
     $           is*js.ne.8.and.          !no helix
     $           iabs(ires-jres).gt.4
     $           ) then
               ehbij(ires,jres)=ehbij(ires,jres)+0.5
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
      function getnhb(hbmap,L)
      implicit none
      integer nmaxres
      real zero
      parameter(nmaxres=2000) !maximum number of residues
      parameter(zero=-0.00001)
c     input variables
      real hbmap(nmaxres,nmaxres)
      integer L
c     local variables
      integer getnhb
      integer i,j
      getnhb=0
      do i=1,L-1
         do j=i+1,L
            if(hbmap(i,j).lt.zero) getnhb=getnhb+1
         enddo
      enddo
      end !end of function getnhb
c     ################################################


c     ################################################
      subroutine savesnapshot(traoutp,header,snapshot,L)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000) !maximum number of residues
c     input variables
      integer traoutp
      character *512 header
      character *512 snapshot(nmaxres)
      integer L
c     local variables
      integer lnblnk
      integer iL
      character*512 line
      write(traoutp,'a')header(1:lnblnk(header))
      do iL=1,L
         write(line,'a') snapshot(iL)
         write(traoutp,'a')line(1:lnblnk(line))
      enddo
      end !end of subroutine savesnapshot
c     ################################################
