*************************************************************************
*     pgf77 -g -Mbounds -Mchkfpstk -Mextend  -Wl,-static -o hbnet.x hbnet.f
*     To debug, type in fec02: pgdbg hbnet.x
*     pgf77 -Mextend -O -s -fast -Wl,-static -o hbnet.x hbnet.f
*
*     This program outputs weighted contact map where each contact represents a
*     hydrogen bond.

      program hbnet
      implicit none
      integer nmaxres
      parameter(nmaxres=2000) !maximum number of residues
      logical readcommand
      character*512 pdbf !input pdb file
      character*512 outf !output contact map filename
      character*512 secf !seq.dat file
      real xyz(3,nmaxres)     !coordinates array
      integer L               !sequence length
      integer resseq(nmaxres) !residue index vector
      integer sec(nmaxres)
      real bondvec(3,nmaxres),bonduni(3,nmaxres)
      real bisector(3,nmaxres), hbvec(3,nmaxres)
      real hbmap(nmaxres,nmaxres)
c     manage input command arguments
      if(readcommand(pdbf,secf,outf).eq..true.) call exit(1)
c     read input PDB structure and scale coordinates by 0.87
      call readPDBfile(pdbf,xyz,L,resseq)
c     read input secondary structure assignments
      call readSeqDatFile(secf,L,sec)
c     prepare bond, bisector, and hydrogen-bond vectors
      call prepare_vectors(xyz,bondvec,bonduni,bisector,hbvec,L)
c     calculate hydrogen bond contact map
      call hydrogenbondmap(xyz,bondvec,bonduni,bisector,hbvec,L,
     $    sec,hbmap)
      call outputMap(outf,hbmap,L)
      end !end of program hbnet
c     ################################################


c     ################################################
      function readcommand(pdbf,secf,outf)
c     read command line arguments
      logical readcommand
      character*512 fnam
      character*512 secf
      character*512 pdbf        !input pdb file
      character*512 outf !output contact map filename
      integer ireq,nreq !nreq: number of Required arguments
      integer narg      !number of words in the command line
      integer i
      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h'
     $     .or.fnam.eq.'-help')then
 100     call system('clear')
         write(*,*)'Usage: hbnet [options]:'
         write(*,*)' Required arguments:'
         write(*,*)' -pdbf pdb file'
         write(*,*)' -secf seq.dat file'
         write(*,*)' -outf out file name storing contact map'
         readcommand=.true. !signal error
         return
      endif
      i=0
      nreq=3
      ireq=0
      narg=iargc()
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-pdbf')then
         i=i+1
         call getarg(i,fnam)
         write(pdbf,*) fnam(1:lnblnk(fnam))
         ireq=ireq+1
      elseif(fnam.eq.'-secf')then
         i=i+1
         call getarg(i,fnam)
         write(secf,*) fnam(1:lnblnk(fnam))
         ireq=ireq+1
      elseif(fnam.eq.'-outf')then
         i=i+1
         call getarg(i,fnam)
         write(outf,*) fnam(1:lnblnk(fnam))
         ireq=ireq+1
      endif
      if(i.lt.narg) goto 115
      if(ireq.lt.nreq) goto 100 !print welcome message and output program

      end !end of subroutine readcommand
c     ################################################


c     ###################################################
      subroutine readPDBfile(pdbf,xyz,L,resseq)
c     initialize xyz,L,resseq,xxx
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
      character*512 pdbf        !pdb file name
      character*100 s,junk,junk2
      real xyz(3,nmaxres)       !coordinates array
      integer i,L               !L: sequence length
      integer resseq(nmaxres)   !residue index vector

 104  format(A100)
 9000 format(A22,i4,A4,3F8.3)
      L=0
      open(unit=10,file=pdbf,status='old',err=666)
 101  read(10,104,end=102) s
         if(s(1:3).eq.'TER' .or. s(1:3).eq.'END') goto 102
         if(s(1:4).eq.'ATOM')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &           eq.'  CA')then
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  L=L+1
                  if(L.ge.nmaxres)goto 102 !too long sequence. Consider only first nmax residues
                  read(s,9000)junk,resseq(L),junk2,
     $                 xyz(1,L),xyz(2,L),xyz(3,L) 
                  do i=1,3
                     xyz(i,L)=xyz(i,L)/0.87 !from Angstroms to lattice units
                  enddo
               endif
            endif
         endif
      goto 101
 102  close(10)
      return
 666  call exit(1)
      end
c     ################################################


c     ###################################################
      subroutine readSeqDatFile(secf,L,sec)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)   !maximum number of residues
c     input variables
      character*512 secf
      integer L,sec(nmaxres)
c     local variables
      integer i,j,iL
      character*3 res
      open(unit=10,file=secf,status='old',err=666)
      do iL=1,L
         read(10,*)i,res,sec(iL),j
      enddo
      return
 666  call exit(1)
      end
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

      do ires=1,L
         do jres=1,L
            hbmap(ires,jres)=0.0
         enddo
      enddo
      write(la,*)'48 0.45 0.1 0.00'
      read(la,*)cr2a,acut_bb,acut_cc,acut_vv
      write(lb,*)'50 0.25 0.4 0.35'
      read(lb,*)cr2b,bcut_bb,bcut_cc,bcut_vv
      call set_ehb(sec,ehbij,L)    
      do ires=1,L-1
         do jres=1,L-1
            if(ires.eq.jres) goto 10
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
c     alpha-helix bond is probed by the (i,i+3) pair instead of the (i,i+4) pair.
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
c               write(*,*)ires,jres,hbmap(ires,jres)
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
c               write(*,*)ires,jres,hbmap(ires,jres)
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
               fact=abs(bb)*cc       !bb->1,cc->1
               hbmap(ires,jres)=energyHBb(ires,jres,fact,relvec,hbvec,
     $              ehbij)
c               write(*,*)ires,jres,hbmap(ires,jres)
            endif
 10      enddo
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


c     #################################################
      subroutine outputMap(outf,hbmap,L)
      implicit none
      integer nmaxres
      parameter(nmaxres=2000)
c     input/output variables
      character*512 outf        !output contact map filename
      character*100000 buf
      character*512 buf2
      integer L
      real hbmap(nmaxres,nmaxres)
c     local variables
      integer ires,jres,lnblnk,numhb  
      real tote,eperbond !total-energy,energy-per-bond
      numhb=0
      tote=0.0
      write(buf,'a')''
      do ires=1,L
         do jres=1,L
            if(ires.ne.jres.and.hbmap(ires,jres).lt.0.0)then
               numhb=numhb+1
               tote=tote+hbmap(ires,jres)
               write(buf2,'i4xxi4xxf6.3')ires,jres,hbmap(ires,jres)
               write(buf,'a')buf(1:lnblnk(buf))//buf2(1:lnblnk(buf2))//
     $              '\n'
            endif
         enddo
      enddo
      eperbond=tote/numhb
      write(buf2,'af7.2ai4af6.3')'#E=',tote,' N=',numhb,' e=',eperbond
      write(buf,'a')buf2(1:lnblnk(buf2))//'\n'//buf(1:lnblnk(buf))
      open(unit=1,file=outf,status='unknown')
      if(lnblnk(buf).gt.1) write(1,'a') buf(1:lnblnk(buf)-1)
      close(1)
      end !end of subroutine outputMap
c     ################################################
