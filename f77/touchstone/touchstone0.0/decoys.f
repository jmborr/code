
*****************************************************************************
*     This program is to generate decoy structures by PHS Monte Carlo       *
*     simulations under CABS lattice model.                                 *
*     Comments and bug report should be addressed to: zhang6@buffalo.edu    *
*     if there is gap, proper number of GLY will be added into the gap,     *
*     RMSD and output will only based on the residues in CA                 *
*****************************************************************************
*     do 5555 i_start=1,niter
*      do 1111 icycle=1,ncycle
*       do 2222 itemp=1,N_rep
*        do 3333 iphot=1,phot
*         do 4444 i_lch=1,Lch
*           +++
*         enddo 4444
*         record decoys
*        enddo 3333
*        if(done) goto 6666
*       enddo 2222
*      enddo 1111
*      adjust er1,er2,er3,er4
*     enddo 5555
*     6666

c        1         2         3         4         5         6         7 !
c 345678901234567890123456789012345678901234567890123456789012345678901234567890
      program Cabs
      implicit integer(i-z)
      parameter(ndim=1000)       !maximum length of chain-length
      parameter(nrep=100)      !maximum number of replicas
      parameter(nvec=312)       !number of vectors
      parameter(ndcy=100)       !maximum number of decoys in each bin
      parameter(nbin=20)        !maximum number of bins
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/three/angle(nvec,nvec)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
      common/hopp/eonehw(0:19)
      common/looks/exc,exc1,exc2
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc

      common/CA/nCA,axca(ndim),ayca(ndim),azca(ndim)
      common/initial/x0(ndim),y0(ndim),z0(ndim)
      common/addition/i_gap,s_gap(100),n_add(100),iLch(ndim),inCA(ndim)
      common/dcyx/dex(nbin,ndcy,ndim)
      common/dcyy/dey(nbin,ndcy,ndim)
      common/dcyz/dez(nbin,ndcy,ndim)
      common/dcyxg/gex(nbin,ndcy,ndim)
      common/dcyyg/gey(nbin,ndcy,ndim)
      common/dcyzg/gez(nbin,ndcy,ndim)
      dimension armsd_r(nbin,ndcy),drmsd_r(nbin,ndcy)
      dimension dexa(nbin,ndcy,ndim)
      dimension deya(nbin,ndcy,ndim)
      dimension deza(nbin,ndcy,ndim)

      common/arandom/  aarand,abrand,acrand,adrand
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12) !safe when vr^2<30
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/fr/frga(ndim),frgb(ndim)
      COMMON/RES/ ER3, Mcom(ndim),Kcom(ndim,50)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/maxi/maxin,vect1,vect2
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/forpreparemove4/ asrr(0:19,0:19,-12:12)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/distres/er4,er2,es3c
      common/rmsdrange/nca1,nca2
      common/msichores/msicho
      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB2,EHB3,EHB4,EHB5,EHB6
      common/ehbenergy1/EHB5a,EHB5b
      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT3a,ESHORT4
      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT9a
      common/eshortenergy4/ESHORT9b,ESHORT9c
      common/otherenergy/E_cord,E_cnum
      common/resnumber/Ncom,Ndis,accur
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/outputxyz/fxyz(3,ndim)
      character*3 sequ
      character*100 sss
      common/aminoacid/sequ(ndim)
      dimension Mphot(1000)

      common/temperature/itemp,atemp
      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
      common/pair1/eh2,eh1b
      dimension E_s(nrep),E_ss(nrep)
      common/paircut/ash

      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/pairmap1/apba(ndim,ndim),apbp(ndim,ndim),apbm(ndim,ndim)
      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap

      common/moveretio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc
      common/icgg/ icg(ndim), EH6  
      common/rs/i_thr0
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/nrepfile/n_repf

      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/mng/m_g(100)
      common/acct/accept0
      character*6 mname
      character fn
      common/movename/mname(100)
      common/readinitial/m_initial
      common/eall/N_sum(100),energ_sum(100),energ_sum2(100),E_min

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/chainm/mv(ndim)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) !SG
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) !cc
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) !hb

      dimension arm_low(nbin),arm_high(nbin),arm_cut(nbin),n_dec(nbin)


ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /1000*1.0/
ccc   

ccccccccccccccccccccccccc common input files cccccccccccccccccccccccccccccc
      open(unit=1,file='centro.comm',  status='old') !eoinp(A,dis)
      open(unit=2,file='profile3_new_new.comm',status='old') !envir
      open(unit=3,file='quasi3_new.comm',  status='old') !pairwise
      open(unit=4,file='sidecent_yang.comm',status='old') !for Sc position.
      open(unit=5,file='r13.comm',    status='old') !E_short of (i,i+2)
      open(unit=12,file='r14.comm',   status='old') !E_short of (i,i+3)
      open(unit=7,file='r14h.comm',   status='old') !stress helical-stru.
      open(unit=8,file='r14e.comm',   status='old') !stress extended-stru.
      open(unit=9,file='r15.comm',   status='old') !E_short(i,i+4)
      open(unit=10,file='r15h.comm',   status='old')
      open(unit=11,file='r15e.comm',   status='old')
      open(unit=22,file='s1234h.comm', status='old') !E_short for Sg
      open(unit=23,file='s1234e.comm', status='old')
      open(unit=25,file='concut.comm', status='old') !cutoff of contact predi.

ccccccccccccccc sequence specified input files cccccccccccccccc///////////
      open(unit=13,file='CA',status='old') !call read_inital
      open(unit=14,file='seq.dat',     status='old')
      open(unit=16,file='comb_native.dat',status='old') !from threading
      open(unit=17,file='dist_native.dat',status='old') !threading

      open(unit=19, file='in.dd',     status='old')
      open(unit=20, file='out.d',     status='unknown')
cccccccccccccccc for E-t cccccccccccccccccccccccccccccccccccccccccccccccccc
c     open(unit=91, file='swepa.d',    status='unknown') !$$
c     open(unit=92, file='swepb.d',    status='unknown') !$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      read(19,*) random,niter,ncycle,phot,N_rep
      read(19,*) h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
      read(19,*) atemp2,atemp1,exc,Mend
      read(19,*) N_decoys,N_bins

      read(19,*) eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh !alpha-type HB
      read(19,*) eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh !alpha-type HB

      read(19,*) eh1,eh2,eh3,eh4  !for EHBs (6)
      read(19,*) eh1a,eh1b
      read(19,*) es1,es2,es3,es4,es5,es6     !for ESHORTs (6)
      read(19,*) es3a,es3b,es3c,es7,es7a,es7b,es7c
      read(19,*) en1,en2,en3                 !for ensemble (3)
      read(19,*) er1,er2,er3,er4             !for restrains (4)

ccccccccccccccctrajectory files cccccccccccccccccccccccccccccccccccccc
      do i=1,N_bins
         arm_low(i)=i-1
         arm_high(i)=i
         arm_cut(i)=0.5*((arm_low(i)+arm_high(i))/2) !dist between structures
         n_dec(i)=0
         if(i.lt.10)then
            fn=char(48+i)
            open(unit=30+i,file=fn//'.PDB',status='unknown')
         else
            fn=char(48+(i-10))
            open(unit=30+i,file='1'//fn//'.PDB',status='unknown')
         endif
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call prepare_vectors      !prepare all possible bond-vectors
      call prepare_neighbors    !define goodc(i,j), angle(i,j), prod(i,j)
      call read_initial         !read initial (x,y,z) from 'CA'

      call set_common           !set common parameters
      call read_seq             !read seq(i),sec(i) from 'seq.dat'
      call read_centro          !read eoinp(i,dis) from 'centro.comm'
      call read_profile         !read envir(ia,im,ip,i,j) from 'profile3.comm'
      call read_E13             !read 1-3 short-range E from 'r13.comm'
      call read_E14             !read 1-4 potential from 'r14*.dat'
      call read_E15             !read 1-5 potential from 'r15*.dat'
      call read_pair            !read 'quarsi3.comm' and 'pair3.dat'
      call read_concut          !read cut-off for contact prediction
      call read_contactrestrain !read contact restrains from 'comb.dat'
      call read_distantrestrain !read distant restrains from 'dist.dat'
      call set_temperature      !set temperature for different replic
      call set_EHB              !set structure-specitic H-bond, EHBIJ(i,j)

      call prepare_beta         !define C_beta, C_group, and hydrogen-bond
      call prepare_frg          !compute the secondary fragment biases

      call get_acorder          !calculate contact order
      call write_parameter      !print out initial parameters

      call set_move_retio       !set movement percentage
      call prepare_move2        !2-bond move, num2=26784

      do i=1,100
         bNSa(i)=0              !aceptance for swep
         bNSt(i)=0

         bNa(i)=0               !acceptance for move2,3,4,5,6,7
         bNt(i)=0

         bNNa(i)=0              !acceptance for different temperature.
         bNNt(i)=0

         N_sum(i)=0
         energ_sum(i)=0         !<E_tot>
         energ_sum2(i)=0        !<E_tot^2>
      enddo
      E_min=10000

ccc   reset phot ------------------------------------>
      n=N_rep/3.5
      do i=1,N_rep
         if(i.lt.n)then
            Mphot(i)=phot/2+phot/3*abs(i-n)**0.6
         else
            Mphot(i)=phot/2+phot/3*abs(i-n)**0.6
         endif
         if(i.gt.N_rep/4.5.and.i.lt.N_rep*2/3)then
            Mphot(i)=phot/2
         endif
         write(20,*)i,Mphot(i),aT_rep(i)
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc       The main cycle start from here !                         ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 5555 i_start=1,niter
ccc   starting from CA in each cycle:
         do i=1,N_rep
            do j=1,Lch
               xrep(j,i)=x0(j)
               yrep(j,i)=y0(j)
               zrep(j,i)=z0(j)
            enddo
         enddo
         do 1111 icycle=1,ncycle
            do 2222 itemp=1,N_rep !iterate for all the replicas
               atemp=aT_rep(itemp) !current temperature
               call set_current !get current (x,y,z,ica)
               call initial_move !update center, axis, energy
ccc
               do 3333 iphot=1,Mphot(itemp) !iterate at fixed temperature
                  do 4444 i_lch=1,Lch
                     fff=aranzy(no)
                     if(fff.le.bh2)then
                        call move2
                     elseif(fff.le.bh3s)then
                        call move3s
                     elseif(fff.le.bh3d)then
                        call move3d
                     elseif(fff.le.bh4s)then
                        call move4s
                     elseif(fff.le.bh4d)then
                        call move4d
                     elseif(fff.le.bh5s)then
                        call move5s
                     elseif(fff.le.bh5d)then
                        call move5d
                     elseif(fff.le.bh6)then
                        call move6
                     elseif(fff.le.bhendn)then
                        call move_n_end
                     else
                        call move_c_end
                     endif
 4444             continue
**************************************************************************
ccc   check and write decoys -------------------------->
                  do i=1,nCA
                     r_1(1,i)=axca(i)
                     r_1(2,i)=ayca(i)
                     r_1(3,i)=azca(i)
                     j=iLch(i)
                     r_2(1,i)=x(j)*0.87
                     r_2(2,i)=y(j)*0.87
                     r_2(3,i)=z(j)*0.87
                  enddo
                  call u3b(w,r_1,r_2,nCA,0,rms,u,t,ier)
                  aaa=dsqrt(rms/float(nCA)) !RMSD12
                  do i=1,N_bins
                     if(aaa.ge.arm_low(i).and.aaa.lt.arm_high(i))then
                        if(n_dec(i).ge.N_decoys)goto 41
                        do k=1,n_dec(i) !check distance with others
                           do m=1,nCA
                              j=iLch(m)
                              r_1(1,m)=x(j)*0.87
                              r_1(2,m)=y(j)*0.87
                              r_1(3,m)=z(j)*0.87
                              r_2(1,m)=dex(i,k,m)
                              r_2(2,m)=dey(i,k,m)
                              r_2(3,m)=dez(i,k,m)
                           enddo
                           call u3b(w,r_1,r_2,nCA,0,rms,u,t,ier)
                           armsd=dsqrt(rms/float(nCA)) !RMSD12
                           if(armsd.lt.arm_cut(i)) goto 41
                        enddo
ccc   record decoys:
                        n_dec(i)=n_dec(i)+1 !number of decoys in each bins
                        do m=1,nCA
                           j=iLch(m)
                           dex(i,n_dec(i),m)=x(j)*0.87
                           dey(i,n_dec(i),m)=y(j)*0.87
                           dez(i,n_dec(i),m)=z(j)*0.87
                           j1=ica(j-1)
                           j2=ica(j)
                           jm=seq(j)
                           gex(i,n_dec(i),m)=(x(j)+GX(j1,j2,jm))*0.87
                           gey(i,n_dec(i),m)=(y(j)+GY(j1,j2,jm))*0.87
                           gez(i,n_dec(i),m)=(z(j)+GZ(j1,j2,jm))*0.87
                        enddo
                        call cal_drmsd(i,n_dec(i),drmsd)
                        armsd_r(i,n_dec(i))=aaa
                        drmsd_r(i,n_dec(i))=drmsd
                        do m=1,Lch
                           dexa(i,n_dec(i),m)=x(m)*0.87
                           deya(i,n_dec(i),m)=y(m)*0.87
                           deza(i,n_dec(i),m)=z(m)*0.87
                        enddo
ccc   
                     endif
                  enddo
 41               continue
c                  do i=1,N_bins
c                     write(*,*)i_start,icycle,itemp,i,n_dec(i),Lch
c                  enddo
**************************************************************************
 3333          continue
               E_rep(itemp)=energy_tot() !whole energy
               if(E_rep(itemp).lt.E_min) E_min=E_rep(itemp)
               do i=1,Lch
                  xrep(i,itemp)=x(i)
                  yrep(i,itemp)=y(i)
                  zrep(i,itemp)=z(i)
               enddo

ccc   check number of decoys:
               n=0
               do i=1,N_bins
                  if(n_dec(i).ge.N_decoys)then
                     n=n+1
                  endif
               enddo
               if(n.ge.N_bins)goto 6666
ccc   
 2222       continue

ccccccccccccccccc<RMSD>, <E> cccccccccccccccccccccccccc
            do i=1,N_rep
               energ_sum(i)=energ_sum(i)+E_rep(i)
               energ_sum2(i)=energ_sum2(i)+E_rep(i)*E_rep(i)
               N_sum(i)=N_sum(i)+1
            enddo

ccccccccccccccccc swap replicas cccccccccccccccccccccccccccccccc
            if(icycle.eq.icycle/2*2)then
               do i=1,N_rep-1,2 !swap odd replicas
                  call swap(i,i+1)
               enddo
            else
               do i=2,N_rep-1,2
                  call swap(i,i+1) !swap even replicas
               enddo
            endif
 1111    continue
         do i=1,N_bins
            write(*,*)i_start,i,n_dec(i)
         enddo
ccc   adjust restraints ----------->
         if(n_dec(1).ge.N_decoys.and.n_dec(N_bins).lt.N_decoys)then
            er1=er1/2           !too many low RMSD decoys
            er2=er2/2
            er3=er3/2
            er4=er4/2
            icycle=icycle+10
         endif
         if(n_dec(1).lt.N_decoys.and.n_dec(N_bins).gt.N_decoys)then
            er1=er1*2           !too many high RMSD decoys
            er2=er2*2
            er3=er3*2
            er4=er4*2
            icycle=icycle/2
            if(icycle.le.3)icycle=3
         endif
 5555 continue
c--------------------------Main cycle ended here!!---------------
 6666 continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ----------- print out all the decoys ---------->
      write(20,*)'    i      N_decoys(i)'
      do k=1,N_bins
         write(20,*)k,n_dec(k)
         do j=1,N_dec(k)
            nf=30+k
            write(nf,1234)j,armsd_r(k,j),drmsd_r(k,j)
 1234       format('HEADER   i=',I5,'  RMSD=',f8.3,'  dRMSD=',f8.3)
            do i=1,nCA
               m=iLch(i)
               write(nf,1237)i,sequ(m),i,
     $              dex(k,j,i),dey(k,j,i),dez(k,j,i)
c               write(nf,1238)i,sequ(m),i,
c     $              gex(k,j,i),gey(k,j,i),gez(k,j,i)
            enddo
         enddo
      enddo
 1237 format('ATOM  ',I5,'  CA',A5,I6,4X,3f8.3)
 1238 format('ATOM  ',I5,'  SC',A5,I6,4X,3f8.3)
c^^^^^^^^^^^^^^^^^ decoys finished ^^^^^^^^^^^^^^^^^^^^^^^^^

      i7=77                     
      if(i7.eq.88) goto 303     !without output full atomic model
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,N_bins
         if(i.lt.10)then
            fn=char(48+i)
            open(unit=50+i,file=fn//'.PDB_all',status='unknown')
         else
            fn=char(48+(i-10))
            open(unit=50+i,file='1'//fn//'.PDB_all',status='unknown')
         endif
      enddo
c     ----------------- run jackel --------------------->
      ienvironment=putenv('JACKALDIR=/nfs/users/yzhang/bin/'//
     $     'jackal/bin/jackal.dir')
      ienvironment=putenv('PATH=$PATH:/nfs/users/yzhang/bin/jackal/bin')
c     ----------- print out all the decoys ---------->
      do k=1,N_bins
         do j=1,N_dec(k)
            nf=50+k
            write(nf,1234)j,armsd_r(k,j),drmsd_r(k,j)
c     full-chain model:
            open(unit=99,file='tmp.CA',status='unknown')
            do i=1,Lch
               write(99,1237)i,sequ(i),i,dexa(k,j,i),deya(k,j,i),deza(k,j,i)
            enddo
            close(99)
c     run jackel:
            ioutput=system("/bin/sleep 1")
            ioutput=system("ctrip -prm 2 -k 0 tmp.CA 2> /dev/null")
            ioutput=system("/bin/sleep 1")
            mm=0
            open(unit=98,file='tmp.CA_fix.pdb',status='old')
            do while (.true.)
               read(98,'A100',end=1000)sss
               read(sss,9000)du,idu,atom,seqT,ii,bbx,bby,bbz
 9000          format(A6,i5,A5,1x,A3,2x,i4,4x,3F8.3)
               do m=1,nCA
                  if(ii.eq.iLch(m))then
                     mm=mm+1
                     write(nf,9000)'ATOM  ',mm,atom,seqT,m,bbx,bby,bbz
                     goto 304
                  endif
               enddo
 304           continue
            enddo
 1000       continue
            close(98)
            write(nf,'A3')'TER'
         enddo
      enddo
c^^^^^^^^^^^^^^^^^ full atom decoys finished ^^^^^^^^^^^^^^^^^^^^^^^^^
 303  continue

      write(20,*)
      write(20,*)'final er1234=',er1,er2,er3,er4
      write(20,*)'escape i_start=',i_start
      write(20,*)

      call test_neighbor        !check all the neighboring residues
      call test_overlap         !test the overlap of C_a and C_b

      write(20,*)'E_final=',energy_tot()

cccccccccccccccccccccccc Na/Nt cccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*) 'i_move    move   Na(i)  Nt(i)   Na(i)/Nt(i)'
      do i=2,15
         if(bNt(i).gt.1)then
            write(20,5004) i,mname(i),bNa(i),bNt(i),bNa(i)/bNt(i)
         else
            write(20,5004) i,mname(i),bNa(i),bNt(i)
         endif
      enddo
 5004 format(I4,A9,2f15.1,f11.6)
      
ccccccccccccccccccccccccccE_final, NSa/NSt ccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'-------------- E_final, Na_swap/Nt_swap ---------'
      WRITE(20,*) 'i  T(i) final_E(i)  Nsa(i)  Nst(i)  Nsa(i)/Nst(i)'
      do i=1, n_rep
         if(bNSt(i).gt.1)then
            write(20,5005) i,aT_rep(i),E_rep(i),
     $           bNSa(i),bNSt(i),bNSa(i)/bNSt(i)
         else
            write(20,5005) i,aT_rep(i),E_rep(i),bNSa(i),bNSt(i)
         endif
      enddo
 5005 format(I4,f7.2,f8.1,2f15.1,f11.6)

ccccccccccccccccccccccc <E>, NNa/NNt ccccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'------------ <energy>, Na(i)/Nt(i) ----------------'
      write(20,*)'i_rep  T   <E>   NNa(i_temp)  NNt(i_temp)  Na/Nt'
      do i=1,N_rep
         energ_sum(i)=energ_sum(i)/(N_sum(i)+0.0001)
         energ_sum2(i)=energ_sum2(i)/(N_sum(i)+0.0001)
         if(bNNt(i).gt.1)then
            cheat=(energ_sum2(i)-energ_sum(i)**2)/(aT_rep(i)**2)
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i),bNNa(i)/bNNt(i)
         else
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i)
         endif
      enddo
 5006 format(I4,f7.2,f8.1,f15.3,2f12.1,f11.6)
      write(20,*)'E_min=',E_min
      
      STOP
      END
ccccccc=======================================================cccccccccc
cc          The main program ended!
ccccccc=======================================================cccccccccc




































cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prepare possible vector v=(vx,vy,vz) satisfy |v*v|=17,18,19,20,21,22
c    and vector(-5:5,-5:5,-5:5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_vectors
      implicit integer(i-z)
                parameter(nvec=312)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension n(100)

      do i=10,30
         n(i)=0
      enddo

      nwmax=0
      aaa=0
      nn=5
      do x=-nn,nn
         do y=-nn,nn
            do z=-nn,nn
               vector(x,y,z)=0
               r=x*x+y*y+z*z
               if(r.ge.14.and.r.le.25) then
                  nwmax=nwmax+1
                  vx(nwmax)=x
                  vy(nwmax)=y
                  vz(nwmax)=z
                  vector(x,y,z)=nwmax
c                  write(*,*)nwmax,vx(nwmax),vy(nwmax),vz(nwmax),r,
c     $                 sqrt(float(r)*0.87*0.87)
                  n(r)=n(r)+1
                  aaa=aaa+sqrt(float(r)*0.87*0.87)
               endif
            enddo
         enddo
      enddo
      aaa=aaa/float(nwmax)
      write(*,*)'n1_all=',nwmax,'  <vr>=',aaa

c      do i=10,30
c         write(*,*)i,n(i),sqrt(float(i)*0.87*0.87)
c      enddo

c     i=1,5
c           10          24   2.751182    
c           11          24   2.885463    
c           12           8   3.013768    
c           13          24   3.136830    
c           14          48   3.255242    x
c           15           0   3.369496    
c           16           6   3.480000    x
c           17          48   3.587102    x
c           18          36   3.691097    x
c           19          24   3.792242    x
c           20          24   3.890758    x
c           21          48   3.986841    x
c           22          24   4.080662    x
c           23           0   4.172373    
c           24          24   4.262112    x
c           25          30   4.350000    x
c           26          72   4.436147    
c           27          32   4.520653    
c           28           0   4.603607    
c           29          72   4.685093    
c           30          48   4.765186    
c nwmax=         616  <vr>=   3.982909    
c nwmax=         312  <vr>=   3.809868    

c      stop
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c define goodc(i,j), angle(i,j), prod(i,j)
c     prod(i,j): v(i)*v(j).
c     angle(i,j): angle of neighbor bonds, i,j--->(1:nvec)
c     goodc(i,j): ture, when angle(i,j) in [60,160]; false, otherwise.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_neighbors
      implicit integer(i-z)
                parameter(nvec=312)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/three/angle(nvec,nvec)
      common/lengths/Lch,Lch1,Lch2

      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      max_m12=0
      do i=1,nvec
         m12(i)=0
      enddo

      mmm=0
      nnn=0
      kkk=0
      do i=1,nvec
      do j=1,nvec
         u21(i,j)=0
         a2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
         b2=vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)
         c2=(vx(i)+vx(j))**2+(vy(i)+vy(j))**2+(vz(i)+vz(j))**2
         cosangle=(a2+b2-c2)/(2*sqrt(a2*b2))
         angle(i,j)=acos(cosangle)*180/3.1415926
c     in database, angle is in [65,165];
         if(angle(i,j).gt.65.and.angle(i,j).lt.165)then
            goodc(i,j)=.true.
            mmm=mmm+1
            ijx=vx(i)+vx(j)
            ijy=vy(i)+vy(j)
            ijz=vz(i)+vz(j)
            do k=1,nvec
               if(vx(k).eq.ijx.and.vy(k).eq.ijy.and.vz(k).eq.ijz)then
                  kkk=kkk+1
                  u21(i,j)=k    !vi+vj=vk
                  m12(k)=m12(k)+1
                  u1(k,m12(k))=i
                  u2(k,m12(k))=j
                  if(max_m12.lt.m12(k))max_m12=m12(k)
c     write(*,*)i,j,k,angle(i,j),m12(k)
                  goto 10
               endif
            enddo
 10         continue
        else
            goodc(i,j)=.false.
         endif
         nnn=nnn+1
c         write(*,*)i,j,mmm,nnn,angle(i,j),goodc(i,j)
      enddo
      enddo

      n=0
      do i=1,nvec
         r=vx(i)**2+vy(i)**2+vz(i)**2
         if(m12(i).gt.0)n=n+1
c     if(r.gt.17)write(*,*)i,r,m12(i)
      enddo
c      write(*,*)'n2_all=',nnn
c      write(*,*)'n2good=',mmm,'  n21=',kkk
c      write(*,*)'n1_all=',nvec,'  n12=',n

c      stop
      return
      end

ccc   read (x,y,z) from 'CA' for ab initio cccccccc
ccc   we did not check the excluded volumn of initial chains!!!
      subroutine read_initial
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      character c1*4,c2*2,c3*3
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/readinitial/m_initial
      common/looks/exc,exc1,exc2

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/echain2/egx(ndim),egy(ndim),egz(ndim)

      character aline*100,du
      dimension ax(1000),ay(1000),az(1000)

      common/CA/nCA,axca(ndim),ayca(ndim),azca(ndim)
      common/initial/x0(ndim),y0(ndim),z0(ndim)
      common/addition/i_gap,s_gap(100),n_add(100),iLch(ndim),inCA(ndim)

      character*3 sequ
      common/aminoacid/sequ(ndim)

ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /1000*1.0/
ccc   

ccc   read data from CA file:
      i=0
      do while (.true.)
         read(13,'A100',end=1010) aline
         i=i+1
         read(aline,'A30,3F8.3')du,axca(i),ayca(i),azca(i)
      enddo
 1010 continue
 9000 format(A30,3F8.3)
      nCA=i
      write(*,*)'nCA=',nCA

ccc   check possible gaps in CA, if gaps, insert GLY
      j=0
      i_gap=0
      do i=1,nCA-1
         j=j+1
         iLch(i)=j
         inCA(j)=i
         ax(j)=axca(i)
         ay(j)=ayca(i)
         az(j)=azca(i)
         dis=sqrt((axca(i)-axca(i+1))**2+(ayca(i)-ayca(i+1))**2+
     $        (azca(i)-azca(i+1))**2)
         if(dis.gt.4.5)then
            i_gap=i_gap+1
            s_gap(i_gap)=i
            n_add(i_gap)=nint(dis/3.3)-1 !number of inserted residues
            do k=1,n_add(i_gap)
               j=j+1
               ax(j)=ax(j-1)+(axca(i+1)-axca(i))/float(n_add(i_gap)+1)
               ay(j)=ay(j-1)+(ayca(i+1)-ayca(i))/float(n_add(i_gap)+1)
               az(j)=az(j-1)+(azca(i+1)-azca(i))/float(n_add(i_gap)+1)
            enddo
         endif
      enddo
      i=nCA
      j=j+1
      iLch(i)=j
      inCA(j)=i
      ax(j)=axca(i)
      ay(j)=ayca(i)
      az(j)=azca(i)
      Lch=j
      write(*,*)i,iLch(i),Lch
c^^^^^^^^^^i_gap, s_gap(i), n_add(i) obtained ^^^^^^^^^^^^^^^^^^^^

ccc   set up initial lattice structure:
      do i=1,Lch
         ax(i)=ax(i)/0.87       !C_alpha scaled by 0.87
         ay(i)=ay(i)/0.87
         az(i)=az(i)/0.87
c         write(*,*)i,ax(i),ay(i),az(i)
      enddo

c     convert (ax,ay,az) into (x,y,z) without excluded volumn----------->
      x0(1)=nint(ax(1))
      y0(1)=nint(ay(1))
      z0(1)=nint(az(1))
      px=x0(1)
      py=y0(1)
      pz=z0(1)
      do 101 i=2,Lch
         armin=1000000.         !minimum distance between c_alpha and lattice
         jj=0                   !chosen vector
         do 100 j=1,nvec
            if(i.ge.3)then      !check good neighbor
               if(.not.goodc(jjjj,j)) goto 100
            endif
            jx=px+vx(j)
            jy=py+vy(j)
            jz=pz+vz(j)
            bx=float(jx)-ax(i)
            by=float(jy)-ay(i)
            bz=float(jz)-az(i)
            ar=bx*bx+by*by+bz*bz
            if(ar.lt.armin)then
               jj=j
               mx=jx
               my=jy
               mz=jz
               armin=ar
            endif
 100     continue
         if(jj.lt.1)then
            write(*,*)'UNSOLVABLE STERIC PROBLEM, exc=',exc,i
            stop
         endif
         jjjj=jj
         x0(i)=mx
         y0(i)=my
         z0(i)=mz
         px=mx
         py=my
         pz=mz
 101  continue

******prepare movement for normal movement *****************
      nfl=Lch
      do i=1,nfl
         ras(i)=i
      enddo
      call move_point           !decide movement point for notmal run
      nfr=0                     !number of frozen fragments.

c      do i=1,nfl
c         write(*,*)i,ras(i),ras2(i),ras3(i),ras4(i),ras5(i)
c      enddo

      do i=1,nCA
         r_1(1,i)=axca(i)
         r_1(2,i)=ayca(i)
         r_1(3,i)=azca(i)
         j=iLch(i)
         r_2(1,i)=x0(j)*0.87
         r_2(2,i)=y0(j)*0.87
         r_2(3,i)=z0(j)*0.87
c         write(*,*)i,j
      enddo
      call u3b(w,r_1,r_2,nCA,1,rms,u,t,ier) !u rotate r_1 to r_2
      armsd=dsqrt(rms/float(nCA)) !RMSD12
      write(20,*)'initial RMSD=',armsd

      goto 99
ccc   initial model ----------------->
      open(71,file='str',status='unknown')
      do i=1,Lch
         a1=x0(i)*0.87
         a2=y0(i)*0.87
         a3=z0(i)*0.87
         write(71,1037)i,'GLY',i,a1,a2,a3
      enddo
 1037 format('ATOM  ',i5,'  CA  ',a3,I6,4X,3F8.3)
      do i=1,Lch-1
         write(71,1239)i,i+1
      enddo
 1239 format('CONECT',I5,I5)
      close(71)

ccc   initial model ----------------->
      open(71,file='str1',status='unknown')
      do i=1,Lch
         a1=ax(i)*0.87
         a2=ay(i)*0.87
         a3=az(i)*0.87
         write(71,1037)i,'GLY',i,a1,a2,a3
      enddo
      do i=1,Lch-1
         write(71,1239)i,i+1
      enddo
      close(71)
 99   continue

      if(i_gap.gt.0)then
         write(20,*)
         write(20,*)'-----gaps in native CA------------'
         do i=1,i_gap
            write(20,*)i,s_gap(i),iLch(s_gap(i)),n_add(i)
         enddo
      else
         write(20,*)'there is no gap in native CA'
      endif
      write(20,*)

c      stop
c^^^^^^^^^^^^ read initial chain finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccc set common used parammeters cccccccccccc
      subroutine set_common
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)
      character protein*10
      common/lengths/Lch,Lch1,Lch2
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/maxdis2/maxdis2(ndim)
      common/arandom/  aarand,abrand,acrand,adrand
      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/sw1/aT_rep(nrep),E_rep(nrep)
      character*6 mname
      common/movename/mname(100)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/distres/er4,er2,es3c
      COMMON/RES/ ER3, Mcom(ndim),Kcom(ndim,50)
      common/zscore/azscore
      common/excluded/vvv(ndim,ndim)

cccccccccccccc set the random generator cccccccccccccccccccccccc
      no=random
      if(no.gt.0)no=-no
      firstrandom=aranzy(no)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Lch1=Lch-1
      Lch2=Lch-2
      anvec=nvec-0.00001
      contt=1.5*float(Lch)      !TARGET NUMBER OF CONTACTS   1.5*N
      do i=1,Lch
         maxdis2(i)=20*i*i      !maximum distance of walk in i steps
      enddo
      write(20,*)'Length:',Lch

ccccccccccccccccccccccc Temperature ccccccccccccccccccccccccccc
c     [80,130] is the standard:
      if(Lch.lt.80)then
         atemp1=atemp1*0.97
         atemp2=atemp2*0.91
         if(Lch.lt.55)then
            atemp1=atemp1*0.97
            atemp2=atemp2*0.91
         endif
      endif
      if(Lch.gt.130)then
         atemp1=atemp1*1.05
         atemp2=atemp2*1.2
         if(Lch.gt.165)then
            atemp1=atemp1*1.05
            atemp2=atemp2*1.2
            if(Lch.gt.200)then
               atemp1=atemp1*1.05
               atemp2=atemp2*1.2
               if(Lch.gt.300)then
                  atemp2=atemp2*1.2
                  if(Lch.gt.400)then
                     atemp2=atemp2*1.2
                     if(Lch.gt.500)then
                        atemp2=atemp2*1.2
                        if(Lch.gt.600)then
                           atemp2=atemp2*1.2
                           if(Lch.gt.700)then
                              atemp2=atemp2*1.2
                              if(Lch.gt.800)then
                                 atemp2=atemp2*1.2
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif

ccccccccccccc Number of replicas #################
      if(Lch.gt.165)then        !50
         N_rep=N_rep+10
      endif
      if(Lch.gt.240)then        !60
         N_rep=N_rep+10
      endif
      if(Lch.gt.300)then        !70
         N_rep=N_rep+10
      endif
      if(Lch.gt.400)then      !80
         N_rep=N_rep+10
      endif
      if(N_rep.gt.80)N_rep=80

ccccccccccccccc movement name cccccccccccccccccccccccc
      mname(2)='move2a'
      mname(3)='move3s'
      mname(4)='move3d'
      mname(5)='move4s'
      mname(6)='move4d'
      mname(7)='move8'
      mname(8)='move5s'
      mname(9)='move5d'
      mname(10)='move6'
      mname(11)='move_n'
      mname(12)='move_c'
      mname(13)='move7a'
      mname(14)='move7b'
      mname(15)='move9'
      mname(16)='tran_N'
      mname(17)='tran_M'
      mname(18)='tran_C'
      mname(19)='rot_N' !no
      mname(20)='rot_M'
      mname(21)='rot_C' !no
      mname(22)='trot_N'
      mname(23)='trot_M'
      mname(24)='trot_C'

c^^^^^^^^^^^^^^^^ common parameters finished ^^^^^^^^^^^^^^^^^^^^^^^      
      return
      end

ccccccccccccccccccc read sequence ccccccccccccccccccccccc
      subroutine read_seq
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      character*3 aa(-1:20), NAME,sequ
      common/seqe/seq(ndim),sec(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/icgg/ icg(ndim), EH6  
      common/aminoacid/sequ(ndim)

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     &     'PRO','MET','ASP','ASN','LEU',
     &     'LYS','GLU','GLN','ARG',
     &     'HIS','PHE','TYR','TRP','CYX'/
c     sequence SEQ()=0 means GLY,  =19 means TRP

      common/CA/nCA,axca(ndim),ayca(ndim),azca(ndim)
      common/initial/x0(ndim),y0(ndim),z0(ndim)
      common/addition/i_gap,s_gap(100),n_add(100),iLch(ndim),inCA(ndim)

      do i=1,Lch
         sec(i)=1   !no secondary structure
         seq(i)=1   !GLY
         sequ(i)='GLY'
      enddo

      do 121 i=1,nCA
         read(14,'i5,3x,a3,i5') k,NAME,mm
         do j=0,20
            if(NAME.eq.aa(j)) then
               m=iLch(i)        !order number on whole chain
               SEC(m)=mm
               SEQ(m)=j
               sequ(m)=name
               go to 121
            endif
         enddo
 121  continue

c^^^^^^^^^^^^^^^^^^^^^ read sequence finished ^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read centrosymmetric potential cccccccccccccccccccc
c     eoinp(A,dis) controls centrosymmetric potential of C_a;
c     eonekd(A) controls centrosymmetric potential of C_g.
c
      subroutine read_centro
      implicit integer(i-z)
      parameter(nvec=312)
      character*3 NAME
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/lengths/Lch,Lch1,Lch2
      common/hopp/eonehw(0:19)

c     read hydrophobic potential for Sg, positive for hydrophobic residue--->
      data eonekd /-0.4, 1.8, -0.8, 2.5, 4.2, -0.7, 4.5,
     &     -1.6, 1.9,  -3.5, -3.5, 3.8,
     &     -3.9, -3.5, -3.5, -4.5,
     &     -3.2, 2.8, -1.3, -0.9/
c            ^          ^     ^     !contradict with 'centro.comm'
c     read hydrophilic potential for Sg, positive for hydrophilic residue--->
      data eonehw /0.0, -0.5, 0.3, -1.0, -1.5, -0.4, -1.8,
     &     0.0, -1.3, 3.0, 0.2, -1.8,
     &     3.0, 3.0, 0.2, 3.0,
     &     -0.5, -2.5, -2.3, -3.4/

c     read centrosymmetric potential for C_a------->
      do i=0,19
         read(1,3009) name, (eoinp(i,j), j=1,4)
         eoinp(i,0)=eoinp(i,1)
         eoinp(i,5)=eoinp(i,4)+0.25 		
         do j=6,100
            eoinp(i,j)=eoinp(i,j-1)+0.25
         enddo
      enddo
 3009 FORMAT(A3,4F6.2)

c     expected gyration radius:
      acrit=2.2*exp(0.38*alog(float(Lch)))/0.87 !gyrat-radius~2.2*l^0.38
*     Defination of gyration-radius: acrit=sqrt(<(r-r0)^2>)

c^^^^^^^^^^^^^^^^^ read centersymmetric potential finished ^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read environment cccccccccccccccccccccccccccc
c         G  A  V  L  I  S  T  C  M  P  D  N  E  Q  K  R  H  F  Y  W
c (0,0,0) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (0,0,1) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (0,0,2) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c   ...                    ...
c (4,4,3) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (4,4,4) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_profile
      implicit integer(i-z)
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1

      do kkk=1,4
      do i=0,19
      do im=0,15
      do ip=0,15
      do ia=0,15
         envir(im,ip,ia,i,kkk)=2.0
      end do
      end do
      end do
      end do
      enddo

c     PROFILE3 potential =envir(#of antiparalel contacts,
c     #of orthogonal, # of parallel, aminoacid's type)
c     ia,im, ip - taken modulo 2, i.e. 0-1 contact, 2-3,...
c     profile3.comm is a score table, i.e. envir(envir_class,A)
c     here environment class is number of contacts on residue A.

      do i=0,19                 !from column
         read(2,*)
         do im=0,8
         do ia=0,8
            read(2,*)(envir(ia,im,ip,i,3),ip=0,8) !this is used
         end do
         read(2,*)
         end do
         read(2,*)
      enddo

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,1),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,2),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,4),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c^^^^^^^^^^^^^^^^^^^^^^^^ read profile finished ^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read 1-3 potential cccccccccccccccccccc
c     interaction between 1'th CA and 3'th CA
      subroutine read_E13
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      common/short2/codevsum,didevsum,csr(ndim,2)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      dimension csre(0:19,0:19,2)

c     R13 potential - two bins only (helical and expanded)
c     r2<48, E=csre(i,j,1); r2>48, E=csre(i,j,2)
      do i=0,19
         do j=0,19
            read(5,*)
            read(5,*) (csre(i,j,k),k=1,2)
         enddo
      enddo

      do i=1,Lch2
         do k=1,2
            csr(i,k)=2.0*csre(seq(i),seq(i+2),k)
         enddo
      enddo

c^^^^^^^^^^^^^^^^^ read E13 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read 1-4 potential cccccccccccccccccccc
c     interaction between 1'th CA and 4'th CA
c     the aim is to obtain IBIN(r14), asr(i,IBIN)
      subroutine read_E14
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/forpreparemove4/asrr(0:19,0:19,-12:12)
      DIMENSION asrh(0:19,0:19,-12:12),asre(0:19,0:19,-12:12)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

ccccccccc read asrr,asrh,asre------------------>
      do i=0,19                 !asrr(Ai,Bi,dis) from 'r14.comm'
         do j=0,19
            read(12,*)
            read(12,*) (asrr(i,j,k),k=-12,-5)
            read(12,*) (asrr(i,j,k),k=-4,3) !without k=4
            read(12,*) (asrr(i,j,k),k=5,12)
            do k=4,1,-1
               asrr(i,j,k)=asrr(i,j,k-1) !without k=0
            enddo
         enddo
      enddo
      do i=0,19                 !asrh(Ai,Bi,dis) from 'r14h.comm'
         do j=0,19
            read(7,*)
            read(7,*) (asrh(i,j,k),k=-12,-5)
            read(7,*) (asrh(i,j,k),k=-4,3)
            read(7,*) (asrh(i,j,k),k=5,12)
            do k=4,1,-1
               asrh(i,j,k)=asrh(i,j,k-1)
            enddo
         enddo
      enddo
      do i=0,19                 !asre(Ai,Bi,dis) from 'r14e.comm'
         do j=0,19
            read(8,*)
            read(8,*) (asre(i,j,k),k=-12,-5)
            read(8,*) (asre(i,j,k),k=-4,3)
            read(8,*) (asre(i,j,k),k=5,12)
            do k=4,1,-1
               asre(i,j,k)=asre(i,j,k-1)
            enddo
         enddo
      enddo
c^^^^^^^^^ read asrr,asrh,asre finished ^^^^^^^^^^^^^^^^^

      do i=1,Lch-3
         do k=-12,12
            asr(i,k)=asrr(seq(i+1),seq(i+2),k) !general
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
               if(sec(i).eq.2) then !helix
                  asr(i,k)=(asr(i,k)+asrh(seq(i+1),seq(i+2),k))/2.0
               endif
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
               if(sec(i).eq.4) then !sheet
                  asr(i,k)=(asr(i,k)+asre(seq(i+1),seq(i+2),k))/2.0
               endif
            endif
         enddo
      enddo
c^^^^^^^^^^^^ asr(i,ibin(r14)) finished ^^^^^^^^^^^^^^^^^^^^^
c     r(i,i+3)=k, E=asr(i,k), 12 bins (24 bins when considering chiral)
      do i=1,300
         kk=int((sqrt(float(i))*0.87))+1
         if(kk.gt.12) kk=12
         IBIN(I) = kk           !convert lattice r^2 into real r
         IBIN(-I)=-kk
      ENDDO
      IBIN(0)=IBIN(1)

c^^^^^^^^^^^^^^^^^ read E14 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read 1-5 potential cccccccccccccccccccc
c     interaction between 1'th CA and 5'th CA
c     the aim is to obtain JBIN(r15), bsr(i,JBIN)
      subroutine read_E15
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      DIMENSION bsrh(0:19,0:19,16)
      dimension bsre(0:19,0:19,16)
      dimension bsrr(0:19,0:19,16)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

cccccccc read bsrr,bsrh,bsre ----------------->
      do i=0,19                 !read bsrr from 'r15.dat'
         do j=0,19
            read(9,*)
            read(9,*) (bsrr(i,j,k),k=1,8)
            read(9,*) (bsrr(i,j,k),k=9,16)
         enddo
      enddo
      do i=0,19                 !read bsrh from 'r15h.dat'
         do j=0,19
            read(10,*)
            read(10,*) (bsrh(i,j,k),k=1,8)
            read(10,*) (bsrh(i,j,k),k=9,16)
         enddo
      enddo	
      do i=0,19                 !read bsre from 'r15e.dat'
         do j=0,19
            read(11,*)
            read(11,*) (bsre(i,j,k),k=1,8)
            read(11,*) (bsre(i,j,k),k=9,16)
         enddo
      enddo	

      do i=1,Lch-4
         do k=1,16
            bsr(i,k)=bsrr(seq(i+1),seq(i+3),k)
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
               bsr(i,k)=(bsr(i,k)+bsrh(seq(i+1),seq(i+3),k))/2.0 !helix
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
               bsr(i,k)=(bsr(i,k)+bsre(seq(i+1),seq(i+3),k))/2.0 !sheet
            endif
         enddo
      enddo
c^^^^^^^^^^^^^^^^^^^^ E_15(Ai,Aj,dis) prepared ^^^^^^^^^^^^^^^^^^^

c     prepare distance bin-------------------->
      do i=0,500
         kk=int((sqrt(float(i))*0.87))+1 !i, lattice-dist; kk, real distance
         if(kk.gt.16) kk=16
         JBIN(I) = kk           !jbin: real distance
      ENDDO

ccccc acops(i,jbin) to enhance the contacts between gragments cccc
      do i=1,Lch-4
         acops(i,1)=(min(bsr(i,1),0.0))/2.0 !acpos<0
         do k=2,15
            acops(i,k)=min(0.0,bsr(i,k-1)+2.0*bsr(i,k)+bsr(i,k+1)) !<0
         enddo
         acops(i,16)=(min(bsr(i,16),0.0))/2.0
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c^^^^^^^^^^^^^^^^^ read E15 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read contact-pair potential cccccccccccccccccccc
c     potential: app(Ai,Aj), apa(Ai,Aj), apm(Ai,Aj)
c     distance range: [arlp,alp], [arla,ala], [arlm,alm]
c     all the general data in 'quarsi3.comm'.
c     the sequence-dependent contact-pair data in 'pair3.dat'.
      subroutine read_pair
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      character*3 NAME
      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/pair1/eh2,eh1b
      common/pairmap1/apba(ndim,ndim),apbp(ndim,ndim),apbm(ndim,ndim)
      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap
      common/paircut/ash

c     Pairwise interactions apablp ... and cut-off parmeters
c     arlp, orientation dependent, pairwise specific, sequence
c     independent

c     read contact-pair potential from 'quarsi3.comm' ------->
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apablp(i,j),j=0,19) !for app
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apablm(i,j),j=0,19) !for apm
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apabla(i,j),j=0,19) !for apa
      enddo
c     read distance-range from 'quarsi3.comm' ------->
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arlp(i,j),j=0,19) !max distance for parallel
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arlm(i,j),j=0,19) !for perpendicular contact
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arla(i,j),j=0,19) !for antiparellel pair
      enddo
 725  format(a3,1x,20f5.1)

c     width of energy-well: [arla,ala]
      ash=0.17
      ash_min=1-ash
      ash_max=1+ash
      do i=0,19
         do j=0,19
            ala(i,j)=(arla(i,j)*ash_max/0.87)**2
            alm(i,j)=(arlm(i,j)*ash_max/0.87)**2
            alp(i,j)=(arlp(i,j)*ash_max/0.87)**2
            arla(i,j)=(arla(i,j)*ash_min/0.87)**2
            arlm(i,j)=(arlm(i,j)*ash_min/0.87)**2
            arlp(i,j)=(arlp(i,j)*ash_min/0.87)**2
         enddo
      enddo
c     E=EH1/2, for r in [0,arlp];
c     E=app-es*fs,  for [0,alp];
c     E=0,     for r in [alp,00].
c^^^^^^^^^^^^^contact interaction range finished ^^^^^^^^^^^^^^^^^

c     combine the data in 'quarsi3.comm' and 'pair3.dat' to get
c     contact potential-------------------->
      do i=1,Lch
         do j=1,Lch
            apa(i,j)=apabla(ii,jj)
            apm(i,j)=apablm(ii,jj)
            app(i,j)=apablp(ii,jj)
         enddo
      enddo
c^^^^^^^^^^^^^^^ pair-potential is obtained ^^^^^^^^^^^^^^^^^^^^^^^^^

c^^^^^^^^^^^^^^^^^ read contact-pair finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read cut-off for contact predictions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_concut
      implicit integer(i-z)
      character*3 NAME
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      common/zscore/azscore

      read(25,*)
      do i=0,19
         read(25,*)NAME,(concut(i,j),j=0,19)
      enddo

      do i=0,19
         do j=0,19
            concut(i,j)=concut(i,j)/100 !real cut-off
            concut(i,j)=concut(i,j)/0.87 !cut-off on lattice
            if(azscore.lt.10)then
               concut(i,j)=concut(i,j)*0.8
            endif
            concut2(i,j)=concut(i,j)**2 !cut-off squared on lattice
         enddo
      enddo

      return
      end

cccccccccccccccc read contact restrains from 'comb.dat' cccccccccc
c     aim: Mcom(i)--->number of contacts on i-residue
c          Kcom(i,ii): ii's contact of i-residue is j=kres
c     contact restrain is only on C_g.
      subroutine read_contactrestrain
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/lengths/Lch,Lch1,Lch2
      common/distres/er4,er2,es3c
      common/res/ er3, Mcom(ndim),Kcom(ndim,50)
      common/resnumber/Ncom,Ndis,accur
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b

      common/addition/i_gap,s_gap(100),n_add(100),iLch(ndim),inCA(ndim)
      DIMENSION r1(4000),r2(4000)

c     READS Side group - side group contacts 
c     (from NMR or therading predictions or clusters)

ccc   pool all the restraints into r1(i),r2(i)--------------->
      read(16,*)ntmp
      i_c=0
      do i=1,ntmp
         i_c=i_c+1
         read(16,*)i1,i2
         r1(i_c)=iLch(i1)
         r2(i_c)=iLch(i2)
      enddo
      Ncom=i_c

ccc   map r1,2(i) into Mcom(i),Kcom(i,Mcom(i))------------>
      do i=1,Lch
         Mcom(i)=0              !number of contacts with 'i'
         do j=1,Ncom
            if(r1(j).eq.i)then
               Mcom(i)=Mcom(i)+1
               Kcom(i,Mcom(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               Mcom(i)=Mcom(i)+1
               Kcom(i,Mcom(i))=r1(j)
            endif
         enddo
      enddo

      colim=1.5*Ncom           !background number for derviation
c     the larger 'colim' is, the weaker the contact restrain is.

ccc   output restraints------------->
      write(20,*)'Number of restraints:',Ncom
      write(20,*)'----------- contact restraints ---------------'
      nnc=0
      do i=1,Lch
         nnc=nnc+Mcom(i)
         write(20,12)i,Mcom(i),(Kcom(i,j),j=1,Mcom(i))
 12      format(i4,'(',i2,'):',20i4)
      enddo
      write(20,*)'Number of contact=',nnc,' Lch=',Lch
      write(20,*)'fc=',float(nnc)/Lch
      write(20,*)

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read distance restrains from 'dist.dat' cccccccccc
      subroutine read_distantrestrain
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nvec=312)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/lengths/Lch,Lch1,Lch2
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/distres/er4,er2,es3c
      common/resnumber/Ncom,Ndis,accur
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
      common/addition/i_gap,s_gap(100),n_add(100),iLch(ndim),inCA(ndim)

      dimension r1(10000),r2(10000),dis(10000),deviation(10000)

      read(17,*)Ndis
      do i=1,Ndis
         read(17,*)i1,i2,nothing,dis(i)
         deviation(i)=0
         r1(i)=iLch(i1)
         r2(i)=iLch(i2)
         dis(i)=dis(i)/0.87
         deviation(i)=deviation(i)/0.87
         deviation(i)=1+deviation(i)**2
      enddo

      do i=1,Lch
         Mdis(i)=0              !number of prediction for 'i'
         do j=1,Ndis
            if(r1(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r2(j) !r2(j) with 'i'
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
            if(r2(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r1(j)
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
         enddo
      enddo
      dilim=1.5*Ndis          !background number for derviation

c      write(20,*)
c      write(20,*)'----------- distant map ---------------'
c      nnc=0
c      do i=1,Lch
c         nnc=nnc+Mdis(i)
c         write(20,12)i,Mdis(i),(Kdis(i,j),dist(i,j)*0.87,j=1,Mdis(i))
c 12      format(i4,':',i3,20(i4,'-'f5.2))
c      enddo
c      write(20,*)'Number of distmap=',nnc,' Lch=',Lch

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccset temperature for different replicas ccccc
      subroutine set_temperature
      implicit integer(i-z)
      parameter(nrep=100)
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)

***********for normal run ***********************************
      do i=1,N_REP
         aT_rep(i)=atemp1*(atemp2/atemp1)**(float(i-1)/(N_rep-1))
c         write(*,*)i,aT_rep(i)
      enddo

c      stop
c^^^^^^^^^^^^^^^^^ set aT_rep(i) finished ^^^^^^^^^^^^^^^
      return
      end

cccccccccccccc set EHBIJ(i,j) ccccccccccccccccccc
      subroutine set_EHB
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nvec=312)
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)

c
c     EHBIJ - set-up secondary structure dependent
c     strength of the hyrogen bond network - stronger for helices
c     and beta-sheets
c

      do i=1,Lch
         is=sec(i)
         do j=1,Lch
            js=sec(j)
            EHBIJ(i,j)=1
            if(iabs(i-j).eq.3.and.is.eq.2.and.js.eq.2)then
               EHBIJ(i,j)=EHBIJ(i,j)+0.5 !helix, enhanced
            endif
            if(is.eq.4.or.js.eq.4) then
               if(is*js.ne.8.and.iabs(i-j).gt.4)then
                  EHBIJ(i,j)=EHBIJ(i,j)+0.5 !beta-beta, enhanced
               endif
            endif
         enddo
      enddo

c^^^^^^^^^^^ set H_bond finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c define hydrogen-bond, C_beta, C_group for all possible neighbors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_beta
      implicit integer(i-z)
                parameter(nvec=312)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/three/angle(nvec,nvec)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)


      do k=0,19
         read(4,*)
         read(4,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
      enddo
      CLOSE(4)                  !sidecent.comm

ccccccccccc define hydrogen-bond, C_beta, C_group for good (i,j)----->
      do 101 i=1,nvec
      do 102 j=1,nvec
         avi=sqrt(float(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)))
         avxi=vx(i)/avi
         avyi=vy(i)/avi
         avzi=vz(i)/avi
         avj=sqrt(float(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)))
         avxj=vx(j)/avj
         avyj=vy(j)/avj
         avzj=vz(j)/avj

         ax=avxi+avxj
         ay=avyi+avyj
         az=avzi+avzj
         aaa=sqrt(ax*ax+ay*ay+az*az)
         ax=ax/aaa                !(vi+vj)/|vi+vj|
         ay=ay/aaa                !(vi+vj)/|vi+vj|
         az=az/aaa                !(vi+vj)/|vi+vj|

         bx=avyi*avzj-avzi*avyj
         by=avzi*avxj-avxi*avzj
         bz=avxi*avyj-avyi*avxj
         bbb=sqrt(bx*bx+by*by+bz*bz)
         bx=bx/bbb             ! vi(x)vj/|vi(x)vj|
         by=by/bbb             ! vi(x)vj/|vi(x)vj|
         bz=bz/bbb             ! vi(x)vj/|vi(x)vj|

c         write(*,*)'ini bbb=',bbb
c         write(*,*)bx**2+by**2+bz**2

         cx=avxi-avxj
         cy=avyi-avyj
         cz=avzi-avzj
         ccc=sqrt(cx*cx+cy*cy+cz*cz)
         cx=cx/ccc              !(vi-vj)/|vi-vj|
         cy=cy/ccc              !(vi-vj)/|vi-vj|
         cz=cz/ccc              !(vi-vj)/|vi-vj|
         cax(i,j)=cx            !(vi-vj)/|vi-vj|
         cay(i,j)=cy            !(vi-vj)/|vi-vj|
         caz(i,j)=cz            !(vi-vj)/|vi-vj|

c     H-bond (unit vector):
         hbx(i,j)=bx
         hby(i,j)=by
         hbz(i,j)=bz

c     side-chain coordinate from C_a to side-chain ---------------->
         do k=0,19
            if(angle(i,j).lt.105) then ! alpha-helix or turn like
               gx(i,j,k)=(axalf(k)*ax+ayalf(k)*bx+azalf(k)*cx)/0.87
               gy(i,j,k)=(axalf(k)*ay+ayalf(k)*by+azalf(k)*cy)/0.87
               gz(i,j,k)=(axalf(k)*az+ayalf(k)*bz+azalf(k)*cz)/0.87
            else                ! beta-sheet
               gx(i,j,k)=(axbet(k)*ax+aybet(k)*bx+azbet(k)*cx)/0.87
               gy(i,j,k)=(axbet(k)*ay+aybet(k)*by+azbet(k)*cy)/0.87
               gz(i,j,k)=(axbet(k)*az+aybet(k)*bz+azbet(k)*cz)/0.87
            endif
         enddo
 102  continue
 101  continue

      return
      end

cccccccccc Compute the secondary fragment biases cccccccccccccccc
c     check local secondary structure from 'seq.dat'
c     if it is beta-structure in [i,i+6],  frga(i)=19.1/0.87;
c     if it is alpha-structure in [i,i+7], frgb(i)=10.5/0.87;
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_frg
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      common/fr/frga(ndim),frgb(ndim)
      common/seqe/seq(ndim),sec(ndim)
      common/lengths/Lch,Lch1,Lch2

      do i=1,Lch
         frga(i)=0.0
         frgb(i)=0.0
      enddo

      do i=1,Lch-7
         q=0
         do j=i,i+7
            if(sec(j).eq.2) q=q+2 !helix structure.
         enddo
         if(q.eq.16)then        !8 continue alpha-residues
            frga(i)=10.5/0.87   !distance for 7 alpha-bonds
         endif
      enddo

      do i=1,Lch-6
         q=0
         do j=i+1,i+5
            if(sec(j).eq.4) q=q+4 !beta structure
         enddo
         if(q.eq.20)then        !5 continue beta-residues
            if(sec(i).ne.2.and.sec(i+6).ne.2)then
               frgb(i)=19.1/0.87 !distance for 6 beta-bonds
            endif
         endif
      enddo

c      do i=1,Lch
c         write(*,*)i,seq(i),sec(i),frga(i),frgb(i)
c      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prepare 2-bond move, i.e. calculate v21(tx,ty,tz,i), v22(tx,ty,tz,i).
c v21(tx,ty,tz,i) --- 1th vector of i'th path from (0,0,0) to (tx,ty,tz)
c v22(tx,ty,tz,i) --- 2th vector of i'th path from (0,0,0) to (tx,ty,tz)
c it will be dangerous if number of used variable is more than 3,000,000
c i.e. the usage of memory of CUP can not beyond 90%.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_move2
      implicit integer(i-z)
                parameter(nvec=312)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,46)
      common/w2b/w22(-10:10,-10:10,-10:10,46)

      max=0
      nnn=0
      mmm=0
      do 101 i=1,nvec
         do 102 j=1,nvec
            if(goodc(i,j))then
               ijx=vx(i)+vx(j)
               ijy=vy(i)+vy(j)
               ijz=vz(i)+vz(j)
               Np2(i,j)=0
               do 103 ii=1,nvec
                  jx=ijx-vx(ii)
                  jy=ijy-vy(ii)
                  jz=ijz-vz(ii)
                  jr=jx*jx+jy*jy+jz*jz
                  if(jr.ge.14.and.jr.le.25)then
                     jj=vector(jx,jy,jz)
                     if(goodc(ii,jj))then
                        Np2(i,j)=Np2(i,j)+1 !based on i,j
                        v21(i,j,Np2(i,j))=ii
                        v22(i,j,Np2(i,j))=jj
                        mmm=mmm+1 !number of total memory occupied by vx21.
                     endif
                  endif
 103           continue
               if(max.lt.Np2(i,j))max=Np2(i,j)
               nnn=nnn+1        !number of possible pairs
            endif
 102     continue
 101  continue

ccc among all 312*312=97344 pairs, 67272 pairs are legal (~70%).
ccc all Np2(i,j)>=2, i.e. there are at least one other path for any pair.
ccc <Np2(i,j)>=27.

      write(*,*)'maximum of Np2(i,j)=',max
      write(*,*)'number of possible pair=',nnn
      write(*,*)'sum of Np2(i,j), total memory=',mmm
cccc  the following is the cases without goodc on prefabracated pairs:
c              N_v  N_pare N_pp   mmm
c     [14,25], 312, 97344, 67272, 1830000 (25% memory, within memory)

cccc the following is cases without goodc limitation:
c     [14,25], 306, 93636, 3872214 (90% memory, beyond memory)
c     [14,24], 282, 79524, 2923758 (80% memory, beyond memory)
c      stop

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prepare 3-bond move, i.e. calculate v31(tx,ty,tz,i), v32(tx,ty,tz,i).
c v31(tx,ty,tz,i) --- 1th vector of i'th path from (0,0,0) to (tx,ty,tz)
c v32(tx,ty,tz,i) --- 2th vector of i'th path from (0,0,0) to (tx,ty,tz)
c v33(tx,ty,tz,i) --- 3th vector of i'th path from (0,0,0) to (tx,ty,tz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_move3
      implicit integer(i-z)
      parameter(nvec=312)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/move31/v31(-15:15,-15:15,-15:15,6)
      common/move32/v32(-15:15,-15:15,-15:15,6)
      common/move33/v33(-15:15,-15:15,-15:15,6)
      common/move34/Np3(-15:15,-15:15,-15:15)

c      dimension v31(-15:15,-15:15,-15:15,100)
c      dimension v32(-15:15,-15:15,-15:15,100)
c      dimension v33(-15:15,-15:15,-15:15,100)
c      dimension Np3(-15:15,-15:15,-15:15)

      do i=-15,15
         do j=-15,15
            do k=-15,15
               Np3(i,j,k)=0
            enddo
         enddo
      enddo


c      write(*,*)'1111111'
      num3=0
      max=0
      do 101 i=1,nvec
         do 102 j=1,nvec
            if(goodc(i,j))then
               do 103 k=1,nvec
                  if(goodc(j,k))then
                     num3=num3+1
                     rx=vx(i)+vx(j)+vx(k)
                     ry=vy(i)+vy(j)+vy(k)
                     rz=vz(i)+vz(j)+vz(k)
                     Np3(rx,ry,rz)=Np3(rx,ry,rz)+1
                     v31(rx,ry,rz,Np3(rx,ry,rz))=i
                     v32(rx,ry,rz,Np3(rx,ry,rz))=j
                     v33(rx,ry,rz,Np3(rx,ry,rz))=k
                     if(max.le.Np3(rx,ry,rz)) max=Np3(rx,ry,rz)
                     write(*,*)i,j,k,max
                  endif
 103           continue
            endif
 102     continue
 101  continue

      write(*,*)'number of move3=',num3,max

      stop
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c calculate contact order according to the length and secondary structure
c	  SIZE     H       E       H/E
c  	2   40  0.116   0.324      0.252
c   	3   60  0.119   0.357      0.230
c  	4   80  0.115   0.280      0.212
c   	5  100  0.105   0.259      0.198
c   	6  120  0.132   0.269      0.168
c   	7  140  0.105   0.272      0.176
c   	8  160  0.114   0.186      0.183
c   	9  180  0.116   0.197      0.160
c      10  200  0.107   0.184      0.134
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_acorder
      implicit integer(i-z)
		parameter(ndim=1000)
                parameter(nvec=312)
      common/seqe/seq(ndim),sec(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/lengths/Lch,Lch1,Lch2
      
c     data struct /' coil','helix',' turn',' beta','    ?'/

************* number of predicted structures *************
      n_H=0                     !number of Helix
      n_E=0                     !number of Extension
      do i=1,Lch
         if(sec(i).eq.2) n_H=n_H+1
         if(sec(i).eq.4) n_E=n_E+1
      enddo

      if(n_H+n_E.lt.2)then      !use alpha=H/E
         if(Lch.lt.50)then
            alph=0.252
         else if(Lch.lt.70)then
            alph=0.230
         else if(Lch.lt.90)then
            alph=0.212
         else if(Lch.lt.110)then
            alph=0.198
         else if(Lch.lt.130)then
            alph=0.168
         else if(Lch.lt.150)then
            alph=0.176
         else if(Lch.lt.170)then
            alph=0.183
         else if(Lch.lt.190)then
            alph=0.160
         else
            alph=0.134
         endif
      else                      !use alpha=aH+bE
         a1=float(n_H)/float(n_H+n_E)
         a2=float(n_E)/float(n_H+n_E)
         if(Lch.lt.50)then
            alph=0.116*a1+0.324*a2
         else if(Lch.lt.70)then
            alph=0.119*a1+0.357*a2
         else if(Lch.lt.90)then
            alph=0.115*a1+0.280*a2
         else if(Lch.lt.110)then
            alph=0.105*a1+0.259*a2
         else if(Lch.lt.130)then
            alph=0.132*a1+0.269*a2
         else if(Lch.lt.150)then
            alph=0.105*a1+0.272*a2
         else if(Lch.lt.170)then
            alph=0.114*a1+0.186*a2
         else if(Lch.lt.190)then
            alph=0.116*a1+0.197*a2
         else
            alph=0.107*a1+0.184*a2
         endif
      endif

      acorder=alph*Lch
c^^^^^^^^^^^^^^^^^^ contact order done ^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   produce initial structures randomly:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine random_initial
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
                parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension ip(ndim)
      common/sw3/icarep(ndim,nrep)

      do 102 k=1,N_REP
 88      x(1)=0
         y(1)=0
         z(1)=0
         m=0
         do 101 i=2,Lch
 99         ip(i)=int(aranzy(no)*nvec)+1
            m=m+1
            if(m.gt.1000000)then
               write(*,*) 'UNSOLVABLE STERIC PROBLEM > EXIT_2'
               goto 88          !unsolvable steric problem
            endif
            if(i.gt.2)then      !check neighbor
               if(.not.goodc(ip(i-1),ip(i)))goto 99
            endif
            x(i)=x(i-1)+vx(ip(i))
            y(i)=y(i-1)+vy(ip(i))
            z(i)=z(i-1)+vz(ip(i))
            do j=1,i-1          !check excluded volumn for Ca
               ir=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
               if(ir.lt.14) goto 99
            enddo
 101     continue

         do i=1,Lch
            xrep(i,k)=x(i)
            yrep(i,k)=y(i)
            zrep(i,k)=z(i)
         enddo
 102  continue

c^^^^^^^^^^^^ initial chains finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     get q(i), cx0(i) for consensus segments of top-2 templates.
c     if can not find consensus segment, back to top-1.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_consensus
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/consensus1/cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)
      dimension cx1(1000),cy1(1000),cz1(1000),q1(1000),ip1(1000)
      dimension cx2(1000),cy2(1000),cz2(1000),q2(1000),ip2(1000)
      dimension qq(1000)
      real xx,yy,zz
ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /1000*1.0/
ccc   

******************************************************************
***   read templates------------->
      rewind(24)
      read(24,*)N_tmp
      if(N_tmp.lt.2)then
         write(20,*)'There is only one template'
         write(20,*)'Top-1 template is used'
         return                 !use the first template
      else
***   read first template:
         read(24,*)N_al
         do i=1,N_al
            read(24,1237)text,ii,text,a1,a2,a3
            cx1(ii)=a1
            cy1(ii)=a2
            cz1(ii)=a3
            q1(ii)=1
         enddo
         read(24,*)text
***   read second template:
         read(24,*)N_al
         do i=1,N_al
            read(24,1237)text,ii,text,a1,a2,a3
            cx2(ii)=a1
            cy2(ii)=a2
            cz2(ii)=a3
            q2(ii)=1
         enddo
         read(24,*)text
***
      endif
 1237 format(A22,I4,A4,3F8.3)

******************************************************************
***   decided qq(i)------------->
      do i=1,Lch
         qq(i)=0
      enddo
      k=0
      do i=1,Lch
         if(q1(i).eq.1.and.q2(i).eq.1)then
            k=k+1
            ip1(k)=i
            r_1(1,k)=cx1(i)
            r_1(2,k)=cy1(i)
            r_1(3,k)=cz1(i)
            r_2(1,k)=cx2(i)
            r_2(2,k)=cy2(i)
            r_2(3,k)=cz2(i)
         endif
      enddo
      write(*,*)"#common aligned=",k
      if(k.lt.10)then
         write(20,*)'There is less than 10 common aligned residues'
         write(20,*)'Top-1 template is used'
         return                 !no common aligned points, using template1
      endif
      call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2
      armsd=dsqrt(rms/k)        !RMSD12
      write(20,*)'RMSD1=',armsd,k !RMSD between template1,2 for common align
      kk=0
      do j=1,k
         xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
         yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
         zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
         dis=sqrt((xx-r_2(1,j))**2+(yy-r_2(2,j))**2+(zz-r_2(3,j))**2)
c         write(*,*)j,ip1(j),dis
         if(dis.lt.5)then
            kk=kk+1
            qq(ip1(j))=1
         endif
      enddo
      if(kk.lt.10)then
         write(20,*)'There is less than 10 close common residues'
         write(20,*)'Top-1 template is used'
         return                 !no consensus points, using template1
      endif

******************************************************************
***   q(i)=qq(i), cx0(i)=(cx1+cx2)/2--------------->
      do i=1,Lch
         q(i)=qq(i)
         cx0(i)=1000000.        !for checking excluded volumn
         cy0(i)=1000000.
         cz0(i)=1000000.
      enddo
      k=0
      do i=1,Lch
         if(q(i).eq.1)then
            k=k+1
            ip2(k)=i
            r_1(1,k)=cx1(i)
            r_1(2,k)=cy1(i)
            r_1(3,k)=cz1(i)
            r_2(1,k)=cx2(i)
            r_2(2,k)=cy2(i)
            r_2(3,k)=cz2(i)
         endif
      enddo
      call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,k
         xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
         yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
         zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
         cx0(ip2(j))=(xx+r_2(1,j))/2
         cy0(ip2(j))=(yy+r_2(2,j))/2
         cz0(ip2(j))=(zz+r_2(3,j))/2
      enddo
********************************************************************
      L_cut=2                   !small segment to be shrowed away.

c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sort ras, recalculate nfl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sort_ras
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/chain0/ras(ndim),nfl

      dimension ras0(ndim)

      ras_min=10000
      ras_max=-10000
      do i=1,nfl
         if(ras(i).gt.ras_max)ras_max=ras(i)
         if(ras(i).lt.ras_min)ras_min=ras(i)
      enddo

      nfl0=1
      ras0(nfl0)=ras_min
      do 1 while(ras0(nfl0).lt.ras_max)
         ras_min=10000
         do i=1,nfl
            if(ras(i).gt.ras0(nfl0))then
               if(ras(i).lt.ras_min)then
                  ras_min=ras(i)
               endif
            endif
         enddo
         nfl0=nfl0+1
         ras0(nfl0)=ras_min
 1    continue

      nfl=nfl0
      do i=1,nfl
         ras(i)=ras0(i)
      enddo

*^^^^^^^^^^^^^^ sort ras done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     (x,y,z) --> (r,thita,phi)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine triangle(bx,by,bz,br,bthita,bphi)
      parameter(bpi=3.1415926)

      br=sqrt(bx*bx+by*by+bz*bz)
      bthita=acos(bz/br)
      if(abs(bx).gt.0.00001)then
         bphi=atan(abs(by/bx))
      else
         bphi=0
      endif
      if(bx.gt.0)then
         if(by.lt.0)bphi=bphi+bpi*1.5
      else
         if(by.gt.0)then
            bphi=bphi+bpi*0.5
         else
            bphi=bphi+bpi
         endif
      endif

c^^^^^^^^^^^^^^^^^^ triangle done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccc decide movement point from ras(i) cccccccccccccccccccc
c     nfl----total number of moveable points on model-chain
c     nfl2---total number of 2-bond-movement
c     nfl3---total number of 3-bond-movement
c     nfl4---total number of 4-bond-movement
c     nfl5---total number of 5-bond-movement
c     ras(i)---position of i-th moveable residues
c     ras2(i)--initial position of i-th 2-bond movement.
c     ras3(i)--initial position of i-th 3-bond movement.
c     ras4(i)--initial position of i-th 4-bond movement.
c     ras5(i)--initial position of i-th 5-bond movement.
c     ras6(i)--initial position of i-th 6-bond movement.
cccccccccccccc
      subroutine move_point
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/chainm/mv(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

      do i=1,Lch+10
         mv(i)=-1               !frozen point
      enddo

      do i=1,nfl
         mv(ras(i))=1           !moveable point
      enddo

c     re-decide the movement range of tremendicy:
      Mend_N=0                  !Mend_N points from 1 are moveable
      k=0
      do i=1,Lch
         k=k+1
         if(k.gt.Mend.or.mv(i).lt.0) goto 111
         Mend_N=Mend_N+1
      enddo
 111  continue
      Mend_C=0                  !Mend_C points from Lch are moveable
      k=0
      do i=Lch,1,-1
         k=k+1
         if(k.gt.Mend.or.mv(i).lt.0) goto 222
         Mend_C=Mend_C+1
      enddo
 222  continue

c     find 6-bond-move position -------------->
      nfl6=0                    !total number of 6-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then
                     if(mv(i+4).gt.0)then
                        if(mv(i+5).gt.0)then
                           if(mv(i+6).gt.0)then !i+6: fixed border
                              nfl6=nfl6+1
                              ras6(nfl6)=i
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 5-bond-move position -------------->
      nfl5=0                    !total number of 5-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then
                     if(mv(i+4).gt.0)then
                        if(mv(i+5).gt.0)then !i+5: fixed border
                           nfl5=nfl5+1
                           ras5(nfl5)=i
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 4-bond-move position -------------->
      nfl4=0                    !total number of 4-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then
                     if(mv(i+4).gt.0)then !i+4: fixed border
                        nfl4=nfl4+1
                        ras4(nfl4)=i
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 3-bond-move position -------------->
      nfl3=0                    !total number of 3-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then !i+3 fixed border
                     nfl3=nfl3+1
                     ras3(nfl3)=i
                  endif
               endif
            endif
         endif
      enddo

c     find 2-bond-move position -------------->
      nfl2=0                    !total number of 2-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then !i+2 fixed border
                  nfl2=nfl2+1
                  ras2(nfl2)=i
               endif
            endif
         endif
      enddo

c^^^^^^^^^^^^^^^^^^^^ move_point finished ^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc connect gaps in [i1,i2] for (bx,by,bz):
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine connect(i1,i2,pass)
      implicit integer(i-z)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/initialstr/cx(1000),cy(1000),cz(1000)

      pass=1
      if(i1.eq.1)then           !!!!!N-terminal or whole structure, random walk
         do i=i2,1,-1
            n_check=0
 10         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i+1)+axx
            cy(i)=cy(i+1)+ayy
            cz(i)=cz(i+1)+azz
            n_check=n_check+1
            if(n_check.gt.100000)then
               pass=3
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 10 !excluded volumn
         enddo
      elseif(i2.eq.Lch)then     !!!!!!!!!!!!!!!!!!C_terminal,
         do i=i1,Lch
            n_check=0
 11         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i-1)+axx
            cy(i)=cy(i-1)+ayy
            cz(i)=cz(i-1)+azz
            n_check=n_check+1
            if(n_check.gt.100000)then
               pass=4
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 11
         enddo
      else                      !!!!!!!!!!!!!interval
         diss=di(cx(i2+1),cy(i2+1),cz(i2+1),cx(i1-1),cy(i1-1),cz(i1-1))
         adis=diss/float((i2+1)-(i1-1))
***   linear connect for big gap---->
         if(adis.ge.3.5)then
            dex=3.5*(cx(i1-1)-cx(i2+1))/diss !3.5*cos(thita)
            dey=3.5*(cy(i1-1)-cy(i2+1))/diss
            dez=3.5*(cz(i1-1)-cz(i2+1))/diss
            do j=i2,i1,-1
               cx(j)=cx(j+1)+dex
               cy(j)=cy(j+1)+dey
               cz(j)=cz(j+1)+dez
            enddo
            return              !end of connection
         endif
***   random walk from i1 to i2--------->
         bdis0=3.5
         m_check=0              !try 2000 times of whole walk
 13      m_check=m_check+1
         do 14 i=i1,i2
            n_check=0           !each point try 2000 times
 12         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i-1)+axx
            cy(i)=cy(i-1)+ayy
            cz(i)=cz(i-1)+azz
            n_check=n_check+1
            if(n_check.gt.1000)bdis0=3.7
            if(n_check.gt.2000)goto 13
            if(m_check.ge.2000)then !can not pass the connection
               pass=5
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 12 !check excluded V
            aaa=float(i2+1-i)
            bdis=di(cx(i),cy(i),cz(i),cx(i2+1),cy(i2+1),cz(i2+1))/aaa
            if(i.lt.i2.and.bdis.ge.bdis0) goto 12 !check remain steps
            if(i.eq.i2)then
               if(bdis.gt.4.2.or.bdis.lt.3.4) goto 12 !last step
            endif
            bdis0=3.5
 14      continue
      endif

c^^^^^^^^^^^^^^^^ connect of gap [i1,i2] finished ^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  produce a vector of length=3.8A:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_bond(axx,ayy,azz,al)
      implicit integer(i-z)
      athita=acos(1.-2.*aranzy(no)) !thita angle in random, [0,pi]
      aphi=2.*3.1415926*aranzy(no) !phi angle in random, [0,2pi]
      axx=al*sin(athita)*cos(aphi)
      ayy=al*sin(athita)*sin(aphi)
      azz=al*cos(athita)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  check weak excluded volumn for cx():
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function mcheck(i0,bx0,by0,bz0)
      implicit integer(i-z)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/initialstr/cx(1000),cy(1000),cz(1000)

      mcheck=1
      do i=1,Lch
         if(i0.ne.i)then
            dis=di(bx0,by0,bz0,cx(i),cy(i),cz(i)) !distance
            if(dis.le.3.1)then
               mcheck=3
               return
            endif
         endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function di(x1,y1,z1,x2,y2,z2)
      di=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function di2(x1,y1,z1,x2,y2,z2)
      di2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function aax(i)
      implicit integer(i-z)
      parameter(ndim=1000)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aax=x(i)
      else
         aax=ex(i)
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function aay(i)
      implicit integer(i-z)
      parameter(ndim=1000)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aay=y(i)
      else
         aay=ey(i)
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function aaz(i)
      implicit integer(i-z)
      parameter(ndim=1000)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aaz=z(i)
      else
         aaz=ez(i)
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     for output of trajectory
c     i-->residue;
c     k-->replica;
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_fxyz(k)
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      common/chainm/mv(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
      common/outputxyz/fxyz(3,ndim)

      do i=1,Lch
         if(mv(i).gt.0)then
            fxyz(1,i)=xrep(i,k)
            fxyz(2,i)=yrep(i,k)
            fxyz(3,i)=zrep(i,k)
         else
            fxyz(1,i)=exrep(i,k)
            fxyz(2,i)=eyrep(i,k)
            fxyz(3,i)=ezrep(i,k)
         endif
      enddo

c^^^^^^^^^^^^^^^^ fxyz finished ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccc print out initial parameters cccccccccccccccccc
      subroutine write_parameter
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nvec=312)
      parameter(nrep=100)       !number of replicas
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/lengths/Lch,Lch1,Lch2
      COMMON/RES/ ER3, Mcom(ndim),Kcom(ndim,50)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/distres/er4,er2,es3c
      common/resnumber/Ncom,Ndis,accur
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/looks/exc,exc1,exc2

      write(20,*)'Protein real length..................',Lch
      write(20,*)
      write(20,*)'Ncycle...............................',ncycle
      write(20,*)
      write(20,*)'number of N_rep......................',N_REP
      write(20,*)'N_swap...............................',phot
      write(20,*)'ncycle*N_swap........................',ncycle*phot
      write(20,*)'local moves each replica:',ncycle*phot*float(Lch)
      write(20,*)'total moves in MC..',ncycle*phot*float(Lch)*N_rep
      write(20,*) 
      write(20,*)'.......N_dist........................',Ndis
      write(20,*)'.......N_comb........................',Ncom
      write(20,*)'.......N_dist/length........',Ndis/float(Lch)
      write(20,*)'.......N_comb/length........',Ncom/float(Lch)
      write(20,*) 
      write(20,*)'maximum temperture...................',atemp2
      write(20,*)'minimum temperture...................',atemp1
      write(20,*)
      write(20,*)'Excluded volumn parameter=',exc,sqrt(exc)*.87
      write(20,*)
      write(20,*)'................er1........',er1
      write(20,*)'................er2........',er2
      write(20,*)'................er3........',er3
      write(20,*)'................er4........',er4
      write(20,*)
      write(20,*)'initial random number seed...........',random
      write(20,*)'contact order........................',acorder
      write(20,*)
      write(20,*)'the first arandom number.............',aranzy(no)
      write(20,*) 
      write(20,*)'*******************************************'
      if(switch.eq.2)then
         write(20,*)'This is RS running of threading structure!'
      else
         write(20,*)'This is a normal simulation.'
      endif
      write(20,*)'*******************************************'
      write(20,*) 

c ^^^^^^^^^^ write parameter finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc set movement percentage cccccccccccccccccccc
cccccccccccccccccc CPU time of each movement ccccccccccccccccc
c                CPU        acceptance rate
c     move2      58s           11.5%
c     move3s     63s            7.5%
c     move3d     55s            6.0%
c     move4s     48s            3.3%
c     move4d     59s            2.6%
c     move4p     48s            3.3%
c     move5s     49s            1.8%
c     move5d     56s            1.3%
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_move_retio !set movement percentage
      implicit integer(i-z)
      common/commonuse2/atemp1,atemp2,N_rep,phot

      common/moveretio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc

      hsum=h2+h3s+h3d+h4s+h4d+h5s+h5d+h6+hend
      bh2=h2/hsum
      bh3s=bh2+h3s/hsum
      bh3d=bh3s+h3d/hsum
      bh4s=bh3d+h4s/hsum
      bh4d=bh4s+h4d/hsum
      bh5s=bh4d+h5s/hsum
      bh5d=bh5s+h5d/hsum
      bh6=bh5d+h6/hsum
      bhendn=bh6+hend/2./hsum
      bhendc=bhendn+hend/2./hsum

c      write(*,1)h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
c      write(*,1)bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6,
c     $     bhendn,bhendc
c      write(*,1)bh2,bh3s-bh2,bh3d-bh3s,bh4s-bh3d,bh4d-bh4s,bh5s-bh4d,
c     $     bh5d-bh5s,bh6-bh5d,bhendn-bh6,bhendc-bhendn
c 1    format(13f6.3)
c      stop

c ^^^^^^^^^^^^^^^^^^^^^^^ set retio finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set current (x,y,z), ica, T, from itemp's replica
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_current
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/echain2/egx(ndim),egy(ndim),egz(ndim)
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)
      common/excluded/vvv(ndim,ndim)
      common/looks/exc,exc1,exc2

*     put the itemp'th replica conformation as current conformation------>
*     moveable points:
      do i=1,Lch
         x(i)=xrep(i,itemp)     !initial coordinate of itemp's replica
         y(i)=yrep(i,itemp)
         z(i)=zrep(i,itemp)
      enddo

      do i=1,Lch1
         j=i+1
         wx=x(j)-x(i)
         wy=y(j)-y(i)
         wz=z(j)-z(i)
         ica(i)=vector(wx,wy,wz) !identify order number of each bond-vector
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

ccc   check what residues should be considered excluded volumn:
      do i=1,Lch
         do j=1,Lch
            vvv(i,j)=1          !every pair should be checked
            r2=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
            if(r2.lt.exc)then   !initially violated
               vvv(i,j)=-2
               if(abs(i-j).gt.2)then
c                  write(*,*)i,j,r2,vvv(i,j)
               endif
            endif
         enddo
      enddo
c^^^^^^^^^^^   vvv(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c^^^^^^^^^^^ set current (x,y,z,ica) finished ^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set current (x,y,z), ica, T, from itemp's replica
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_current_RS
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) !SG
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) !cc
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) !hb

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/chainm/mv(ndim)

      do i=1,Lch
         x(i)=xrep(i,itemp)
         y(i)=yrep(i,itemp)
         z(i)=zrep(i,itemp)
         ica(i)=icarep(i,itemp)
         if(mv(i).lt.0)then
            ex(i)=exrep(i,itemp) !Ca
            ey(i)=eyrep(i,itemp)
            ez(i)=ezrep(i,itemp)
            egx(i)=egxrep(i,itemp) !SG
            egy(i)=egyrep(i,itemp)
            egz(i)=egzrep(i,itemp)
            ecx(i)=ecxrep(i,itemp) !cc
            ecy(i)=ecyrep(i,itemp)
            ecz(i)=eczrep(i,itemp)
            ebx(i)=ebxrep(i,itemp) !hb
            eby(i)=ebyrep(i,itemp)
            ebz(i)=ebzrep(i,itemp)
         endif
      enddo
      ica(0)=ica(2)             !only useful for moveable pieces
      ica(Lch)=ica(Lch2)        !only useful for moveable pieces

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, translate (x,y,z) to calculate (amx,amy,amz)
c     2, calculate initial energy, 
c     3, prepare initial parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_move
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/seqe/seq(ndim),sec(ndim)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev

      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/chainm/mv(ndim)

      dimension axx(1000),ayy(1000),azz(1000)

*     move coordinates of the chain to its center of mass ---->
      call get_center
*     calculate (amx,amy,amz) -------------->
      call get_axis

cccccccc prepare all the initial parameters of this replica ccccccccc
      icnto=0                   !need when ISTAT>0
      sumcto=0                  !need when ISTAT>0
      do i=1,Lch
         nop(i)=0               !number of parallel contacts
         noa(i)=0
         nom(i)=0
      enddo

      energ=EHB(1,Lch,1)+ESHORT(1,Lch,10) !initialize all below parameters
ccc   eprofo,eprofn: calculated each time;
ccc   icnt: always start from icnto; (need initial)
ccc   sumct: always start from sumcto; (need initial)
ccc   nop(): start from 0, or nopp; (need initial)
ccc   noa(): start from 0, or noaa; (need initial)
ccc   nom(): start from 0, or nomm; (need initial)
ccc   codevsum: only used when istat=-1, i.e. in Enew.
ccc   didevsum: only used when istat=-1, i.e. in Enew.
ccc   afs(): keep in store.

cccc  backup all the parameters, calculated in EHB() and ESHORT() cccccccccc
      eprofo=0.0
      do k=1,Lch
         is=seq(k)
         ia=noa(k)              !number of antiparallel contact-apirs
         ip=nop(k)              !number of parallel contact-apirs ^^
         im=nom(k)              !number of orgonal contact-apirs
         eprofo=eprofo+envir(ia,im,ip,is,3)
      enddo

      icnto=icnt                !backup of contact-order
      sumcto=sumct              !backup of contact-order
      do i=1,Lch
         nopp(i)=nop(i)         !backup of number of contact-pair
         noaa(i)=noa(i)         !backup of number of contact-pair
         nomm(i)=nom(i)         !backup of number of contact-pair
      enddo
      codevsum=conew            !backup  panelity for total deviation of comb
      didevsum=dinew            !backup  panelity for total deviation of dist
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c^^^^^^^^^^^^^^^^^ initialization finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     move molecule to the orginal point.
c    we should calculate whole energy after each time when update the
c     center. Otherwise, the status of centersymmetric and bury of some
c     residue is unchanged with the update of center point. Thus
c    dE concerned with these off-dated residue is incorrect.
c     On the other hand, if one do not update center point for too long
c     time, the centrosymmetric potential may be not precise.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_center
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/initialinput/switch,k_cycle,k_phot,N_ann

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/chainm/mv(ndim)

      if(switch.eq.1)then       !!!!!!!!!normal run
         sx=0
         sy=0
         sz=0
         do i=1,Lch
            sx=sx+x(i)
            sy=sy+y(i)
            sz=sz+z(i)
         enddo
         sx=nint(float(sx)/float(Lch)) ! 4.4->4; 4.6->5
         sy=nint(float(sy)/float(Lch))
         sz=nint(float(sz)/float(Lch))
         do i=1,Lch
            x(i)=x(i)-sx
            y(i)=y(i)-sy
            z(i)=z(i)-sz
         enddo
      else                      !!!!!!!!template based run
***** calculate center of mass (sx,sy,sz)------->
         asx=0
         asy=0
         asz=0
         do i=1,Lch
            if(mv(i).gt.0)then
               asx=asx+x(i)
               asy=asy+y(i)
               asz=asz+z(i)
            else
               asx=asx+ex(i)
               asy=asy+ey(i)
               asz=asz+ez(i)
            endif
         enddo
         sx=nint(asx/float(Lch)) ! 4.4->4; 4.6->5
         sy=nint(asy/float(Lch))
         sz=nint(asz/float(Lch))
****   move coordinate to center of mass --------------->
         do i=1,Lch
            if(mv(i).gt.0)then
               x(i)=x(i)-sx
               y(i)=y(i)-sy
               z(i)=z(i)-sz
            else
               ex(i)=ex(i)-sx   !Ca
               ey(i)=ey(i)-sy
               ez(i)=ez(i)-sz
               egx(i)=egx(i)-sx !SG
               egy(i)=egy(i)-sy
               egz(i)=egz(i)-sz
            endif
         enddo
      endif

c^^^^^^^^^^^^^^^^^ centerization finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate (amx,amy,amz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_axis
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain2/egx(ndim),egy(ndim),egz(ndim)

      common/chainm/mv(ndim)
      common/initialinput/switch,k_cycle,k_phot,N_ann

      asx=0
      asy=0
      asz=0
      do i=1,Lch
         if(mv(i).gt.0)then
            asx=asx+(x(i)+gx(ica(i-1),ica(i),seq(i)))**2
            asy=asy+(y(i)+gy(ica(i-1),ica(i),seq(i)))**2
            asz=asz+(z(i)+gz(ica(i-1),ica(i),seq(i)))**2
         else
            asx=asx+egx(i)**2
            asy=asy+egy(i)**2
            asz=asz+egz(i)**2
         endif
      enddo
      amx=asx/Lch+0.00001
      amy=asy/Lch+0.00001
      amz=asz/Lch+0.00001
***   amx = sqrt(3)*x, x is the principle axis of the ellispoid.
***   amy = sqrt(3)*y, y is the principle axis of the ellispoid.
***   amz = sqrt(3)*z, z is the principle axis of the ellispoid.

c^^^^^^^^^^^^^^^^^ (amx,amy,amz) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c check whether passage of [jj,kk] overlap with other parts of chain.
c look =.ture., without overlap
c look =.false., with overlap
c only C_alpha is checked.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION LOOK(jj,kk)
      IMPLICIT INTEGER(I-Z)
      LOGICAL LOOK
      parameter(ndim=1000)
      parameter(nvec=312)
      common/seqe/seq(ndim),sec(ndim) 	 	
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2
      common/chainm/mv(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/excluded/vvv(ndim,ndim)

c     Ca-Ca r2>13;   r>sqrt(13)*0.87=3.14A
c     Ca-Ca r2>14;   r>sqrt(14)*0.87=3.26A
c     Ca-Ca r2>15;   r>sqrt(15)*0.87=3.37A
c     Ca-Ca r2>16;   r>sqrt(16)*0.87=3.48A
c     Ca-Ca r2>17;   r>sqrt(17)*0.87=3.59A
c     Ca-Ca r2>18;   r>sqrt(18)*0.87=3.69A
c     Ca-Ca r2>19;   r>sqrt(19)*0.87=3.79A
c     Ca-Ca r2>20;   r>sqrt(20)*0.87=3.89A
c     Ca-Ca r2>21;   r>sqrt(21)*0.87=3.99A
c     Ca-Ca r2>22;   r>sqrt(22)*0.87=4.08A
c     Ca-Ca r2>23;   r>sqrt(23)*0.87=4.17A *
c     Ca-Ca r2>24;   r>sqrt(24)*0.87=4.26A
c     Ca-Ca r2>25;   r>sqrt(25)*0.87=4.35A
c     Ca-Ca r2>26;   r>sqrt(26)*0.87=4.44A
c     Ca-Ca r2>27;   r>sqrt(27)*0.87=4.52A

      if(switch.eq.1)then       !normal lattice running
         do k=jj,kk
            i1=k-3
            i2=max(k+3,kk+1)
            do i=1,Lch
               if(i.le.i1.or.i.ge.i2)then
                  ir2=(x(k)-x(i))**2+(y(k)-y(i))**2+(z(k)-z(i))**2
                  if(ir2.lt.exc) then
                     LOOK=.FALSE.
                     RETURN
                  endif
               endif
            enddo
         enddo
      else                      !off-lattice
         do k=jj,kk
            if(mv(k).gt.0)then
               axk=x(k)
               ayk=y(k)
               azk=z(k)
            else
               axk=ex(k)
               ayk=ey(k)
               azk=ez(k)
            endif
            i1=k-3
            i2=max(k+3,kk+1)
            do i=1,Lch
               if(i.le.i1.or.i.ge.i2)then
                  if(vvv(i,k).gt.0)then
                     if(mv(i).gt.0)then
                        axi=x(i)
                        ayi=y(i)
                        azi=z(i)
                     else
                        axi=ex(i)
                        ayi=ey(i)
                        azi=ez(i)
                     endif
                     ar2=(axk-axi)**2+(ayk-ayi)**2+(azk-azi)**2
                     if(ar2.lt.exc) then
                        LOOK=.FALSE.
                        RETURN
                     endif
                  endif
               endif
            enddo
         enddo
      endif
      LOOK=.TRUE.

c ^^^^^^^^^^ Look finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, pairwise interaction;
c     2, H-bond;
c     3, Electric;
c     4, contact number;
c     5, contact order.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION EHB(jjjj,kkkk,ISTAT)

      IMPLICIT INTEGER(I-Z)
      parameter(ndim=1000)
                parameter(nvec=312)
      common/seqe/seq(ndim),sec(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/icgg/ icg(ndim),EH6  		
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3 
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)	
      common/pair1/eh2,eh1b
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB2,EHB3,EHB4,EHB5,EHB6
      common/distres/er4,er2,es3c
      common/chainm/mv(ndim)
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      common/ehbenergy1/EHB5a,EHB5b

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb


      EHB1=0                    !+1/r of Ca-SC
      EHB1a=0                   !+1/r for non-parallel contact of Ca-Ca
      EHB1b=0                   !excluded volume of SC-SC

      EHB2=0                    !pairwise for SC-SC
      EHB3=0                    !enhance good-piece contacts
      EHB4=0                    !-1/r for parallel contact of Ca-Ca

      EHB5a=0                    !H-bond energy for alpha-type
      EHB5b=0                    !H-bond energy for beta-type

      if(ISTAT.gt.0)THEN        !Enew
         ICNT=ICNTO
         SUMCT=SUMCTO
      ENDIF

c     coupling of secondary structure with pairwise interactions
c     included - thus expanded range
      jj=jjjj-1
      if(jj.lt.1)jj=1
      kk=kkkk+1
      if(kk.gt.Lch)kk=Lch
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
********** circle for movement involved window ******************
*****************************************************************
      DO 1002 k=jj,kk
         kseq=seq(k)
         if(mv(k).gt.0)then     !moveable point
            axk=x(k)            !Ca
            ayk=y(k)
            azk=z(k)
            nv1=ica(k-1)
            nv2=ica(k)
            agxk=axk+GX(nv1,nv2,kseq) !Cg_k
            agyk=ayk+GY(nv1,nv2,kseq)
            agzk=azk+GZ(nv1,nv2,kseq)
            bxk=HBX(nv1,nv2)    !H-bond direction
            byk=HBY(nv1,nv2)
            bzk=HBZ(nv1,nv2)
            cxk=CAX(nv1,nv2)    !vertical vector
            cyk=CAY(nv1,nv2)
            czk=CAZ(nv1,nv2)	
         else
            axk=ex(k)           !Ca
            ayk=ey(k)
            azk=ez(k)
            agxk=egx(k)         !SG
            agyk=egy(k)
            agzk=egz(k)
            bxk=ebx(k)          !Hb
            byk=eby(k)
            bzk=ebz(k)
            cxk=ecx(k)          !cc
            cyk=ecy(k)
            czk=ecz(k)
         endif
         km2=k-2
         kp2=k+2
         if(km2.ge.1.and.kp2.le.Lch)then
            xxx=nint((aax(km2)-aax(kp2))**2+
     $           (aay(km2)-aay(kp2))**2+(aaz(km2)-aaz(kp2))**2) !real dist
            if(xxx.gt.500)xxx=500
            ek5=acops(km2,jbin(xxx)) !k2-residue, jbin. es<0.
         else
            ek5=0
         endif
*****************************************************************
********************** start i-cicle ****************************
*****************************************************************
         do 1001 i=1,Lch
            iend=max(k+1,kk)
            if(i.ge.k-1.and.i.le.iend)goto 1001 !to avoid repeat
            if(mv(i).gt.0)then
               axi=x(i)
               ayi=y(i)
               azi=z(i)
            else
               axi=ex(i)
               ayi=ey(i)
               azi=ez(i)
            endif
            axki=(axk-axi)      !r_k-r_i
            ayki=(ayk-ayi)
            azki=(azk-azi)
            Cr2=axki**2+ayki**2+azki**2+0.000001 !(r_k-r_i)**2
            if(Cr2.lt.120)then  !120-->9.53, there is pairwise interaction
               idist=iabs(i-k)
               iseq=seq(i)
               if(mv(i).gt.0)then
                  nv1=ica(i-1)
                  nv2=ica(i)
                  agxi=axi+GX(nv1,nv2,iseq) !Cg_k
                  agyi=ayi+GY(nv1,nv2,iseq)
                  agzi=azi+GZ(nv1,nv2,iseq)
                  bxi=HBX(nv1,nv2) !H-bond direction
                  byi=HBY(nv1,nv2)
                  bzi=HBZ(nv1,nv2)
                  cxi=CAX(nv1,nv2) !outside 
                  cyi=CAY(nv1,nv2)
                  czi=CAZ(nv1,nv2)	
               else
                  agxi=egx(i)   !SG
                  agyi=egy(i)
                  agzi=egz(i)
                  bxi=ebx(i)    !H-bond
                  byi=eby(i)
                  bzi=ebz(i)
                  cxi=ecx(i)    !bisector vector
                  cyi=ecy(i)
                  czi=ecz(i)
               endif

c     1/r excluded for Ca-SC pair ----------------->
               if(kseq.gt.0) then !not GLY, we have SG
                  aks=(agxk-axi)**2+(agyk-ayi)**2+(agzk-azi)**2+.0001 !SG_kCa_i
                  if(aks.lt.36) then !36-->5.22A
                     if(aks.lt.13.0) aks=13.0 !13-->3.17A
                     EHB1=EHB1+((13.0/aks)-0.5)
                  endif
               endif
               if(iseq.gt.0) then
                  aks=(axk-agxi)**2+(ayk-agyi)**2+(azk-agzi)**2+.0001 !Ca_kSG_i
                  if(aks.lt.36.0) then
                     if(aks.lt.13.0) aks=13.0
                     EHB1=EHB1+((13.0/aks)-0.5)
                  endif
               endif

c     quarsi3 for SC-SC pair, 1/r for Ca-Ca ----------------->
               Gr2=(agxi-agxk)**2+(agyi-agyk)**2+(agzi-agzk)**2 !SG_i,SG_k
               cc=cxi*cxk+cyi*cyk+czi*czk !c_i*c_k
c               write(*,*)i+1,k+1,cc,Gr2
               IF(cc.gt.0.5)THEN !c_i//c_k
                  if(Cr2.lt.60)EHB4=EHB4-(30/(max(30.,Cr2))-0.5) !6.74A
                  if(Gr2.lt.alp(iseq,kseq))then
                     EHB3=EHB3-ek5*ei5(i,idist)
                     NOP(k)=NOP(k)+ISTAT
                     NOP(i)=NOP(i)+ISTAT
                     ICNT=ICNT+istat*idist
                     SUMCT=SUMCT+ISTAT
                     if(Gr2.gt.arlp(iseq,kseq))THEN !excluded volume
                        EHB2=EHB2+app(i,k)
                     else
                        EHB1b=EHB1b+1
                     endif
                  endif
               ELSE
                  if(Cr2.lt.33)EHB1a=EHB1a+(16.0/Cr2-0.5) !5A
                  IF(cc.lt.-0.5) THEN !antiparallel pair-wise of (i,k)
                     if(Gr2.lt.ala(iseq,kseq))then
                        EHB3=EHB3-ek5*ei5(i,idist)
                        NOA(k)=NOA(k)+ISTAT
                        NOA(i)=NOA(i)+ISTAT
                        ICNT=ICNT+istat*idist
                        SUMCT=SUMCT+ISTAT
                        if(Gr2.gt.arla(iseq,kseq))THEN
                           EHB2=EHB2+apa(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
                  ELSE          !neither parallel nor antiparallel
                     if(Gr2.lt.alm(iseq,kseq)) then
                        EHB3=EHB3-ek5*ei5(i,idist)
                        NOM(k)=NOM(k)+ISTAT
                        NOM(i)=NOM(i)+ISTAT
                        ICNT=ICNT+istat*idist !distance of pairs
                        SUMCT=SUMCT+ISTAT !number of contact pairs
                        if(Gr2.gt.arlm(iseq,kseq))THEN
                           EHB2=EHB2+apm(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
                  ENDIF
               ENDIF
c^^^^^^^^^^^^^ interaction of pair-wise finished ^^^^^^^^^^^^^^^^^

ccccccccccccccccccccccccc Hydrogen-bond cccccccccccccccccccccccccc
***** 1, cc; 2, bb; 3, distance; 4, vij; 5, sec(i).
c               if(idist.eq.3)then
c                  alpha-H-bond
c               elseif(idist.lt.20)then
c                  anti-parallel-sheet-H-bond
c               else
c                  if(bb<0)then
c                     antiparallel sheet-H-bond
c                  else
c                     parallel sheet-H-bond
c                  endif
c               endif
*HHHHH***************** alpha ******************************************
***   bb   bb_an  cc    cc_an vv   vv_an hbl r2i  r2k  Cr2  Cr2m  Hdis
***   0.82  34.9  0.42  65.3  0.52  58.8 5.7 3.29 3.53 35.9 31.62 4.19 (alpha)
***   0.81  35.3  0.38  67.3  0.48  61.0 5.7 3.55 3.39 35.5 30.45 4.10 (alpha1)
***   0.85  31.9  0.50  59.5  0.58  54.0 5.7 3.03 2.92 34.1 30.29 4.13 (alpha2)
************************************************************************
*********************** beta *******************************************
***   bb    bb_an   cc  cc_an vv12 vv12_a vv12 vv12_a hbl r2i r2k Cr2 Cr2m
***   -0.81 144.8 0.79  36.7 -0.94 160.6 -0.95 162.4 5.3 8.38 1.37 30 30 (r)
***   -0.84 149.1 0.83  31.8 -0.92 161.1 -0.94 162.7 5.3 3.91 4.43 33 31 (r1)
***   -0.64 135.0 0.46  61.8 -0.75 139.6 -0.74 137.6 5.1 7.81 6.55 33 31 (r2)
***    0.98   8.6 0.92  22.2  0.92  22.5  0.93  18.8 5.3 2.49 0.99 31 30 (p)
************************************************************************
               if(i.eq.1.or.i.eq.Lch.or.k.eq.1.or.k.eq.Lch)goto 1003
               if(idist.eq.3)then !!!!!!!!!!!!!!!!!!!!!!!!!!!
***   alpha-helix:
                  if(sec(k).ne.4.and.sec(i).ne.4)then
                  if(Cr2.lt.Cr2a.and.cc.gt.acut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
                  if(bb.gt.acut_bb)then
                  av11=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                  if(av11.gt.acut_vv)then
                  av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                  if(av22.gt.acut_vv)then
***   ->
                     fact=(1-abs(cc-0.4))*(1-abs(bb-0.815))
                     EHB5a=EHB5a+energyHBa(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                  endif
                  endif
                  endif
                  endif
                  endif
               elseif(idist.gt.4.and.idist.lt.20)then !!!!!!!!!!!!!!!!!!!!!!!
***   antiparallel-sheet
                  if(sec(k).ne.2.and.sec(i).ne.2)then
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
                  if(bb.lt.-bcut_bb)then !antiparallel
                  av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                  if(av12.lt.-bcut_vv)then
                  av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                  if(av21.lt.-bcut_vv)then
***   ->
                     fact=abs(bb)*cc !bb->-1,cc->1
                     EHB5b=EHB5b+energyHBb(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                  endif
                  endif
                  endif
                  endif
                  endif
               elseif(idist.ge.20)then !!!!!!!!!!!!!!!!!!!!!!!!
                  if(sec(k).ne.2.and.sec(i).ne.2)then
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk) 
                  if(bb.lt.-bcut_bb)then !antiparallel
***   antiparallel-sheet:
                     av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                     if(av12.lt.-bcut_vv)then
                     av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                       aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                     if(av21.lt.-bcut_vv)then
***   ->
                        fact=abs(bb)*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                     endif
                     endif
                  elseif(bb.gt.bcut_bb)then !bb>0, parallel H-bond
***   parallel-sheet:
                     av11=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                     if(av11.gt.bcut_vv)then
                     av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                       axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                     if(av22.gt.bcut_vv)then
***   ->
                        fact=bb*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                     endif
                     endif
                  endif
                  endif
                  endif
               endif            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 1003          continue
c^^^^^^^^^^^^^H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            endif               !Ca-Ca<120
 1001    continue               !i -> [1,Lch]
 1002 continue                  !k -> [jj,kk]

c     ICNT/SUMCT: average distance of residues in each contact-pairs.
c     ICNT/SUMCT are calculated when call Enew
c     ICNTO/SUMCTO are not changed when call Eold
      if(istat.lt.0) then       !Eold
         b=abs(ICNT/(0.000001+float(SUMCT))-acorder)
         a=abs(ICNTO/(0.000001+float(SUMCTO))-acorder)
         d=abs(float(SUMCT)-contt) !deviation of contact number on new conform
         c=abs(float(SUMCTO)-contt) !deviation of contact number on new confor
         dord=en2*(b-a)+en3*(d-c) !not included in EHB
      endif

      EHB
     $     =eh1*EHB1            !+1/r of Ca-SC
     $     +eh1a*EHB1a          !+1/r for non-parallel of Ca-Ca
     $     +eh1b*EHB1b          !excluded volumn of SC-SC
     $     +eh2*EHB2            !pairwise for SC-SC
     $     +eh3*EHB3            !enhance good piece
     $     +eh4*EHB4            !-1/r for parallel contact of Ca-Ca
     $     +eh5a*EHB5a          !H-bond energy (alpha)
     $     +eh5b*EHB5b          !H-bond energy (beta)

c ^^^^^^^^^^ EHB finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     H-bond for alpha-helix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energyHBa(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk)
      implicit integer(i-z)
      parameter(ndim=1000)      !maximum length of chain-length
      parameter(nvec=312)       !number of vectors
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh

c     because helix is right-hand, bi=v(i-1)(x)v(i) is always
c     same direction as v(i---k), if (k>i);
c     reverse direction as v(i---k), if (k<i);

      energyHBa=0

      bxk=5.7*dxk
      byk=5.7*dyk
      bzk=5.7*dzk
      bxi=5.7*dxi
      byi=5.7*dyi
      bzi=5.7*dzi
      if(k.gt.i)then            !bki,bxk,axki are in same directory
         br=(bxk-axki)**2+(byk-ayki)**2+(bzk-azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2)) !br->3.4
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
         br=(bxi-axki)**2+(byi-ayki)**2+(bzi-azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
      else                      !bki,bxk, and axki are in reverse directory
         br=(bxk+axki)**2+(byk+ayki)**2+(bzk+azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
         br=(bxi+axki)**2+(byi+ayki)**2+(bzi+azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
      endif

c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     H-bond for beta-sheet
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energyHBb(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk)
      implicit integer(i-z)
      parameter(ndim=1000)      !maximum length of chain-length
      parameter(nvec=312)       !number of vectors
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh

      energyHBb=0

      bxk=5.3*dxk
      byk=5.3*dyk
      bzk=5.3*dzk
      bxi=5.3*dxi
      byi=5.3*dyi
      bzi=5.3*dzi
      br=(bxk-axki)**2+(byk-ayki)**2+(bzk-azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)           !br->0
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxk+axki)**2+(byk+ayki)**2+(bzk+azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxi-axki)**2+(byi-ayki)**2+(bzi-azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxi+axki)**2+(byi+ayki)**2+(bzi+azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif

c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate v(12).v(34)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function avv(ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3,ax4,ay4,az4)

      ax12=ax2-ax1
      ay12=ay2-ay1
      az12=az2-az1
      ax34=ax4-ax3
      ay34=ay4-ay3
      az34=az4-az3
      avv=(ax12*ax34+ay12*ay34+az12*az34)/(0.000001+
     $     sqrt((ax12**2+ay12**2+az12**2)*(ax34**2+ay34**2+az34**2)))

c^^^^^^^^^^^^^ v.v finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate ei5=acops(i,jbin(r15))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ei5(i,idist)
      implicit integer(i-z)
      parameter(ndim=1000)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/lengths/Lch,Lch1,Lch2
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)

      ei5=0
      if(idist.gt.5)then
         im2=i-2
         ip2=i+2
         if(im2.ge.1.and.ip2.le.Lch)then
            if(mv(im2).gt.0)then
               axm=x(im2)
               aym=y(im2)
               azm=z(im2)
            else
               axm=ex(im2)
               aym=ey(im2)
               azm=ez(im2)
            endif
            if(mv(ip2).gt.0)then
               axp=x(ip2)
               ayp=y(ip2)
               azp=z(ip2)
            else
               axp=ex(ip2)
               ayp=ey(ip2)
               azp=ez(ip2)
            endif
            xxx=nint(axm-axp)**2+(aym-ayp)**2+(azm-azp)**2
            if(xxx.gt.500)xxx=500
            ei5=acops(im2,jbin(xxx)) !i2-residue, jbin. es<0.
         endif
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, centrosymmetric energy of (x,y,z) and Cg. (centro.comm).
c     2, bury potential of SC.
c     3, distmap and contact restrains.
c     4, bias to (prodicted) protein-like structure, panality on crumpling.
c     5, E13,E15,E15.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION ESHORT(iiii,jjjj,ISTAT)
      IMPLICIT INTEGER(I-Z)
      parameter(ndim=1000)	
                parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim) 	  
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/fr/frga(ndim),frgb(ndim)
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/three/angle(nvec,nvec)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      COMMON/RES/ ER3, Mcom(ndim),Kcom(ndim,50)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/msichores/msicho
      common/distres/er4,er2,es3c
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT3a,ESHORT4
      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT9a
      common/eshortenergy4/ESHORT9b,ESHORT9c
      common/hopp/eonehw(0:19)
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      common/chainm/mv(ndim)

      common/hom1/asr2(ndim,4),asr3(ndim,4),asr4(ndim,14),asr5(ndim,8)
      common/hom2/ibb2(0:999),ibb3(0:999),ibb4(-999:999),ibb5(0:999)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/zscore/azscore

cc      ESHORT=0.0
c      ESHORT1=0                 !centrosymmetric potential for Ca
      ESHORT2=0                 !bury potential for SC

      ESHORT3=0  !distance restrain from both threading (for C_a)
      ESHORT3a=0 !deviation of distance restrain
      ESHORT4=0  !contact restrain from threading (for SC)
      ESHORT4a=0 !deviation of contact restrain

      ESHORT5=0  !bias2,3: v(i)-v(i+4) anti/parallel; c(i)-c(i+2) anit/paralel
      ESHORT5a=0 !crumpling
      ESHORT5b=0 !bias4 to predicted alpha/beta structure.
      ESHORT5c=0 !bias1 to possible alpha/beta structure. 

      ESHORT6=0  !correlation of E13 of Ca
      ESHORT7=0  !correlation of E14, from both common and 2th specific data
      ESHORT8=0  !correlation of E15, from both common and 2th specific data

      if(istat.lt.0)THEN       !Eold 
         diold=0.0
         coold=0.0
      else                      !Enew
         conew=0.0   !penality for incorrect contact (from contact restrain)
         dinew=0.0   !penality for incorrect distance (from distant restain)
      endif

      do 1 i=iiii,jjjj
         iseq=seq(i)
         if(mv(i).gt.0)then
            axi=x(i)
            ayi=y(i)
            azi=z(i)
            nv1=ica(i-1)
            nv2=ica(i)
            agxi=axi+GX(nv1,nv2,iseq) !Cg_k
            agyi=ayi+GY(nv1,nv2,iseq)
            agzi=azi+GZ(nv1,nv2,iseq)
         else
            axi=ex(i)           !Ca
            ayi=ey(i)
            azi=ez(i)
            agxi=egx(i)         !SG
            agyi=egy(i)
            agzi=egz(i)
         endif
***** Centrosymmetric energy of Ca ---------->
         aract2=axi**2+ayi**2+azi**2+0.01
         aract=sqrt(aract2)     !radius: r
         ff=aract/acrit+0.000001 !ff=r/r0
c         ika=int(ff*3.0)        !=0, r<r0/3; 1, r0/3<r<2r0<3.
c         ESHORT1=ESHORT1+eoinp(iseq,ika) !k: ID, ika: distance from center
***** Bury energy of SG --------------------->
         fff=agxi**2/amx+agyi**2/amy+agzi**2/amz-3.0
         if(fff.lt.0.0) then    !C_g inside ellipsoid
            if(fff.lt.-1.0) fff=-1.0 !C_g inside core of ellipsoid
            eka=eonekd(iseq)
            eha=eonehw(iseq)
            ESHORT2=ESHORT2+fff*eka !Kyte-Doolittle, + hydrophobic
            ESHORT2=ESHORT2-fff*eha !Hopp-Woods, + hydrophilic
         endif
c^^^^^^^^^^^^^^^ centrosymmetric energy finished ^^^^^^^^^^^^^^^^^^

         ff=min(1.0, 1.0/ff)
         afs(i)=max(ff,0.5)     !=1, r<r0; [0.5,1], r0<r<2r0; 0.5, r>2r0
c     afs(i) burial factors used for modulation of chain stifness
c     rigid inside, more flexible on the surface
c     afs(i) in [0.5,1.0]. 0.5, outside circle; 1.0 inside circle.
c         if(r<r0)then
c            afs(i)=1
c         elseif(r<2*r0)then
c            afs(i)=r0/r
c         else
c            afs=0.5
c         endif

****  restrain from 'dist.dat'---------------------->
         N_dis=Mdis(i)          !number of distance restraints on 'i'
         do k=1,N_dis
            j=kdis(i,k)         !i-j restraints
            if(j.lt.i.or.j.gt.jjjj) then !to avoid repeat
               dij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
               dij=sqrt(dij2)
               err=abs(dij-dist(i,k)) !arca: predicted dis
               if(err.gt.dev(i,k)) then !dev: deviation for arca
                  ESHORT3=ESHORT3+1
                  if(istat.lt.0) then !Eold
                     diold=diold+err/dev(i,k) !accumlation of err
                  else
                     dinew=dinew+err/dev(i,k)
                  endif
               endif
            endif
         enddo
c^^^^^^^^^^^^^^^^^ distant restrain finished ^^^^^^^^^^^^^^^^^^^^^^^^^
1      continue

ccc   Contact restrain from 'comb.dat' -------------------------------->
       if(azscore.gt.10)then    !exact restrains
          do 11 i=iiii,jjjj
             N_com=Mcom(i)      !number of contacts on i
             if(N_com.ge.1)then
                iseq=seq(i)
                if(mv(i).gt.0)then
                   agxi=x(i)+GX(ica(i-1),ica(i),iseq) !SG
                   agyi=y(i)+GY(ica(i-1),ica(i),iseq) !SG
                   agzi=z(i)+GZ(ica(i-1),ica(i),iseq) !SG
                else
                   agxi=egx(i)  !SG
                   agyi=egy(i)
                   agzi=egz(i)
                endif
                
                do k=1,N_com
                   j=Kcom(i,k)  !k'th contact with i
                   if(j.lt.i.OR.j.gt.jjjj)then
                      jseq=seq(j)
                      if(mv(j).gt.0)then
                         agxj=x(j)+gx(ica(j-1),ica(j),jseq)
                         agyj=y(j)+gy(ica(j-1),ica(j),jseq)
                         agzj=z(j)+gz(ica(j-1),ica(j),jseq)
                      else
                         agxj=egx(j)
                         agyj=egy(j)
                         agzj=egz(j)
                      endif
                      aij2=(agxi-agxj)**2+(agyi-agyj)**2+(agzi-agzj)**2
                 
                      if(aij2.gt.concut2(iseq,jseq))then
                         ESHORT4=ESHORT4+2 !no contact
                         if(istat.lt.0)then
                            coold=coold+sqrt(aij2)-concut(iseq,jseq) !penality
                         else
                            conew=conew+sqrt(aij2)-concut(iseq,jseq) !panelity
                         endif
                      endif
                   endif
                enddo
             endif
 11       continue
       else                     !!!!!!! +-1 fuzzy restrains
        i1=iiii-1
        if(i1.lt.1)i1=1
        i2=jjjj+1
        if(i2.gt.Lch)i2=Lch
        do 22 i=i1,i2
        N_com=Mcom(i)      !number of contacts on i
        if(N_com.ge.1)then
cccccc SG_i:
           if(mv(i).gt.0)then
              agxi=x(i)+GX(ica(i-1),ica(i),seq(i)) !SG_i
              agyi=y(i)+GY(ica(i-1),ica(i),seq(i))
              agzi=z(i)+GZ(ica(i-1),ica(i),seq(i))
           else
              agxi=egx(i)
              agyi=egy(i)
              agzi=egz(i)
           endif
ccccccSG_i-1:
           if(i.gt.1)then
              if(mv(i-1).gt.0)then
                 agxim=x(i-1)+GX(ica(i-2),ica(i-1),seq(i-1))
                 agyim=y(i-1)+GY(ica(i-2),ica(i-1),seq(i-1))
                 agzim=z(i-1)+GZ(ica(i-2),ica(i-1),seq(i-1))
              else
                 agxim=egx(i-1) !SG
                 agyim=egy(i-1)
                 agzim=egz(i-1)
              endif
           else
              agxim=agxi
              agyim=agyi
              agzim=agzi
           endif
ccccccSG_i+1:
           if(i.lt.Lch)then
              if(mv(i+1).gt.0)then
                 agxip=x(i+1)+GX(ica(i),ica(i+1),seq(i+1))
                 agyip=y(i+1)+GY(ica(i),ica(i+1),seq(i+1))
                 agzip=z(i+1)+GZ(ica(i),ica(i+1),seq(i+1))
              else
                 agxip=egx(i+1) !SG
                 agyip=egy(i+1)
                 agzip=egz(i+1)
              endif
           else
              agxip=agxi
              agyip=agyi
              agzip=agzi
           endif
             
           do k=1,N_com
              j=Kcom(i,k)       !k'th contact with i
              if(j.lt.i.OR.j.gt.jjjj+1)then
ccccccSG_j:
                 if(mv(j).gt.0)then
                    agxj=x(j)+gx(ica(j-1),ica(j),seq(j))
                    agyj=y(j)+gy(ica(j-1),ica(j),seq(j))
                    agzj=z(j)+gz(ica(j-1),ica(j),seq(j))
                 else
                    agxj=egx(j)
                    agyj=egy(j)
                    agzj=egz(j)
                 endif
cccccc SG_j-1:
                 if(j.gt.1)then
                    if(mv(j-1).gt.0)then
                       agxjm=x(j-1)+gx(ica(j-2),ica(j-1),seq(j-1))
                       agyjm=y(j-1)+gy(ica(j-2),ica(j-1),seq(j-1))
                       agzjm=z(j-1)+gz(ica(j-2),ica(j-1),seq(j-1))
                    else
                       agxjm=egx(j-1) !SG
                       agyjm=egy(j-1)
                       agzjm=egz(j-1)
                    endif
                 else
                    agxjm=agxj
                    agyjm=agyj
                    agzjm=agzj
                 endif
cccccc SG_j+1:
                 if(j.lt.Lch)then
                    if(mv(j+1).gt.0)then
                       agxjp=x(j+1)+gx(ica(i),ica(i+1),seq(j+1))
                       agyjp=y(j+1)+gy(ica(i),ica(i+1),seq(j+1))
                       agzjp=z(j+1)+gz(ica(i),ica(i+1),seq(j+1))
                    else
                       agxjp=egx(j+1)
                       agyjp=egy(j+1)
                       agzjp=egz(j+1)
                    endif
                 else
                    agxjp=agxj
                    agyjp=agyj
                    agzjp=agzj
                 endif

                 ar1=(agxi-agxj)**2+(agyi-agyj)**2+(agzi-agzj)**2
                 ar2=(agxi-agxjm)**2+(agyi-agyjm)**2+(agzi-agzjm)**2
                 ar3=(agxi-agxjp)**2+(agyi-agyjp)**2+(agzi-agzjp)**2

                 ar4=(agxim-agxj)**2+(agyim-agyj)**2+(agzim-agzj)**2
c                 ar5=(agxim-agxjm)**2+(agyim-agyjm)**2+(agzim-agzjm)**2
c                 ar6=(agxim-agxjp)**2+(agyim-agyjp)**2+(agzim-agzjp)**2 !

                 ar7=(agxip-agxj)**2+(agyip-agyj)**2+(agzip-agzj)**2
c                 ar8=(agxip-agxjm)**2+(agyip-agyjm)**2+(agzip-agzjm)**2 !
c                 ar9=(agxip-agxjp)**2+(agyip-agyjp)**2+(agzip-agzjp)**2
c                 ar=min(ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9)
c                 ar=min(ar1,ar2,ar3,ar4,ar5,ar7,ar9) !w/o shift 2
                 ar=min(ar1,ar2,ar3,ar4,ar7) !w/o double shift

                 if(ar1.lt.concut2(seq(i),seq(j)))then
                    ESHORT4=ESHORT4-1
                 else           !panalty
                    if(istat.lt.0)then
                       coold=coold+sqrt(ar1)-concut(seq(i),seq(j)) !penality
                    else
                       conew=conew+sqrt(ar1)-concut(seq(i),seq(j)) !panelity
                    endif
                 endif
                 if(ar.lt.concut2(seq(i),seq(j)))ESHORT4=ESHORT4-1
              endif
            enddo
          endif
 22     continue
      endif
c^^^^^^^^^^^^^^^^^^ contact restrains finished ^^^^^^^^^^^^^^^^^^^^^

*********E13 --------------------------------------->
      i1=max(iiii-1,1)
      i2=min(jjjj-1,Lch-2)
      do i=i1,i2
         ar13=(aax(i)-aax(i+2))**2+(aay(i)-aay(i+2))**2+
     $        (aaz(i)-aaz(i+2))**2
         if(ar13.lt.48) then    !6.03A
            ESHORT6=ESHORT6+csr(i,1)
         else
            ESHORT6=ESHORT6+csr(i,2)
         endif
      enddo
c^^^^^^^^^^^^^^ E_13 finished ^^^^^^^^^^^^^^^^^^^^^^^^^

*********E14, E15, proteinlike bias1 to possible alpha/beta --------------->
      i1=max(iiii-4,1)
      i2=min(jjjj,Lch-4)
      do 2 i=i1,i2
         if(mv(i).gt.0)then
            ax1=x(i)
            ay1=y(i)
            az1=z(i)
            agx1=x(i)+gx(ica(i-1),ica(i),seq(i))
            agy1=y(i)+gy(ica(i-1),ica(i),seq(i))
            agz1=z(i)+gz(ica(i-1),ica(i),seq(i))
         else
            ax1=ex(i)
            ay1=ey(i)
            az1=ez(i)
            agx1=egx(i)
            agy1=egy(i)
            agz1=egz(i)
         endif
         if(mv(i+1).gt.0)then
            ax2=x(i+1)
            ay2=y(i+1)
            az2=z(i+1)
            agx2=x(i+1)+gx(ica(i),ica(i+1),seq(i+1))
            agy2=y(i+1)+gy(ica(i),ica(i+1),seq(i+1))
            agz2=z(i+1)+gz(ica(i),ica(i+1),seq(i+1))
         else
            ax2=ex(i+1)
            ay2=ey(i+1)
            az2=ez(i+1)
            agx2=egx(i+1)
            agy2=egy(i+1)
            agz2=egz(i+1)
         endif
         if(mv(i+2).gt.0)then
            ax3=x(i+2)
            ay3=y(i+2)
            az3=z(i+2)
            agx3=x(i+2)+gx(ica(i+1),ica(i+2),seq(i+2))
            agy3=y(i+2)+gy(ica(i+1),ica(i+2),seq(i+2))
            agz3=z(i+2)+gz(ica(i+1),ica(i+2),seq(i+2))
         else
            ax3=ex(i+2)
            ay3=ey(i+2)
            az3=ez(i+2)
            agx3=egx(i+2)
            agy3=egy(i+2)
            agz3=egz(i+2)
         endif
         if(mv(i+3).gt.0)then
            ax4=x(i+3)
            ay4=y(i+3)
            az4=z(i+3)
            agx4=x(i+3)+gx(ica(i+2),ica(i+3),seq(i+3))
            agy4=y(i+3)+gy(ica(i+2),ica(i+3),seq(i+3))
            agz4=z(i+3)+gz(ica(i+2),ica(i+3),seq(i+3))
         else
            ax4=ex(i+3)
            ay4=ey(i+3)
            az4=ez(i+3)
            agx4=egx(i+3)
            agy4=egy(i+3)
            agz4=egz(i+3)
         endif
         if(mv(i+4).gt.0)then
            ax5=x(i+4)
            ay5=y(i+4)
            az5=z(i+4)
            agx5=x(i+4)+gx(ica(i+3),ica(i+4),seq(i+4))
            agy5=y(i+4)+gy(ica(i+3),ica(i+4),seq(i+4))
            agz5=z(i+4)+gz(ica(i+3),ica(i+4),seq(i+4))
         else
            ax5=ex(i+4)
            ay5=ey(i+4)
            az5=ez(i+4)
            agx5=egx(i+4)
            agy5=egy(i+4)
            agz5=egz(i+4)
         endif
ccccccE14:
         ax=ax2-ax1
         ay=ay2-ay1
         az=az2-az1
         bx=ax3-ax2
         by=ay3-ay2
         bz=az3-az2
         cx=ax4-ax3
         cy=ay4-ay3
         cz=az4-az3
         abx=ay*bz-az*by
         aby=az*bx-ax*bz
         abz=ax*by-ay*bx
         hand=abx*cx+aby*cy+abz*cz !a(x)b.c, chirality, >0, right-hand
         ar14=(ax1-ax4)**2+(ay1-ay4)**2+(az1-az4)**2
         r14=nint(ar14)
         if(r14.gt.300) r14=300
         if(hand.lt.0) r14=-r14 !<0, left-hand three-bond.
         ESHORT7=ESHORT7+asr(i,ibin(r14)) !asr(i,dis(i,i+4)) from r14*.dat
cccccccccE15:
         ar15=(ax1-ax5)**2+(ay1-ay5)**2+(az1-az5)**2
         r15=nint(ar15)
         if(r15.gt.500) r15=500
         ESHORT8=ESHORT8+bsr(i,jbin(r15))

cccccccccbias1: encourage helix/sheet-like structure
         if(ar15.lt.75)THEN     !7.53A: helix
            if(sec(i+1).ne.4.AND.sec(i+2).ne.4.AND.sec(i+3).ne.4)then
               if(ibin(r14).gt.4.AND.ibin(r14).lt.8)then !right-hand,3A->8A
                  dot13=(ax4-ax3)*(ax2-ax1)
     $                 +(ay4-ay3)*(ay2-ay1)+(az4-az3)*(az2-az1)
                  if(dot13.lt.0)then
                     dot24=(ax5-ax4)*(ax3-ax2)
     $                    +(ay5-ay4)*(ay3-ay2)+(az5-az4)*(az3-az2)
                     if(dot24.lt.0)then
                        dot14=(ax5-ax4)*(ax2-ax1)
     $                       +(ay5-ay4)*(ay2-ay1)+(az5-az4)*(az2-az1)
                        if(dot14.gt.0)then
                           ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
                           ESHORT5c=ESHORT5c-2-ff
                        endif
                     endif
                  endif
               endif
            endif
         elseif(ar15.gt.160)then !11A, beta
            if(sec(i+1).ne.2.AND.sec(i+2).ne.2.AND.sec(i+3).ne.2)then
               if(mv(i+1).gt.0)then
                  bx2=hbx(ica(i),ica(i+1))
                  by2=hby(ica(i),ica(i+1))
                  bz2=hbz(ica(i),ica(i+1))
               else
                  bx2=ebx(i+1)  !Hb
                  by2=eby(i+1)
                  bz2=ebz(i+1)
               endif
               if(mv(i+3).gt.0)then
                  bx4=hbx(ica(i+2),ica(i+3))
                  by4=hby(ica(i+2),ica(i+3))
                  bz4=hbz(ica(i+2),ica(i+3))
               else
                  bx4=ebx(i+3)
                  by4=eby(i+3)
                  bz4=ebz(i+3)
               endif
               dot24=bx2*bx4+by2*by4+bz2*bz4 !cos(H1.H3)
               if(dot24.gt.0.71)then !angle<45
                  if(mv(i+2).gt.0)then
                     bx3=hbx(ica(i+1),ica(i+2))
                     by3=hby(ica(i+1),ica(i+2))
                     bz3=hbz(ica(i+1),ica(i+2))
                  else
                     bx3=ebx(i+2)
                     by3=eby(i+2)
                     bz3=ebz(i+2)
                  endif
                  dot23=bx2*bx3+by2*by3+bz2*bz3 !cos(H1.H2)
                  if(dot23.lt.-0.71)then !angle >135
                     ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
                     ESHORT5c=ESHORT5c-2-ff
                  endif
               endif
            endif
         endif
 2    continue
c^^^^^^^^^^E14, E_15, bias1 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^

****  Protein-like bias2: b(i), b(i+4) parallel------->
      i1=max(iiii-4,1)
      i2=min(jjjj-1,Lch-5)
      do i=i1,i2
         ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2 !ff=1, r<r0; =0.25, r>2r0.
         if(mv(i).gt.0)then
            bx1=hbx(ica(i-1),ica(i))
            by1=hby(ica(i-1),ica(i))
            bz1=hbz(ica(i-1),ica(i))
         else
            bx1=ebx(i)
            by1=eby(i)
            bz1=ebz(i)
         endif
         if(mv(i+4).gt.0)then
            bx5=hbx(ica(i+3),ica(i+4))
            by5=hby(ica(i+3),ica(i+4))
            bz5=hbz(ica(i+3),ica(i+4))
         else
            bx5=ebx(i+4)
            by5=eby(i+4)
            bz5=ebz(i+4)
         endif
         b15=bx1*bx5+by1*by5+bz1*bz5
         if(seq(i).eq.2.and.seq(i+2).eq.2.and.seq(i+4).eq.2)then
            if(b15.gt.0.9)then
               ESHORT5=ESHORT5-ff !alpha
            endif
         else
            if(b15.gt.0.5.or.b15.lt.-0.3)then !beta or turn
               ESHORT5=ESHORT5-ff
            endif
         endif
      enddo
c^^^^^^^^^^^^^^^^^^ bias2 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****  Protein-like bias3: c(i),c(i+2) anti/parallel:
      i1=max(iiii-2,1)
      i2=min(jjjj,Lch-2)
      do i=i1,i2
         ff=(afs(i)+afs(i+1)+afs(i+2))/3
         if(mv(i).gt.0)then
            cx1=cax(ica(i-1),ica(i))
            cy1=cay(ica(i-1),ica(i))
            cz1=caz(ica(i-1),ica(i))
         else
            cx1=ecx(i)
            cy1=ecy(i)
            cz1=ecz(i)
         endif
         if(mv(i+2).gt.0)then
            cx3=cax(ica(i+1),ica(i+2))
            cy3=cay(ica(i+1),ica(i+2))
            cz3=caz(ica(i+1),ica(i+2))
         else
            cx3=ecx(i+2)
            cy3=ecy(i+2)
            cz3=ecz(i+2)
         endif
         c13=abs(cx1*cx3+cy1*cy3+cz1*cz3)
         c13=min(0.71,c13)/0.71 !c13 is the same in [0,45]
         ESHORT5=ESHORT5-ff*c13
      enddo
c^^^^^^^^^^^^^^^^^^ bias3 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccccccbias4a: to predicted alpha fragment --------------->
      i1=max(iiii-6,1)
      i2=min(jjjj-1,Lch-7)
      do i=i1,i2
         if(frga(i).gt.1)then
            dis=sqrt((aax(i)-aax(i+7))**2+
     $           (aay(i)-aay(i+7))**2+(aaz(i)-aaz(i+7))**2)
            ESHORT5b=ESHORT5b+abs(dis-frga(i))
         endif
      enddo
ccccccbias4b: to predicted beta fragment --------------->
      i1=max(iiii-5,1)
      i2=min(jjjj-1,Lch-6)
      do i=i1,i2
         if(frgb(i).gt.1)then
            dis=sqrt((aax(i)-aax(i+6))**2+
     $           (aay(i)-aay(i+6))**2+(aaz(i)-aaz(i+6))**2)
            ESHORT5b=ESHORT5b+abs(dis-frgb(i))
         endif
      enddo
c^^^^^^^^^^^^^^^ predicted alpha/beta bias finished ^^^^^^^^^^^^^^^

ccc   penality for crumpling structures ---------------------------->
      i1=max(iiii-11,1)
      i2=min(jjjj-1,Lch-12)
      do i=i1,i2
         if(mv(i).gt.0)then
            ax0=x(i)
            ay0=y(i)
            az0=z(i)
         else
            ax0=ex(i)
            ay0=ey(i)
            az0=ez(i)
         endif
         if(mv(i+4).gt.0)then
            ax4=x(i+4)
            ay4=y(i+4)
            az4=z(i+4)
         else
            ax4=ex(i+4)
            ay4=ey(i+4)
            az4=ez(i+4)
         endif
         if(mv(i+8).gt.0)then
            ax8=x(i+8)
            ay8=y(i+8)
            az8=z(i+8)
         else
            ax8=ex(i+8)
            ay8=ey(i+8)
            az8=ez(i+8)
         endif
         avx1=ax4-ax0
         avy1=ay4-ay0
         avz1=az4-az0
         avx2=ax8-ax4
         avy2=ay8-ay4
         avz2=az8-az4
         aaa=avx1*avx2+avy1*avy2+avz1*avz2
         if(aaa.lt.0)then
            if(mv(i+12).gt.0)then
               ax12=x(i+12)
               ay12=y(i+12)
               az12=z(i+12)
            else
               ax12=ex(i+12)
               ay12=ey(i+12)
               az12=ez(i+12)
            endif
            avx3=ax12-ax8
            avy3=ay12-ay8
            avz3=az12-az8
            bbb=avx1*avx3+avy1*avy3+avz1*avz3
            if(bbb.gt.0)then
               ccc=avx2*avx3+avy2*avy3+avz2*avz3
               if(ccc.lt.0)then
                  ESHORT5a=ESHORT5a+1 !crumpling
               endif
            endif
         endif
      enddo
c^^^^^^^^^^^^^ penality of bizard structure finished (ESC1) ^^^^^^^^^^^^

c     Further penalize deriviation of restrain, if larger than (colim, dilim):
      IF(istat.eq.10) then      !calculate E from the beginning
cc         if(dinew.gt.dilim) ESHORT=ESHORT+er1*(dinew-dilim)
cc         if(conew.gt.colim) ESHORT=ESHORT+er3*(conew-colim)
         if(dinew.gt.dilim) ESHORT3a=dinew-dilim !dilim=number of restrains.
         if(conew.gt.colim) ESHORT4a=conew-colim !colim=number of restrains.
      endif
c     codevsum and didevsum are backuped as total deviation for old 
c     conformation after the first istat=10 are finished

c     the panalty was not counted in Enew, the following is the difference
c     of panlity on new and old conformation, i.e. dE=Enew-Eold. Since
c     this part energy appear only at Eold, the sign is oppsite:
      if(istat.lt.0) then       !return to old energy
         codev=codevsum+conew-coold !total deviation for new conformation
         didev=didevsum+dinew-diold !total deviation for new conformation
c     total penalty-energy beacuse of restrain-deviation on old conformation
cc         if(codevsum.gt.colim) ESHORT=ESHORT+er3*(codevsum-colim)
cc         if(didevsum.gt.dilim) ESHORT=ESHORT+er1*(didevsum-dilim)
         if(didevsum.gt.dilim) ESHORT3a=didevsum-dilim !E_old
         if(codevsum.gt.colim) ESHORT4a=codevsum-colim
c     total penalty-energy beacuse of restrain-deviation on new conformation:
cc         if(codev.gt.colim) ESHORT=ESHORT-er3*(codev-colim)
cc         if(didev.gt.dilim) ESHORT=ESHORT-er1*(didev-dilim)
         if(didev.gt.dilim) ESHORT3a=ESHORT3a-(didev-dilim) !E_old-E_new
         if(codev.gt.colim) ESHORT4a=ESHORT4a-(codev-colim)
      endif

      ESHORT=
c     $     es1*ESHORT1
     $     +es2*ESHORT2
     $     +er1*ESHORT3
     $     +er2*ESHORT3a
     $     +er3*ESHORT4
     $     +er4*ESHORT4a
     $     +es3*ESHORT5
     $     +es3a*ESHORT5a
     $     +es3b*ESHORT5b
     $     +es3c*ESHORT5c
     $     +es4*ESHORT6
     $     +es5*ESHORT7
     $     +es6*ESHORT8

c ^^^^^^^^^^ E_short finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  2-bond movement: prefabricated
c  residues of 'm, m1=m+1,m2=m+2' are involved,
c  but positions of m1 are changed.
c  old ---> new ---> old ----> new
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move2
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      double precision bNa,bNt,bNNa,bNNt
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy(no)*nfl2)+1  ![1,nfl2]
      m=ras2(i)
      m1=m+1
      m2=m+2

cccc  all the pairs have at least one another pair, so nc>=2
      nc=Np2(ica(m),ica(m1))    !number of pathes from m to m2

c     choose p'th new path from (m) to (m2) -------------->
 201  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m))) goto 201
      nn(m1)=v22(ica(m),ica(m1),p)
      if(.not.goodc(nn(m1),ica(m2))) goto 201
c^^^^^^^^^^ new conformation chosen finished ^^^^^^^^^^^^^^^^^^^^^^^^
      if(nn(m).eq.ica(m))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(m)=nn(m)
      ica(m1)=nn(m1)

      if(look(m,m2))then        ! check excluded volumn for passage of [m,m2]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m2,1)+ESHORT(m,m2,1) !icnt,nop are repeated.
         do kkk=m,m2
            afsn(kkk)=afs(kkk)  !afs(i) is useful in ESHORT
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         Eold=EHB(m,m2,-1)+ESHORT(m,m2,-1) !repeated part of icnt, nop removed

c     calculate eprofn while dord was calculated when call EHB(m,m2,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)          
            ia=noa(pp)          
            ip=nop(pp)          !nopp(old)+1(new)-1(old)=new
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(2)=bNa(2)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt          !backup1
            sumcto=sumct        !backup2
            do pp=1,Lch
               nopp(pp)=nop(pp) !backup3
               nomm(pp)=nom(pp) !backup4
               noaa(pp)=noa(pp) !backup5
            enddo
            eprofo=eprofn       !backup6
            codevsum=codev      !backup7
            didevsum=didev      !backup8
            do kkk=m,m2
               afs(kkk)=afsn(kkk) !backup9
            enddo

c            energ=energ+de
c     change the conformation to the new position--------->
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            call get_center     !update (amx,amy,amz)
            call get_axis
         endif                  !for id
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
      endif                     !for look(i,m)
 113  continue
      bNt(2)=bNt(2)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move2 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3-bond movement: single-move
c  residues of 'i,i+1,i+2,i+3' are involved,
c  but positions bond-vectors of i+1, i+2 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move3s
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy(no)*nfl3)+1
      m=ras3(i)
      m1=m+1
      m2=m+2
      m3=m+3

ccccccccccccc 1th 2-bond movement cccccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccc 2th 2-bond movement cccccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      nn(m2)=v22(t1,ica(m2),p)
      if(.not.goodc(nn(m2),ica(m3)))goto 202
c^^^^^^^^ three-bonds movement conformation finished ^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)

      if(look(m,m3))then        ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m3,1)+ESHORT(m,m3,1) !use ifs
         do kkk=m,m3
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         Eold=EHB(m,m3,-1)+ESHORT(m,m3,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(3)=bNa(3)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m3
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
      endif                     !for look(i,m)
 113  continue
      bNt(3)=bNt(3)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move3s finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3-bond movement: double move
c  residues of 'i,i+1,i+2,i+3' are involved,
c  but positions bond-vectors of i+1, i+2 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move3d
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy(no)*nfl3)+1
      m=ras3(i)
      m1=m+1
      m2=m+2
      m3=m+3

ccccccccccccc temporal 1th 2-bond movement cccccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(no)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccc temporal 2th 2-bond movement cccccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(no)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202

ccccccccccccc 1th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 203  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 203
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,t2))goto 203
ccccccccccccc 2th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,t2)
 204  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,t2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 204
      nn(m2)=v22(ttt1,t2,p)
      if(.not.goodc(nn(m2),ica(m3)))goto 204
c^^^^^^^^ three-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1))goto 113


c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)

      if(look(m,m3))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m3,1)+ESHORT(m,m3,1) !use ifs
         do kkk=m,m3
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         Eold=EHB(m,m3,-1)+ESHORT(m,m3,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(4)=bNa(4)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m3
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
      endif                     !for look(i,m)
113   continue
      bNt(4)=bNt(4)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move3d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     4-bond movement: single move
c     residues of 'i,i+1,i+2,i+3,i+4' are involved,
c     but positions bond-vectors of i+1, i+2, i+3 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move4s
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy(no)*nfl4)+1
      m=ras4(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m2)=v21(t2,ica(m3),p)
      if(.not.goodc(nn(m1),nn(m2)))goto 203
      nn(m3)=v22(t2,ica(m3),p)
      if(.not.goodc(nn(m3),ica(m4)))goto 203
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).
     $     and.nn(m2).eq.ica(m2))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)

      if(look(m,m4))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m4,1)+ESHORT(m,m4,1) !use ifs
         do kkk=m,m4
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         Eold=EHB(m,m4,-1)+ESHORT(m,m4,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(5)=bNa(5)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m4
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
      endif                     !for look(i,m)
113   continue
      bNt(5)=bNt(5)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move4s finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     4-bond movement: double move
c     residues of 'i,i+1,i+2,i+3,i+4' are involved,
c     but positions bond-vectors of i+1, i+2, i+3 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move4d
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy(no)*nfl4)+1
      m=ras4(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4

ccccccccccccccc 1th temperor 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(no)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(no)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(no)*nc)+1    ![1,nc]
      tt2=v21(t2,ica(m3),p)
      if(.not.goodc(tt1,tt2))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203

ccccccccccccc first 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 204  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 204
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,tt2))goto 204
ccccccccccccc second 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,tt2)
 205  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,tt2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 205
      ttt2=v22(ttt1,tt2,p)
      if(.not.goodc(ttt2,t3))goto 205
ccccccccccccc third 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt2,t3)
 206  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m2)=v21(ttt2,t3,p)
      if(.not.goodc(nn(m1),nn(m2)))goto 206
      nn(m3)=v22(ttt2,t3,p)
      if(.not.goodc(nn(m3),ica(m4)))goto 206
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).
     $     and.nn(m2).eq.ica(m2))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)

      if(look(m,m4))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m4,1)+ESHORT(m,m4,1) !use ifs
         do kkk=m,m4
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         Eold=EHB(m,m4,-1)+ESHORT(m,m4,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(6)=bNa(6)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m4
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
      endif                     !for look(i,m)
113   continue
      bNt(6)=bNt(6)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move4d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5-bond movement: single move
c     residues of 'i,i+1,i+2,i+3,i+4' are involved,
c     but positions bond-vectors of i+1, i+2, i+3 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move5s
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy(no)*nfl5)+1
      m=ras5(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4
      m5=m+5

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m2)=v21(t2,ica(m3),p)
      if(.not.goodc(nn(m1),nn(m2)))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203
ccccccccccccccc 4th 2-bond movement cccccccccccccccc
      nc=Np2(t3,ica(m4))
 204  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m3)=v21(t3,ica(m4),p)
      if(.not.goodc(nn(m2),nn(m3)))goto 204
      nn(m4)=v22(t3,ica(m4),p)
      if(.not.goodc(nn(m4),ica(m5)))goto 204
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).and.
     $     nn(m2).eq.ica(m2).and.nn(m3).eq.ica(m3))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      ox(m4)=x(m4)              !memory of old path
      oy(m4)=y(m4)              !memory of old path
      oz(m4)=z(m4)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      oo(m4)=ica(m4)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path
      nx(m4)=nx(m3)+vx(nn(m3))  !memory of new path
      ny(m4)=ny(m3)+vy(nn(m3))  !memory of new path
      nz(m4)=nz(m3)+vz(nn(m3))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      x(m4)=nx(m4)
      y(m4)=ny(m4)
      z(m4)=nz(m4)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)
      ica(m4)=nn(m4)

      if(look(m,m5))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m5,1)+ESHORT(m,m5,1) !use ifs
         do kkk=m,m5
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
         Eold=EHB(m,m5,-1)+ESHORT(m,m5,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(8)=bNa(8)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m5
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            x(m4)=nx(m4)
            y(m4)=ny(m4)
            z(m4)=nz(m4)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            ica(m4)=nn(m4)
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
      endif                     !for look(i,m)
113   continue
      bNt(8)=bNt(8)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move5s finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5-bond movement: double move
c     residues of 'i,i+1,i+2,i+3,i+4,i+5' are involved,
c     but positions bond-vectors of i+1, i+2, i+3, i+4 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move5d
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy(no)*nfl5)+1
      m=ras5(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4
      m5=m+5

ccccccccccccccc temperor 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(no)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc temperor 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(no)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc temperor 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(no)*nc)+1    ![1,nc]
      tt2=v21(t2,ica(m3),p)
      if(.not.goodc(tt1,tt2))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203
ccccccccccccccc temperor 4th 2-bond movement cccccccccccccccc
      nc=Np2(t3,ica(m4))
 204  p=int(aranzy(no)*nc)+1    ![1,nc]
      tt3=v21(t3,ica(m4),p)
      if(.not.goodc(tt2,tt3))goto 204
      t4=v22(t3,ica(m4),p)
      if(.not.goodc(t4,ica(m5)))goto 204

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(t,tt1)
 205  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 205
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,tt2))goto 205
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(ttt1,tt2)
 206  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,tt2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 206
      ttt2=v22(ttt1,tt2,p)
      if(.not.goodc(ttt2,tt3))goto 206
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(ttt2,tt3)
 207  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m2)=v21(ttt2,tt3,p)
      if(.not.goodc(nn(m1),nn(m2)))goto 207
      ttt3=v22(ttt2,tt3,p)
      if(.not.goodc(ttt3,t4))goto 207
ccccccccccccccc 4th 2-bond movement cccccccccccccccc
      nc=Np2(ttt3,t4)
 208  p=int(aranzy(no)*nc)+1    ![1,nc]
      nn(m3)=v21(ttt3,t4,p)
      if(.not.goodc(nn(m2),nn(m3)))goto 208
      nn(m4)=v22(ttt3,t4,p)
      if(.not.goodc(nn(m4),ica(m5)))goto 208
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).and.
     $     nn(m2).eq.ica(m2).and.nn(m3).eq.ica(m3))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      ox(m4)=x(m4)              !memory of old path
      oy(m4)=y(m4)              !memory of old path
      oz(m4)=z(m4)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      oo(m4)=ica(m4)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path
      nx(m4)=nx(m3)+vx(nn(m3))  !memory of new path
      ny(m4)=ny(m3)+vy(nn(m3))  !memory of new path
      nz(m4)=nz(m3)+vz(nn(m3))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      x(m4)=nx(m4)
      y(m4)=ny(m4)
      z(m4)=nz(m4)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)
      ica(m4)=nn(m4)

      if(look(m,m5))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m5,1)+ESHORT(m,m5,1) !use ifs
         do kkk=m,m5
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
         Eold=EHB(m,m5,-1)+ESHORT(m,m5,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(9)=bNa(9)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m5
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            x(m4)=nx(m4)
            y(m4)=ny(m4)
            z(m4)=nz(m4)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            ica(m4)=nn(m4)
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
      endif                     !for look(i,m)
113   continue
      bNt(9)=bNt(9)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move5d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     6-bond movement: translation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move6
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

c      nob=6+int(aranzy(no)*6.999) !number of moved bonds, [6,12]
      nob=6

c     choose the position (m) to be moved
      i=int(aranzy(no)*nfl6)+1 ![1,nfl6]
      m=ras6(i)
      m1=m+1
      mnob=m+nob
      mnob1=mnob-1
      mnob2=mnob-2

 201  mx=int(aranzy(no)*2.999999)-1 ![-1,0,1]
      my=int(aranzy(no)*2.999999)-1 ![-1,0,1]
      mz=int(aranzy(no)*2.999999)-1 ![-1,0,1]
      if((mx*mx+my*my+mz*mz).eq.0) goto 201

      wx=vx(ica(m))+mx
      wy=vy(ica(m))+my
      wz=vz(ica(m))+mz
      ir=wx*wx+wy*wy+wz*wz
      if(ir.lt.14)goto 201
      if(ir.gt.25)goto 201
      nn(m)=vector(wx,wy,wz)    !new vector of vertix 'm'
      if(.not.goodc(ica(m-1),nn(m))) goto 202
      if(.not.goodc(nn(m),ica(m1))) goto 202

      wx=vx(ica(mnob1))-mx
      wy=vy(ica(mnob1))-my
      wz=vz(ica(mnob1))-mz
      ir=wx*wx+wy*wy+wz*wz
      if(ir.lt.14) goto 202
      if(ir.gt.25) goto 202
      nn(mnob1)=vector(wx,wy,wz) !new vector of vertix 'm+5'
      if(.not.goodc(ica(mnob2),nn(mnob1))) goto 202
      if(.not.goodc(nn(mnob1),ica(mnob))) goto 202

c     prepare new path and backup old path ------------>
      do i=m1,mnob1
         ox(i)=x(i)             !memory of old path
         oy(i)=y(i)             !memory of old path
         oz(i)=z(i)             !memory of old path
         nx(i)=x(i)+mx          !memory of new path
         ny(i)=y(i)+my          !memory of new path
         nz(i)=z(i)+mz          !memory of new path
      enddo
      oo(m)=ica(m)              !memory of old path
      oo(mnob1)=ica(mnob1)      !memory of old path

c     change conformation to new path -------------->
      do i=m1,mnob1
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
      enddo
      ica(m)=nn(m)
      ica(mnob1)=nn(mnob1)

      if(look(m,mnob))then      ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,mnob,1)+ESHORT(m,mnob,1) !use ifs
         do kkk=m,mnob
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=m1,mnob1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
         enddo
         ica(m)=oo(m)
         ica(mnob1)=oo(mnob1)
         Eold=EHB(m,mnob,-1)+ESHORT(m,mnob,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(10)=bNa(10)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,mnob
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            do i=m1,mnob1
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
            enddo
            ica(m)=nn(m)
            ica(mnob1)=nn(mnob1)
            call get_center     !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         do i=m1,mnob1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
         enddo
         ica(m)=oo(m)
         ica(mnob1)=oo(mnob1)
      endif                     !for look(i,m)

 202  continue
      bNt(10)=bNt(10)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move6 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     random walk in either N-terminal 
c     m residues of '1, 2, ..., m' are relocated by random walk.
c     Maximum number of moving point = Mend_N
c     old ---> new ---> old ----> new
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move_n_end
      implicit integer(i-z)
      parameter(ndim=1000)
                parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      double precision bNa,bNt,bNNa,bNNt
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

      m=int(aranzy(no)*(Mend_N-1))+1 !m=[1,Mend_N-1], [1,m] will be moved.

c     back-up old path---------->
      do i=1,m
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path and move to new path -------------->
      do i=m,1,-1
 111     iv=int(aranzy(no)*anvec)+1
         if(.not.goodc(iv,ica(i+1))) goto 111
         x(i)=x(i+1)-vx(iv)
         y(i)=y(i+1)-vy(iv)
         z(i)=z(i+1)-vz(iv)
         ica(i)=iv
         nx(i)=x(i)             !memory of new conformation
         ny(i)=y(i)             !memory of new conformation
         nz(i)=z(i)             !memory of new conformation
         nn(i)=ica(i)           !memory of new conformation
      enddo
      ica(0)=ica(2)

      if(look(1,m+1))then       ! check excluded volumn for passage of [2,m+1]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(1,m+1,1)+ESHORT(1,m+1,1) !use ifs
         do kkk=1,m+1
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
         
c     return back the conformation and calculate E_old --------->
         do i=1,m
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         Eold=EHB(1,m+1,-1)+ESHORT(1,m+1,-1) !cut off old nop by istat=-1

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(11)=bNa(11)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
c            energ=energ+de
            do i=1,m
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            ica(0)=ica(2)
            do kkk=1,m+1
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev
            call get_center     !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
         do i=1,m
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
      endif                     !for look(i,m)
113   continue
      bNt(11)=bNt(11)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move_n_end finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     random walk in either C-terminal 
c     Maximum number of moving point = Mend_C
c     old ---> new ---> old ----> new
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move_c_end
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      double precision bNa,bNt,bNNa,bNNt
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

      m=int(aranzy(no)*(Mend_C-1))+1 !m=[1,Mend_C-1], m points will be moved
      mm=Lch-m+1                !(mm,Lch] will be relocated
c     mm-1 should be moveable residue, ica(m-2) exists.

c     backup old path----------->
      do i=mm,Lch
         ox(i)=x(i)             !memory of old conformation
         oy(i)=y(i)             !memory of old conformation
         oz(i)=z(i)             !memory of old conformation
         oo(i-1)=ica(i-1)       !memory of old conformation
      enddo
c     prepare new path------------->
      do i=mm,Lch
 111     iv=int(aranzy(no)*anvec)+1
         if(.not.goodc(ica(i-2),iv)) goto 111
         x(i)=x(i-1)+vx(iv)
         y(i)=y(i-1)+vy(iv)
         z(i)=z(i-1)+vz(iv)
         ica(i-1)=iv
         nx(i)=x(i)             !memory of new conformation
         ny(i)=y(i)             !memory of new conformation
         nz(i)=z(i)             !memory of new conformation
         nn(i-1)=ica(i-1)       !memory of new conformation
      enddo
      ica(Lch)=ica(Lch2)

      if(look(mm-1,Lch))then    ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(mm-1,Lch,1)+ESHORT(mm-1,Lch,1) !use ifs
         do kkk=mm-1,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=mm,Lch
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i-1)=oo(i-1)
         enddo
         ica(Lch)=ica(Lch2)
         Eold=EHB(mm-1,Lch,-1)+ESHORT(mm-1,Lch,-1) !nop go to new

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(12)=bNa(12)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
c     energ=energ+de
            do i=mm,Lch
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i-1)=nn(i-1)
            enddo
            ica(Lch)=ica(Lch2)

            do kkk=mm-1,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
         do i=mm,Lch
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i-1)=oo(i-1)
         enddo
         ica(Lch)=ica(Lch2)
      endif                     !for look(i,m)
113   continue
      bNt(12)=bNt(12)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move_c_end finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Sequence down-shift, i12<i21
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move7a
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      mm=12

      do i=1,lenf
         i21=4+int(aranzy(no)*(Lch-6)) ![4,lenf-3]
         if(u21(ica(i21),ica(i21+1)).ne.0)goto 31
      enddo
      goto 33                   !can not find a two-bond can merg into one-bond

 31   i0=max(2,i21-mm)
      do i=1,5
         i12=i0+int(aranzy(no)*((i21-i0)-1.0001)) ![i0,i21-2]
         nc=m12(ica(i12))       !number of ways of extending 1-->2
         if(nc.gt.0)goto 32
      enddo
      goto 33                   !cannot find a one-bond can extend to two-bond

c     ------------- 2 -> 1 ----------------------------
 32   nn(i21+1)=u21(ica(i21),ica(i21+1))
      if(.not.goodc(ica(i21-1),nn(i21+1)))goto33
      if(.not.goodc(nn(i21+1),ica(i21+2)))goto33
c     ------------- 1 -> 2 ----------------------------
      p=int(aranzy(no)*nc)+1
      nn(i12)=u1(ica(i12),p)
      nn(i12+1)=u2(ica(i12),p)
      if(.not.goodc(ica(i12-1),nn(i12)))goto33
      if(.not.goodc(nn(i12+1),ica(i12+1)))goto33

c     backup old path ------------>
      oo(i12)=ica(i12)
      do i=i12+1,i21+1
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path ------------>
      do i=i12+2,i21
         nn(i)=ica(i-1)
      enddo
      nx(i12+1)=x(i12)+vx(nn(i12))
      ny(i12+1)=y(i12)+vy(nn(i12))
      nz(i12+1)=z(i12)+vz(nn(i12))
      do i=i12+2,i21+1
         nx(i)=nx(i-1)+vx(nn(i-1))
         ny(i)=ny(i-1)+vy(nn(i-1))
         nz(i)=nz(i-1)+vz(nn(i-1))
      enddo

c     change conformation to new path -------------->
      ica(i12)=nn(i12)
      do i=i12+1,i21+1
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo

      m1=i12       !fixed
      m2=i21+2     !fixed
      if(look(m1,m2))then       ! check excluded volumn for passage of [m1,m2]
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs
         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(i12)=oo(i12)
         do i=i12+1,i21+1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(13)=bNa(13)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            ica(i12)=nn(i12)
            do i=i12+1,i21+1
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            call get_center     !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(i12)=oo(i12)
         do i=i12+1,i21+1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
      endif                     !for look(i,m)

 33   continue
      bNt(13)=bNt(13)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move7a finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Sequence up-shift, i12>i21
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move7b
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      mm=12

      do i=1,lenf
         i21=2+int(aranzy(no)*(Lch-6)) ![2,lenf-5]
         if(u21(ica(i21),ica(i21+1)).ne.0)goto 31
      enddo
      goto 33                   !can not find a two-bond can merg into one-bond

 31   i0=min(lenf2,i21+mm)
      do i=1,5
         i12=i21+3+int(aranzy(no)*(i0-i21-2.0001)) ![i21+3,i0]
         nc=m12(ica(i12))       !number of ways of extending 1-->2
         if(nc.gt.0)goto 32
      enddo
      goto 33                !can not find a one-bond can extend into two-bond

c     ------------- 2 -> 1 ----------------------------
 32   nn(i21)=u21(ica(i21),ica(i21+1))
      if(.not.goodc(ica(i21-1),nn(i21)))goto33
      if(.not.goodc(nn(i21),ica(i21+2)))goto33
c     ------------- 1 -> 2 ----------------------------
      p=int(aranzy(no)*nc)+1
      nn(i12-1)=u1(ica(i12),p)
      nn(i12)=u2(ica(i12),p)
      if(.not.goodc(ica(i12-1),nn(i12-1)))goto33
      if(.not.goodc(nn(i12),ica(i12+1)))goto33

c     backup old path ------------>
      oo(i21)=ica(i21)
      do i=i21+1,i12
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path ------------>
      do i=i21+1,i12-2
         nn(i)=ica(i+1)
      enddo
      nx(i21+1)=x(i21)+vx(nn(i21))
      ny(i21+1)=y(i21)+vy(nn(i21))
      nz(i21+1)=z(i21)+vz(nn(i21))
      do i=i21+2,i12
         nx(i)=nx(i-1)+vx(nn(i-1))
         ny(i)=ny(i-1)+vy(nn(i-1))
         nz(i)=nz(i-1)+vz(nn(i-1))
      enddo

c     change conformation to new path -------------->
      ica(i21)=nn(i21)
      do i=i21+1,i12
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo

      m1=i21       !fixed
      m2=i12+1     !fixed
      if(look(m1,m2))then       ! check excluded volumn for passage of [m1,m2]
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs
         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(i21)=oo(i21)
         do i=i21+1,i12
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(14)=bNa(14)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            ica(i21)=nn(i21)
            do i=i21+1,i12
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            call get_center     !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(i21)=oo(i21)
         do i=i21+1,i12
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
      endif                     !for look(i,m)

 33   continue
      bNt(14)=bNt(14)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move7b finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  q-bond random walk movement:  (Mran1 <= q <= Mran2)
c  residues of 'i,i+1, ..., i+q,i+q+1' are involved, 
c  but positions bond-vectors of i+1, i+2, ..., i+q are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move8
      implicit integer(i-z)
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/maxdis2/maxdis2(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      common/move31/v31(-15:15,-15:15,-15:15,6)
      common/move32/v32(-15:15,-15:15,-15:15,6)
      common/move33/v33(-15:15,-15:15,-15:15,6)
      common/move34/Np3(-15:15,-15:15,-15:15)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

c     q bonds will be moved, q is in [Mran1,Mran2]
      Mran1=4
      M_ran2=6
      q=int(aranzy(no)*(Mran2-Mran1+1))+Mran1
      q3=q-3

c     choose the position (m1) to be moved
c     m1>1, m2<lenf, otherwise can not do goodc()!
c     [m1+1,m2-1] will be really moved ---------------->
 111  m1=int(aranzy(no)*(lenf-q-2.00001))+2 !m1=[2,lenf-q-1]
      m2=m1+q                   !m2<=lenf-1

      mx=x(m2)-x(m1)
      my=y(m2)-y(m1)
      mz=z(m2)-z(m1)
cccccccccccccccFirst q-3 bonds cccccccccccccccccccccccccc
      nn(m1-1)=ica(m1-1)
      do i=1,q3
         uu1=0
 112     uu=1+int(nvec*aranzy(no)) ! uu in [1,nvec]
         uu1=uu1+1
         if(uu1.gt.nvec)goto 111 ! path maybe wrong!
         xx=mx-vx(uu)
         yy=my-vy(uu)
         zz=mz-vz(uu)

         xx2=xx*xx
         if(xx2.ge.maxdis2(q-i))goto 112
         yy2=yy*yy
         if(yy2.ge.maxdis2(q-i))goto 112
         rr=xx2+yy2+zz*zz   ! rr=xx*xx+yy*yy+zz*zz
         if(rr.ge.maxdis2(q-i))goto 112 !too distant for remaining walks.
         if(rr.lt.12)goto 112 !too close for overlap
         if(.not.goodc(uu,nn(m1+i-2)))goto 112 !check neighbor

         nn(m1+i-1)=uu
         mx=xx
         my=yy
         mz=zz
      enddo
      if(mx.gt.12.or.mx.lt.-12) goto 111 !no defined vector in 3-bond move
      if(my.gt.12.or.my.lt.-12) goto 111 !no defined vector in 3-bond move
      if(mz.gt.12.or.mz.lt.-12) goto 111 !no defined vector in 3-bond move
      
ccccccccccccccccccc Last 3 bonds ccccccccccccccccccccccccc
      nc=Np3(mx,my,mz)          !number of pathes from i to i+m ???
      if(nc.eq.0)then
         write(*,*)'absent q-movement in move5',q,mx,my,mz,m1,m2
         goto 114
      endif
      uu1=0
 113  p=int(aranzy(no)*(nc-0.00001))+1    ![1,nc]
      uu1=uu1+1
      if(uu1.gt.nc)goto 111
      nn(m2-3)=v31(mx,my,mz,p)
      if(.not.goodc(nn(m2-4),nn(m2-3))) goto 113
      nn(m2-1)=v33(mx,my,mz,p)  !???
      if(.not.goodc(nn(m2-1),ica(m2))) goto 113
      nn(m2-2)=v32(mx,my,mz,p)  !???

c     backup old path ------------>
      do i=1,q
         ox(m1+i)=x(m1+i)
         oy(m1+i)=y(m1+i)
         oz(m1+i)=z(m1+i)
         oo(m1+i-1)=ica(m1+i-1)
      enddo
c     prepare new path ------------>
      nx(m1+1)=x(m1)+vx(nn(m1))
      ny(m1+1)=y(m1)+vy(nn(m1))
      nz(m1+1)=z(m1)+vz(nn(m1))
      do i=2,q
         nx(m1+i)=nx(m1+i-1)+vx(nn(m1+i-1))
         ny(m1+i)=ny(m1+i-1)+vy(nn(m1+i-1))
         nz(m1+i)=nz(m1+i-1)+vz(nn(m1+i-1))
      enddo

c     change conformation to new path -------------->
      do i=1,q
         x(m1+i)=nx(m1+i)
         y(m1+i)=ny(m1+i)
         z(m1+i)=nz(m1+i)
         ica(m1+i-1)=nn(m1+i-1)
      enddo

      if(look(m1,m2))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
c            nhbn(pp)=nhbnn(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs
         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=1,q
            x(m1+i)=ox(m1+i)
            y(m1+i)=oy(m1+i)
            z(m1+i)=oz(m1+i)
            ica(m1+i-1)=oo(m1+i-1)
         enddo
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)        !noa is now 
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+eprofn-eprofo
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.3)then        !rejected
            icnt=icnto
            sumct=sumcto 		
         else                   !accepted
            bNa(7)=bNa(7)+1
            bN5a(q)=bN5a(q)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
c               nhbnn(pp)=nhbn(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            do i=1,q
               x(m1+i)=nx(m1+i)
               y(m1+i)=ny(m1+i)
               z(m1+i)=nz(m1+i)
               ica(m1+i-1)=nn(m1+i-1)
            enddo
            call get_center       !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         do i=1,q
            x(m1+i)=ox(m1+i)
            y(m1+i)=oy(m1+i)
            z(m1+i)=oz(m1+i)
            ica(m1+i-1)=oo(m1+i-1)
         enddo
      endif                     !for look(i,m)
 114  continue
      bNt(7)=bNt(7)+1
      bN5t(q)=bN5t(q)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move8 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     globe movement, there is biase.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move9
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1000)
      parameter(nvec=312)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1

      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension fax(ndim),fay(ndim),faz(ndim)
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

******backup old conformation --------------->
      do i=1,lenf
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo

******prepare new confromation ---------------->
c     (x,y,z)->(fax,fay,faz) based on [x(k),y[k],z[k]]:
      k=int(lenf*aranzy(no))+1   ![1,lenf]
      ddxyz=1.5
 99   fdx=(aranzy(no)*ddxyz*2.0)-ddxyz ![-1.5, 1.5]
      fdy=(aranzy(no)*ddxyz*2.0)-ddxyz
      fdz=(aranzy(no)*ddxyz*2.0)-ddxyz
      ar=fdx*fdx+fdy*fdy+fdz*fdz
      if(ar.lt.0.75) go to 99
      do i=1,lenf
         ar=sqrt(float((x(k)-x(i))**2+(y(k)-y(i))**2+(z(k)-z(i))**2))
         fax(i)=x(i)+fdx*(1.0-ar/(1.5*acrit))
         fay(i)=y(i)+fdy*(1.0-ar/(1.5*acrit))
         faz(i)=z(i)+fdz*(1.0-ar/(1.5*acrit))
      enddo

******project the conformation onto lattices, i.e. (fax,fay,faz)->(nx,ny,nz):
      px=nint(fax(1))
      py=nint(fay(1))
      pz=nint(faz(1))
      nx(1)=px
      ny(1)=py
      nz(1)=pz
      DO 101 i=2,lenf
         armin=10000.
         kk=0
         do 1009 k=1,nvec
            if(i.ne.2) then
               if(.not.goodc(kkkk,k))goto 1009
            endif
            kx=px+vx(k)
            ky=py+vy(k)
            kz=pz+vz(k)
            bx=fax(i)-float(kx)
            by=fay(i)-float(ky)
            bz=faz(i)-float(kz)
            ar=bx*bx+by*by+bz*bz
            if(ar.lt.armin) then
               kk=k
               mx=kx
               my=ky
               mz=kz
               armin=ar	
            endif
 1009    continue
         kkkk=kk
         if(kk.EQ.0) then
            write(20,*)' 	ERROR in the CHAIN global move' 
            GO TO  113          !do not move
         endif
         nx(i)=mx
         ny(i)=my
         nz(i)=mz
         px=mx
         py=my
         pz=mz
 101  continue
***************** new ica(i):
      ic=0
      do i=1,lenf1
         j=i+1
         wx=nx(j)-nx(i)
         wy=ny(j)-ny(i)
         wz=nz(j)-nz(i)
         nn(i)=vector(wx,wy,wz)
         if(nn(i).ne.oo(i)) ic=ic+1
      enddo
      if(ic.eq.0) go to 113     !do not move
c^^^^^^^^^^^^^^^^^^ new conformation obtained ^^^^^^^^^^^^^^^^^^^^^^^^^

c     change conformation to new path -------------->
      do i=1,lenf
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo

      if(look(2,lenf1))then
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(2,lenf1,1)+ESHORT(2,lenf1,1)
         do kkk=2,lenf1
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=1,lenf
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         Eold=EHB(2,lenf1,-1)+ESHORT(2,lenf1,-1)
      
c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo
      
c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,atemp,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(15)=bNa(15)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 	
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=2,lenf1
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            do i=1,lenf
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            call get_center     !update (amx,amy,amz)
            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         do i=1,lenf
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
      endif
 113  continue
      bNt(15)=bNt(15)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move9 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccc Try to swap i1-replica and i2-replica ccccccccccccccccc
      subroutine swap(i1,i2)
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1000)	!maximum length of chain-length
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/nswap/bNSa(100),bNSt(100)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/chain0/ras(ndim),nfl
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      
      if(metro_swap(i1,i2).eq.1)then
         bNSa(i1)=bNSa(i1)+1    !number of accepted swaps
         do i=1,Lch
*     x:
            tempx=xrep(i,i1)
            tempy=yrep(i,i1)
            tempz=zrep(i,i1)
            xrep(i,i1)=xrep(i,i2)
            yrep(i,i1)=yrep(i,i2)
            zrep(i,i1)=zrep(i,i2)
            xrep(i,i2)=tempx
            yrep(i,i2)=tempy
            zrep(i,i2)=tempz
         enddo
***   swap E ----------------------------->
         attt=E_rep(i1)         !exchange E_rep
         E_rep(i1)=E_rep(i2)
         E_rep(i2)=attt
***   swap moveable points---------------->
      endif
      bNSt(i1)=bNSt(i1)+1       !number of total swaps.

      return
      end

ccccccccccccccc swap: id=1, accepted; id=3, rejected ccccccccccccccccc
      integer function metro_swap(i,j)
      parameter(nrep=100)
      common/sw1/aT_rep(nrep),E_rep(nrep)

c     w_12=wi(j)wj(i)/wi(i)wj(j) ------------->
      aaa=(1/aT_rep(i)-1/aT_rep(j))*(E_rep(i)-E_rep(j))
      aweight=exp(aaa)

      if(aranzy(no).le.aweight)then
         metro_swap=1		!swap accepted
      else
         metro_swap=3		!swap rejected
      endif

      return
      end

ccccccccccccccc Try to swap i1-replica and i2-replica ccccccccccccccccc
c     swap the replicas with combined chains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine swap_RS(i1,i2)
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1000)	!maximum length of chain-length
      parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/nswap/bNSa(100),bNSt(100)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      common/chain0/ras(ndim),nfl
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/chainm/mv(ndim)

      if(metro_swap(i1,i2).eq.1)then
         bNSa(i1)=bNSa(i1)+1    !number of accepted swaps
         do i=1,Lch
            t_ica=icarep(i,i1)             !ica(i)
            icarep(i,i1)=icarep(i,i2)
            icarep(i,i2)=t_ica
            if(mv(i).gt.0)then
               t_x=xrep(i,i1)              !x(i)
               t_y=yrep(i,i1)
               t_z=zrep(i,i1)
               xrep(i,i1)=xrep(i,i2)
               yrep(i,i1)=yrep(i,i2)
               zrep(i,i1)=zrep(i,i2)
               xrep(i,i2)=t_x
               yrep(i,i2)=t_y
               zrep(i,i2)=t_z
            else
               e_x=exrep(i,i1)           !ex(i)
               e_y=eyrep(i,i1)
               e_z=ezrep(i,i1)
               exrep(i,i1)=exrep(i,i2)
               eyrep(i,i1)=eyrep(i,i2)
               ezrep(i,i1)=ezrep(i,i2)
               exrep(i,i2)=e_x
               eyrep(i,i2)=e_y
               ezrep(i,i2)=e_z
               e_x=egxrep(i,i1)          !egx(i)
               e_y=egyrep(i,i1)
               e_z=egzrep(i,i1)
               egxrep(i,i1)=egxrep(i,i2)
               egyrep(i,i1)=egyrep(i,i2)
               egzrep(i,i1)=egzrep(i,i2)
               egxrep(i,i2)=e_x
               egyrep(i,i2)=e_y
               egzrep(i,i2)=e_z
               e_x=ecxrep(i,i1)          !ecx(i)
               e_y=ecyrep(i,i1)
               e_z=eczrep(i,i1)
               ecxrep(i,i1)=ecxrep(i,i2)
               ecyrep(i,i1)=ecyrep(i,i2)
               eczrep(i,i1)=eczrep(i,i2)
               ecxrep(i,i2)=e_x
               ecyrep(i,i2)=e_y
               eczrep(i,i2)=e_z
               e_x=ebxrep(i,i1)          !ebx(i)
               e_y=ebyrep(i,i1)
               e_z=ebzrep(i,i1)
               ebxrep(i,i1)=ebxrep(i,i2)
               ebyrep(i,i1)=ebyrep(i,i2)
               ebzrep(i,i1)=ebzrep(i,i2)
               ebxrep(i,i2)=e_x
               ebyrep(i,i2)=e_y
               ebzrep(i,i2)=e_z
            endif
         enddo
***   swap E ----------------------------->
         attt=E_rep(i1)         !exchange E_rep
         E_rep(i1)=E_rep(i2)
         E_rep(i2)=attt
***   swap moveable points---------------->
      endif
      bNSt(i1)=bNSt(i1)+1       !number of total swaps.

c^^^^^^^^^^^^^^^^^ swap_RS finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test whether all the neighboring backbones are good neighbors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine test_neighbor
      implicit integer(i-z)
		parameter(ndim=1000)
		parameter(nrep=100)
                parameter(nvec=312)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      do j=1,Lch
         i=j-1
         if(j.ne.1.and.j.ne.Lch) then
            ii=ica(i)
            jj=ica(j)
            if(.not.goodc(ii,jj)) then
           write(20,8112)i,j,vx(ii),vy(ii),vz(ii),vx(jj),vy(jj),vz(jj)
 8112          format(5x,'warning2 -wrong input chain - vectors ',8i4)
               stop
            endif
         endif
      enddo 

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test for overlaps of Ca:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine test_overlap
      implicit integer(i-z)
      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2

      do i=1,Lch
         do j=i+3,Lch
            dis2=(aax(i)-aax(j))**2+(aay(i)-aay(j))**2+
     $           (aaz(i)-aaz(j))**2
            if(dis2.lt.exc)then
               write(*,*)i,j,dis2,'  Ca-Ca overlap'
            endif
         enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate energy from the beginning
c     parameters will be modified after running this subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energy_tot()
      implicit integer(i-z)
		parameter(ndim=1000)
                parameter(nvec=312)
      common/lengths/Lch,Lch1,Lch2
      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/seqe/seq(ndim),sec(ndim)
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/pair1/eh2,eh1b
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      COMMON/RES/ ER3, Mcom(ndim),Kcom(ndim,50)
      common/otherenergy/E_cord,E_cnum

      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB2,EHB3,EHB4,EHB5,EHB6
      common/ehbenergy1/EHB5a,EHB5b

      ICNTO=0
      SUMCTO=0
      DO i=1,Lch
         NOP(i)=0
         NOA(i)=0
         NOM(i)=0
      enddo

      energy_tot=EHB(1,Lch,1)+ESHORT(1,Lch,10)

      eprofo=0.0
      DO k=1,Lch
         is=seq(k)
         ia=NOA(k)
         ip=NOP(k)
         im=NOM(k)
         eprofo=eprofo+envir(ia,im,ip,is,3)
      enddo

      E_cord=abs(ICNT/(0.00001+float(SUMCT))-acorder)
      E_cnum=abs(float(SUMCT)-contt)

      energy_tot=energy_tot+en1*eprofo+en2*E_cord+en3*E_cnum

c     following is for next movement ---------->
      icnto=icnt                !backup of contact-order
      sumcto=sumct              !backup of contact-order
      do i=1,Lch
         nopp(i)=nop(i)         !backup of number of contact-pair
         noaa(i)=noa(i)         !backup of number of contact-pair
         nomm(i)=nom(i)         !backup of number of contact-pair
      enddo
      codevsum=conew            !backup of panelity for total deviation of comb
      didevsum=dinew            !backup of panelity for total deviation of dist

c^^^^^^^^^^^^^^^^^^^^ Energy_tot finished ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc Metropolis: id=1, accepted; id=3, rejected ccccccccccccccccc
      subroutine metro(dE,T,id)
c     double precision weight

      id=1
      if(dE.gt.0)then
c     if(aranzy(no).gt.weight(Et,T)/weight(E,T))then
c     if(aranzy(no).gt.weight12(dE,T))then
         if(aranzy(no).gt.dexp(-dble(dE)/T))then
            id=3
         endif
      endif
      return
      end

ccccccccccccccccweight factor cccccccccccccccccccccccccccccccccccccccccccc
      function weight12(dE,atemp)
      weight12=exp(-arcsh(dE)/atemp) !exp(-dE/T)
      return
      end

ccccccccccccccc flat function ccccccccccccccccccccccccccc
      function arcsh(x)
      arcsh=log(x+sqrt(1+x*x))
      return
      end

cccccccccccCalculate dRMSD between 1 and CA cccccccccccccccccccccccccccccccccc
      subroutine cal_drmsd(nb,nd,drmsd)
      parameter(ndim=1000)      !maximum length of chain-length
      parameter(ndcy=100)       !maximum number of decoys in each bin
      parameter(nbin=20)        !maximum number of bins
      common/CA/nCA,axca(ndim),ayca(ndim),azca(ndim)
      common/dcyx/dex(nbin,ndcy,ndim)
      common/dcyy/dey(nbin,ndcy,ndim)
      common/dcyz/dez(nbin,ndcy,ndim)

      nn=nCA
      drmsd=0
      do i=1,nn
         do j=i+1,nn
            d1=sqrt((dex(nb,nd,i)-dex(nb,nd,j))**2
     $           +(dey(nb,nd,i)-dey(nb,nd,j))**2
     $           +(dez(nb,nd,i)-dez(nb,nd,j))**2)
            d2=sqrt((axca(i)-axca(j))**2
     $           +(ayca(i)-ayca(j))**2
     $           +(azca(i)-azca(j))**2)
            drmsd=drmsd+(d1-d2)**2
         enddo
      enddo
      drmsd=sqrt(drmsd/(nn*(nn-1)))

      return
      end


ccccccccccccccccc calculate RMSD and DRMS cccccccccccccccccccccccccccc
c x_data(i)=r_d(1,i); y_data(i)=r_d(2,i); z_data(i)=r_d(3,i);
c x_model(i)=r_m(1,i); y_model(i)=r_m(2,i); z_model(i)=r_m(3,i);
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zyrmsd(arms)
      implicit integer (i-z)
      parameter(ndim=1000)	!maximum length of chain-length
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/rmsdrange/nca1,nca2
      common/CA1/dx(ndim),dy(ndim),dz(ndim)
      double precision r_m(3,ndim),r_d(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms !armsd is real
      data w /ndim*1.0/

      common/CA/nCA,axca(ndim),ayca(ndim),azca(ndim)
      common/initial/x0(ndim),y0(ndim),z0(ndim)
      common/addition/i_gap,s_gap(100),n_add(100),iLch(ndim),inCA(ndim)


      k=0
      do i=nca1,nca2
         k=k+1
         r_m(1,k)=x(i+1)*0.87   !model conformation
         r_m(2,k)=y(i+1)*0.87
         r_m(3,k)=z(i+1)*0.87
         r_d(1,k)=dx(i)         !native conformation
         r_d(2,k)=dy(i)
         r_d(3,k)=dz(i)
      enddo

      nn=k                      !number of data points
      call u3b(w,r_m,r_d,nn,0,rms,u,t,ier)
      arms=dsqrt(rms/nn)	!RMSD is real, rms is double precision

      return
      end
        
cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(1), x(3, 1), y(3, 1), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     &e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     &, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     &ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     &, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     &ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     &, e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
c 156 "rms.for"
      wc = zero
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do 1 j = 1, 3
      d = zero
      if (i .eq. j) d = one
      u(i,j) = d
      a(i,j) = d
    1 r(i,j) = zero
      ier = -1
c**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
c 170 "rms.for"
      if (n .lt. 1) return 
c 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
    2 yc(i) = yc(i) + (w(m) * y(i,m))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** EIGENVALUES *****************************************
c**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
c 200 "rms.for"
    5 rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j)
     &)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     &) + (rr1 * rr3)) - (rr2 * rr2)) / three
c 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    6 e(i) = spur
c**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
c 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
c 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
c.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
c 224 "rms.for"
      if (mode) 10, 50, 10
c**************** EIGENVECTORS *****************************************
c 226 "rms.for"
    8 if (mode) 30, 50, 30
c 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
   14 d = d + (ss(k) * ss(k))
      if (d .gt. zero) d = one / dsqrt(d)
      do 15 i = 1, 3
   15 a(i,l) = a(i,l) * d
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   17 p = p + (a(i,m1) ** 2)
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
   18 a(i,m1) = a(i,m1) * p
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
c****************** ROTATION MATRIX ************************************
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     &))
c 288 "rms.for"
   31 d = d + (b(i,l) ** 2)
      if (d .gt. zero) d = one / dsqrt(d)
      do 32 i = 1, 3
   32 b(i,l) = b(i,l) * d
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
   33 p = p + (b(i,2) ** 2)
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
   34 b(i,2) = b(i,2) * p
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 39 j = 1, 3
c****************** TRANSLATION VECTOR *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     & * xc(3))
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   51 e(i) = dsqrt(e(i))
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end

****************************************************
*     no=-1 (less than 0)                          *
*     rrr=ranzy(no)                                *
****************************************************
      function aranzy(no)
c   random numbers uniformly distributed between 0 and 1.
c   (code by Y. Zhang, ITP, Acdemia Sincica, 1999.6
c   can be replaced by any other suitable generator for random numbers)
      common/asdfa/ix1,ix2,ix3,rm1,rm2,r(99)
      data ia1,ic1,m1/1279,351762,1664557/
      data ia2,ic2,m2/2011,221592,1048583/
      data ia3,ic3,m3/15551,6150,29101/
      if(no.ge.0) go to 2
      ix1=mod(-no,m1)
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ix1,m2)
      ix1=mod(ia1*ix1+ic1,m1)
      ix3=mod(ix1,m3)
      rm1=1./float(m1)
      rm2=1./float(m2)
      do 1 j=1,99
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
 1    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 2    ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(99*ix3)/m3
      aranzy=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      no=ix1
      return
      end
