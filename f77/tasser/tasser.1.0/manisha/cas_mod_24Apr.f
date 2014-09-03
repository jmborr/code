CM0123456789012345678901234567890123456789012345678901234567890123456789
c     pgf77 -Mextend -O -s -fast -Wl,-static -o cas cas.f

      program TASSER
      implicit integer(i-z)
c     maximum length of chain-length
      parameter(ndim=1500)
c     maximum number of replicas
      parameter(nrep=100)
c     number of vectors on lattice
      parameter(nvec=312)

      integer maxdis2(ndim)
      character*6 mname(100)
      integer sec(ndim),seq(ndim)
      character*3 sequ(ndim)
      integer icg(ndim)
      real eonekd(0:19),eonehw(0:19),acrit
      real envir(0:15,0:15,0:15,0:19,4)
      real csr(ndim,2)
      integer ibin(-300:300)
      real asr(ndim,-12:12)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      real atemp1,atemp2,exc,er1,er3,er4,er5,er6,er7,acorder
      real h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
      real bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6,bhendn
      real colim,dilim,coold,diold,codevsum,conew,didevsum,dinew,contt
      real en1,en2,en3
      real Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      real Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      real dist_CA_cut
      real es2,es3,es3a,es3b,es3c,es4,es5,es6
      real didev,codev
      real anvec
      real eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b
      real ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,ehb4,ehb5a,ehb5b

      integer Lch,ncycle,N_rep,Ndis,Ncom,random,switch,phot,L_cut
      integer Lch1,Lch2
      integer icnto,icnt,sumcto,sumct
      integer t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      integer i_bigbond,teco
      integer nfl2,nfl3,nfl4,nfl5,nfl6
      integer itemp,icycle
      integer Mran1,Mran2,Mend,Mend_N,Mend_C
      integer ichos
      integer n_repf



      character fn
      dimension JBIN(0:500)
      dimension vvv(ndim,ndim)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension Mcom(ndim),Kcom(ndim,100)
      dimension aweig(4000,4000)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension apar(ndim,ndim)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      dimension area(ndim)
      dimension aT_rep(nrep),aTs_rep(nrep),aTTs_rep(nrep)
      dimension ehbij(ndim,ndim)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      logical goodc(nvec,nvec)
      dimension u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)
      dimension angle(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension axalf(0:19),ayalf(0:19),azalf(0:19)
      dimension axbet(0:19),aybet(0:19),azbet(0:19)
      dimension frga(ndim),frgb(ndim)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      dimension mv(ndim)
      dimension cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)
      dimension ras(ndim)
      dimension d_xyz00(ndim),angle00(ndim)
      dimension nfr_i(ndim),nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension N_sum(nrep),energ_sum(nrep),energ_sum2(nrep)
      dimension bNSa(nrep),bNSt(nrep)
      dimension bNa(nrep),bNt(nrep),bNNa(nrep),bNNt(nrep)
      dimension bN5a(nrep),bN5t(nrep)
      dimension E_rep(nrep)
      dimension armsd_sum(nrep),N_rmsd(nrep)
      dimension icarep(ndim,nrep)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      dimension egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      dimension ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      dimension ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension AA(3,3),EE(3),HH(3,3)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension ras2(ndim)
      dimension ras3(ndim)
      dimension ras4(ndim)
      dimension ras5(ndim)
      dimension ras6(ndim)
      dimension cx(1000),cy(1000),cz(1000)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      dimension bsr(ndim,16)
      dimension acops(ndim,16)

      common/seed1/no
CM      common/asdfa/ix1,ix2,ix3,rm1,rm2,r(99)




CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/three/angle(nvec,nvec)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
CM      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
CM      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
CM      common/hopp/eonehw(0:19)
CM      common/looks/exc,exc1,exc2
CM      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
CM      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
CM      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc

CM      common/seed1/no
CM      common/arandom/  aarand,abrand,acrand,adrand
CM      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
CM      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
CM      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
CMc     safe when vr^2<30
CM      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
CM      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
CM      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
CM      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/fr/frga(ndim),frgb(ndim)
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
CM      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
CM      common/maxi/maxin,vect1,vect2
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/forpreparemove4/ asrr(0:19,0:19,-12:12)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
CM      common/distres/er4,es3c
CM      common/rmsdrange/nca1,nca2
CM      common/CA/dx(ndim),dy(ndim),dz(ndim)
CM      common/msichores/msicho
CM      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,EHB5,EHB6
CM      common/ehbenergy1/EHB5a,EHB5b
CM      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT4,ESHORT11
CM      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
CM      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT10
CM      common/otherenergy/E_cord,E_cnum
CM      common/resnumber/Ncom,Ndis,accur
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      common/thrtem/aT_ann(100)
CM      common/outputxyz/fxyz(3,ndim)
CM
CM      common/temperature/itemp,atemp
CM      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
CM      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
CM      common/pair1/eh2,eh1b,eh1c
CM      dimension E_s(nrep),E_ss(nrep)
CM      common/paircut/ash
CM      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
CM
CM      common/nswap/bNSa(100),bNSt(100)
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
CM      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap
CM
CM      common/moveratio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
CM      common/moveratio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
CM      common/moveratio3/bh7a,bh7b,bhendn,bhendc
CM      common/icgg/ icg(ndim), EH6  
CM      common/rs/i_thr0
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
CM      common/nrepfile/n_repf
CM
CM      common/commonuse1/random,ncycle
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/mng/m_g(100)
CM      common/acct/accept0
CM      character*6 mname
CM      character fn
CM      common/movename/mname(100)
CM      common/readinitial/m_initial
CM      common/eall/N_sum(100),energ_sum(100),energ_sum2(100),E_min
CM
CM      common/chain0/ras(ndim),nfl
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CMc     CA
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CMc     SG
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim)
CMc     cc
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
CMc     Hb
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)
CM      common/chainm/mv(ndim)
CM
CM      common/sw1/aT_rep(nrep),E_rep(nrep)
CM      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
CM      common/sw3/icarep(ndim,nrep)
CMc     Ca
CM      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
CMc     SG
CM      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
CMc     cc
CM      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
CMc     hb
CM      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
CMc     chuan denotes the energy scale
CM      common/weight/chuan
CM      common/bigbond/i_bigbond,teco
CM      common/ssp/ssp
CM      common/defoangle/defo_angle
CM      common/fractpair/fract_pair1,fract_pair3
CM      common/zscore/izscore
CM      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
CM      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
CM      common/ichos/ichos
CM      common/nana1/nana

ccccc common input files
c     cutoff of contact predi.
      open(unit=1, file='contact.comm',  status='old')
c     envir
      open(unit=2, file='profile3.comm', status='old')
c     SG-potential, format
      open(unit=3, file='quasi3.comm',   status='old')
c     for Sc position.
      open(unit=4, file='sidechain.comm',status='old')
c     E_short of (i,i+2)
      open(unit=5, file='r13.comm',  status='old')
c     E_short of (i,i+3)
      open(unit=12,file='r14.comm',  status='old')
c     stress helical-stru.
      open(unit=7, file='r14h.comm', status='old')
c     stress extended-stru.
      open(unit=8, file='r14e.comm', status='old')
c     E_short(i,i+4)
      open(unit=9, file='r15.comm',  status='old')
      open(unit=10,file='r15h.comm', status='old')
      open(unit=11,file='r15e.comm', status='old')

ccccc sequence specified input files
      open(unit=21,file='rmsinp',status='old')


ccccc file derived from mkseq.pl
      open(unit=14,file='seq.dat',     status='old')

ccccc files derived from mkdat.pl
c     SG-pair potential
      open(unit=15, file='par.dat',   status='old')
c     SG-contact restraints
      open(unit=16, file='comb.dat',  status='old')
c     CA-distant restraints
      open(unit=17, file='dist.dat',  status='old')
c     CA-contact restraints
      open(unit=18, file='combCA.dat',status='old')
c     long range CA dist restraint
      open(unit=22, file='distL.dat', status='old')
c     surface exposure propensity
      open(unit=27, file='exp.dat',   status='unknown')
c     template coordinates
      open(unit=24, file='chain.dat', status='old')

ccccc files derived from mkpair.pl
c     orientation specific
      open(unit=25, file='pair3.dat', status='unknown')
c     orientation independent
      open(unit=26, file='pair1.dat', status='unknown')

      open(unit=19, file='in.dd',     status='old')
      open(unit=20, file='out.d',     status='unknown')

      read(19,*) random,ncycle,phot,N_rep
      write(*,*) '1=',random,ncycle,phot,N_rep

      read(19,*) h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
      write(*,*) '2=',h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend

      read(19,*) atemp2,atemp1,exc,exc1,exc2,Mend,defo_angle
      write(*,*) '3=',atemp2,atemp1,exc,exc1,exc2,Mend,defo_angle

      read(19,*) d_xyz0,angle0,L_cut,teco
      write(*,*) '4=',d_xyz0,angle0,L_cut,teco

c     i_thr0: 0->consensus; 1->top-1; 2->2th
      read(19,*) switch,i_thr0,ssp
      write(*,*) '5=',switch,i_thr0,ssp

c     alpha-type HB
      read(19,*) eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      write(*,*) '6=',eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh

c     alpha-type HB
      read(19,*) eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      write(*,*) '7=',eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh

c     for EHBs (6)
      read(19,*) eh1,eh2,eh3,eh4
      write(*,*) '8=',eh1,eh2,eh3,eh4

      read(19,*) eh1a,eh1b
      write(*,*) '9=',eh1a,eh1b

c     for ESHORTs (6)
      read(19,*) es2,es3,es4,es5,es6
      write(*,*) '10=',es2,es3,es4,es5,es6

      read(19,*) es3a,es3b,es3c
      write(*,*) '11=',es3a,es3b,es3c

c     for ensemble (3)
      read(19,*) en1,en2,en3
      write(*,*) '12=',en1,en2,en3

      read(19,*) chuan
      write(*,*) '13=',chuan

c     SQ2
      read(19,*) aTs2,aTs1
      write(*,*) '14=',aTs2,aTs1

c     ARCSH
      read(19,*) aTTs2,aTTs1
      write(*,*) '15=',aTTs2,aTTs1

      read(19,*) nana
      write(*,*) '16=',nana


c     set common parameters
      call set_common(random,Lch,Lch1,Lch2,anvec,contt,maxdis2,
     +                      ichos,vvv,mk_pair3,t_comb,t_dist,t_combCA,
     +                      s_comb,s_dist,s_combCA,atemp1,atemp2,N_rep,
     +                      mname,izscore,er1,er3,er4,er5,eh1c,er6,er7,
     +                      chuan,ssp,i_bigbond,fract_pair1,fract_pair3,
     +                      teco)

                  print *,'After set_common NO=',no

c     read seq(i),sec(i) from 'seq.dat'
      call read_seq (Lch,sec,seq,sequ)
cc     read eoinp(i,dis) from 'centro.comm'
      call read_centro(acrit,Lch,eonekd,eonehw)
cc     read envir(ia,im,ip,i,j) from 'profile3.comm'
      call read_profile(envir)
cc     read 1-3 short-range E from 'r13.comm'
      call read_E13(csr,seq,Lch2)
c     read 1-4 potential from 'r14*.dat'
      call read_E14(Lch,asr,sec,seq,ibin)
c     read 1-5 potential from 'r15*.dat'
      call read_E15(Lch,bsr,sec,seq,jbin,acops)
c     read 'quarsi3.comm'
      call read_quarsi3(arla,arlm,arlp,alp,ala,alm,Lch,seq,
     +                  mk_pair3,i_pair3,fract_pair3,apa,apm,app)
c     read 'par.dat', pair-wise potential
      call read_par(Lch,mk_pair3,i_pair3,fract_pair1,chuan,eh1c,
     +              apar)
c     read cut-off for contact prediction
      call read_concut(izscore,concut,concut2)
c     read contact restrains from 'comb.dat'
      call read_contactrestrain(izscore,aweig,Mcom,Kcom,Lch,colim,
     +                          Ncom)
c     read distant restrains from 'dist.dat'
      call read_distantrestrain(Lch,Mdis,kdis,dist,dev,Ndis)
c     read long distant restrain from 'distL.dat'
      call read_longdistantrestrain(Lch,MdisL,kdisL,distL)
c     read CAcontact restrains from 'combCA.dat'
      call read_CAcontact(izscore,dist_CA_cut,aweigCA,Lch,McomCA,KcomCA)
c     read solvent expose prediction
      call read_exp(nana,Lch,area)
c     reset temperature according N_rest
      call reset_temperature(er1,er3,er4,Ncom,Lch,atemp2)
c     set temperature for different replicas
      call set_temperature(N_rep,aT_rep,atemp1,atemp2,aTs_rep,aTs1,
     +                     aTs2,aTTs_rep,aTTs1,aTTs2)
c     set structure-specitic H-bond, EHBIJ(i,j)
      call set_EHB(Lch,sec,ehbij)
c     prepare all possible bond-vectors
      call prepare_vectors(vx,vy,vz,vector)
c     define goodc(i,j), angle(i,j), prod(i,j)
      call prepare_neighbors(m12,u21,vx,vy,vz,angle,goodc,u1,u2)
c     define C_beta, C_group, and hydrogen-bond
      call prepare_beta(axalf,ayalf,azalf,axbet,aybet,azbet,vx,vy,vz,
     +                        hbx,hby,hbz,cax,cay,caz,angle,gx,gy,gz)
c     compute the secondary fragment biases
      call prepare_frg(Lch,frga,frgb,sec)
c     calculate contact order
      call get_acorder(Lch,sec,acorder)
c     print out initial parameters
                  print *,'Before write_para NO=',no
      call write_parameter(Lch,ncycle,N_rep,phot,Ndis,Ncom,atemp2,
     +                     atemp1,exc,er1,er3,er4,er5,er6,er7,eh1c,
     +                     random,acorder,switch)
c     set movement percentage
      call set_move_ratio(h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend,bh2,
     +                    bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6,bhendn)
c     two-bond move, num2=26784
      call prepare_move2(Nw,goodc,vx,vy,vz,Np2,v21,v22,vector,w21,w22)
c     number of output replicas
      n_repf=16                 

      print *,'switch',switch
c     case when some fragments are frozen
      if(switch.gt.1)then
         n_repf=14
c     initial model from templates
      call template_initial2(exc,i_thr0,Lch,q,cx0,cy0,cz0,ssp,M_a,
     +                             M_i,M_f,L_a,L_cut,n_i,n_f,cx,cy,cz,
     +                             ras,nfr_i,nfr_f,mv,angle00,angle0,
     +                             d_xyz00,d_xyz0,N_rep,exrep,eyrep,
     +                             ezrep,ecxrep,ecyrep,eczrep,ebxrep,
     +                             ebyrep,ebzrep,seq,axalf,ayalf,azalf,
     +                             axbet,aybet,azbet,egxrep,egyrep,
     +                             egzrep,ex,ey,ez,ax,ay,az,x,y,z,goodc,
     +                             vx,vy,vz,xrep,yrep,zrep,sec,
     +                             icarep,vector,sequ,vvv,ex0,ey0,ez0,
     +                             i_bigbond,
     +                             nfl,Mend_N,Mend,Mend_C,nfl6,
     +                             ras6,nfl5,ras5,nfl4,ras4,nfl3,ras3,
     +                             nfl2,ras2,Lch2)
      endif

      print *,"After template_initial2 nfl2=",nfl2

ccccc open trajectory files
      if(n_repf.gt.N_rep)n_repf=N_rep
      do i=1,n_repf
         if(i.lt.10)then
            fn=char(48+i)
            open(unit=30+i,file='rep'//fn//'.tra',status='unknown')
         else
            fn=char(48+(i-10))
            open(unit=30+i,file='rep1'//fn//'.tra',status='unknown')
         endif
      enddo

c     case when some fragments are frozen
      if(switch.gt.1)then
CM         call template_simulation !simulate the templates
      endif


c     read initial (x,y,z) from 'chain.dat'
CM      call read_initial2        
      call read_initial2(Lch,q,cx0,cy0,cz0,cx,cy,cz,x,y,z,vx,vy,
     +                   vz,xrep,yrep,zrep,N_rep,nfl,ras,nfr,
     +                         L_cut,exc,amcheck_dis,goodc,mv,Mend_N,
     +                         Mend,Mend_C,nfl6,ras6,nfl5,ras5,nfl4,
     +                         ras4,nfl3,ras3,nfl2,ras2)

      print *,"After read_initial2 nfl2=",nfl2

      do i=1,100
c     number of accepted replica swaps for replica "i"
         bNSa(i)=0              
c     total number of replica swaps for replica "i"
         bNSt(i)=0
c     number of accepted Monte Carlo moves for move-type "i"
         bNa(i)=0    
c     total number of Monte Carlo moves for move-type "i"
         bNt(i)=0
c     number of accepted Monte Carlo moves for replica "i"
         bNNa(i)=0              !acceptance for different temperature.
c     total number of Monte Carlo moves for replica "i"
         bNNt(i)=0
c     number of times we calculate the total energy of replica "i"
         N_sum(i)=0
c     current average energy of replica "i", updated every cycle
         energ_sum(i)=0
c     current average of the square of the energy of replica "i",
c     updated every cycle
         energ_sum2(i)=0       
      enddo


c     current number of snapshots stored in the trajectories
      i_tr=0  
c     absolute energy minimum observed
      E_min=10000

ccccc The main cycle start from here
      do 1111 icycle=1,ncycle
c     iterate for all the replicas
         do 2222 itemp=1,N_rep  
c     temperature of the replica
            atemp=aT_rep(itemp)	
c     temperature for the square energy transformation
            aTs=aTs_rep(itemp)
c     temperature for the arcSh energy transformation
            aTTs=aTTs_rep(itemp)
c     load replica coordinates into (x(i),y(i),z(i)) and peptide bond
c     indexes into ica(i)
CM            call set_current
            call set_current(x,y,z,xrep,yrep,zrep,ica,vector,Lch,itemp,
     +                     Lch1)
       
           print *,'before initial_move',exc,dist_CA_cut

           call initial_move(icnto,sumcto,Lch,nop,noa,nom,eprofo,seq,
     +                        envir,icnt,sumct,nopp,noaa,nomm,
     +                        codevsum,conew,didevsum,dinew,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2,
     +                hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb,vvv,exc,acut_cc,acut_bb,acut_vv,bcut_cc,
     +                bcut_bb,bcut_vv,dord,en2,en3,contt,
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,switch,
     +             Cr2a,Cr2b,eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,AA,
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,colim,codev,
     +                didev,jbin)

CM      do ijk = 0,500
CM        print *,'After initial_move jbin=',ijk,jbin(ijk)
CM      enddo
CM      CALL EXIT(-1)
            print *,'Before loop 3333=',phot,Lch
ccc
CM            do 3333 iphot=1,phot !N_swap, iterate at fixed temperature
            do iphot=1,phot !N_swap, iterate at fixed temperature

            print *,'Before loop 4444=',no,fff
            print *,'BH values=',bh2,bh3s,bh3d,bh4s,bh4d,bh5s,
     +                bh5d,bh6,bhendn
CM               do 4444 i_lch=1,Lch
               do i_lch=1,Lch
                  fff=aranzy()
CM                  print *,'NO=',no
CM                  print *,'4444 loop=',i_lch,fff
                  if(fff.le.bh2)then
                       print *,'fff le bh2'
CM                     call move2
                       print *,'before move2 nfl2=',nfl2,itemp
                       call move2(nfl2,ras2,Np2,ica,v21,goodc,v22,x,y,z,
     +                      vx,vy,vz,Lch,nop,nopp,nom,nomm,noa,noaa,
     +                      afsn,afs,seq,envir,dord,en1,eprofo,
     +                      bNa,bNNa,itemp,icnto,icnt,sumcto,sumct,
     +                      codevsum,codev,didevsum,didev,energ,bNt,
     +                      bNNt,jbin,es2,er1,er3,er4,es3,es3a,es3b,
     +                      es3c,es4,es5,es6,er5,er6,er7,conew,colim,
     +                      coold,dinew,diold,eprofn,ehb1,ehb1a,ehb1b,
     +                      ehb1c,ehb2,ehb3,ehb4,ehb5a,enb5b,mv,gx,gy,
     +                      gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +                      egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +                      apar,alp,arlp,app,ala,arla,apa,alm,arlm,apm,
     +                      sec,concut2,acut_cc,acut_bb,acut_vv,bcut_cc,
     +                      bcut_bb,bcut_vv,en2,en3,contt,Cr2a,Cr2b,eh1,
     +                      eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +                      acut_hh,bcut_hh,EHBIJ,hh,ee,area,mdis,n_dis,
     +                      kdis,dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                      McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                      frga,frgb,ESHORT2,ESHORT3,ESHORT4,ESHORT4a,
     +                      ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c,ESHORT6,
     +                      ESHORT7,ESHORT8,ESHORT9,ESHORT10,ESHORT11,
     +                      cex,cey,cez,acrit,distL,aweig,dist_CA_cut,
     +                      dis,switch,AA,exc,vvv,atemp,ichos,aTs) 
                  elseif(fff.le.bh3s)then
                       print *,'fff le bh3s'
CM                     call move3s
                       call move3s(nfl3,ras3,Np2,ica,v21,goodc,v22,
     +                  x,y,z,vx,vy,vz,nop,nopp,nom,nomm,noa,noaa,
     +                  afsn,afs,seq,envir,eprofo,bNa,bNNa,
     +                  icnto,icnt,sumcto,sumct,Lch,eprofn,codevsum,
     +                  codev,didevsum,didev,bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
                  elseif(fff.le.bh3d)then
                       print *,'fff le bh3d'
CM                     call move3d
                          call move3d(nfl3,ras3,Np2,ica,v21,goodc,v22,
     +                  x,y,z,vx,vy,vz,Lch,nop,nopp,nom,nomm,noa,noaa,
     +                  afsn,afs,seq,envir,energ,bNa,bNNa,icnto,
     +                  icnt,sumct,sumcto,eprofo,eprofn,codevsum,
     +                  codev,didevsum,didev,bNt,bNNt,itemp,
     +                 jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
                  elseif(fff.le.bh4s)then
                       print *,'fff le bh4s'
CM                     call move4s
                          call move4s(nfl4,ras4,Np2,ica,v21,v22,goodc,
     +                  x,y,z,vx,vy,vz,Lch,nop,nopp,nom,nomm,noa,
     +                  noaa,afsn,afs,eprofn,seq,envir,energ,bNa,
     +                  bNNa,icnto,icnt,sumcto,sumct,eprofo,codevsum,
     +                  codev,didevsum,didev,bNt,bNNt,itemp,     
     +                 jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
                  elseif(fff.le.bh4d)then
                       print *,'fff le bh4d'
CM                     call move4d
                        call move4d(nfl4,ras4,Np2,ica,v21,v22,goodc,x,
     +                  y,z,vx,vy,vz,Lch,nop,nopp,nom,nomm,noa,noaa,
     +                  afsn,afs,eprofn,seq,envir,bNa,bNNa,icnto,
     +                  icnt,sumcto,sumct,eprofo,codevsum,codev,
     +                  didevsum,didev,bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
                  elseif(fff.le.bh5s)then
                       print *,'fff le bh5s'
CM                     call move5s
                        call move5s(nfl5,ras5,Np2,ica,v21,goodc,v22,x,y,
     +                  z,vx,vy,vz,Lch,nop,nopp,nom,nomm,noa,noaa,afsn,
     +                  afs,eprofn,envir,bNa,bNNa,icnto,icnt,sumcto,
     +                  sumct,eprofo,codevsum,codev,didevsum,didev,
     +                  bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,seq,AA,atemp,ichos,aTs) 
                  elseif(fff.le.bh5d)then
                       print *,'fff le bh5d'
CM                     call move5d
                          call move5d(nfl5,ras5,Np2,ica,v21,goodc,v22,x,
     +                  y,z,vx,vy,vz,
     +                  Lch,nop,nopp,nom,nomm,noa,noaa,afsn,afs,eprofn,
     +                  seq,envir,bNa,bNNa,icnto,icnt,sumcto,sumct,
     +                  eprofo,codevsum,codev,didevsum,didev,bNt,bNNt,
     +                  itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
                  elseif(fff.le.bh6)then
                       print *,'fff le bh6'
CM                     call move6
                       call move6(nfl6,ras6,vx,vy,vz,ica,vector,goodc,x,
     +                 y,z,Lch,
     +                 nop,nopp,nom,nomm,noa,noaa,afsn,afs,seq,eprofn,
     +                 envir,eprofo,bNa,bNNa,icnto,icnt,sumcto,sumct,
     +                 codevsum,codev,didevsum,didev,bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
                  elseif(fff.le.bhendn)then
                       print *,'fff le bhendn'
CM                     call move_n_end
                          call move_n_end(Mend_N,x,y,z,ica,anvec,goodc,
     +                      vx,vy,vz,Lch,
     +                      nop,nopp,nom,nomm,noa,noaa,afsn,afs,seq,
     +                      eprofn,envir,bNa,bNNa,icnto,icnt,sumcto,
     +                      sumct,eprofo,codevsum,codev,didevsum,didev,
     +                      bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
                  else
                       print *,'fff le move_c_end'
CM                     call move_c_end
                          call move_c_end(Mend_C,Lch,x,y,z,ica,anvec,
     +                      goodc,vx,
     +                      vy,vz,nop,nopp,nom,nomm,noa,noaa,afsn,afs,
     +                      eprofn,seq,envir,bNa,bNNa,icnto,icnt,
     +                      sumcto,sumct,eprofo,codevsum,codev,didevsum,
     +                      didev,
     +                      bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,Lch2,exc,vvv,atemp,ichos,aTs) 
                  endif
CM 4444          continue
           enddo
CM 3333       continue
       enddo
ccc   
ccccccccccrecord energy and (x,y,z) cccccccccccccccccccc
CM            E_rep(itemp)=energy_tot() !whole energy
CM              E_rep(itemp)=energy_tot(icnto,sumcto,Lch,nop,noa,nom,
CM     +                                seq,envir,
CM     +                    acorder,contt,en1,en2,en3,icnt,
CM     +                    sumct,codevsum,conew,didevsum,dinew,
CM     +                    nopp,noaa,nomm,
CM     +                    eprofo,icnto,sumcto) !whole energy

             E_rep(itemp)=energy_tot(icnto,sumcto,Lch,nop,noa,nom,seq,
     +                    envir,
     +                    acorder,contt,en1,en2,en3,icnt,
     +                    sumct,codevsum,conew,didevsum,dinew,nopp,
     +                    noaa,nomm,
     +                    eprofo,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin, 
     +                hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,colim,codev,
     +                didev)

            if(E_rep(itemp).lt.E_min) E_min=E_rep(itemp)
            do i=1,Lch
               xrep(i,itemp)=x(i)
               yrep(i,itemp)=y(i)
               zrep(i,itemp)=z(i)
            enddo
 2222    continue

cccccccccccccccccc print out 'swep.d' cccccccccccccccccccccccccccc
c         write(91,91)icycle,(E_rep(i),i=1,20)       !$$
c         write(92,91)icycle,(E_rep(i),i=21,N_rep)   !$$
c 91      format(i10,21f9.1)

ccccccccccccccccc<RMSD>, <E> cccccccccccccccccccccccccc
         do i=1,N_rep            
            energ_sum(i)=energ_sum(i)+E_rep(i)
            energ_sum2(i)=energ_sum2(i)+E_rep(i)*E_rep(i)
            N_sum(i)=N_sum(i)+1
         enddo

ccccccccccccccccccccc snapshots of E(1), E(N_rep) ccccccccccccc
         if(icycle.eq.icycle/1*1)then
            i_tr=i_tr+1
            do k=1,n_repf
                 print *,'To file=',Lch,E_rep(k),i_tr,icycle
CM               write(30+k,'i8,1x,f10.1,2i8')Lch,E_rep(k),i_tr,icycle
CM               write(30+k,*)Lch,E_rep(k),i_tr,icycle
               ijk=30+k
               print *,'ijk=',ijk
               write(ijk,91)Lch,E_rep(k),i_tr,icycle
 91            format(i8,1x,f10.1,1x,i8,1x,i8)
               do i=1,Lch
                  abx=xrep(i,k)*0.87
                  aby=yrep(i,k)*0.87
                  abz=zrep(i,k)*0.87
CM                  write(30+k,'f10.3,1x,f10.3,1x,f10.3')abx,aby,abz
CM                  write(30+k,*)abx,aby,abz
                  write(30+k,92)abx,aby,abz
 92               format(f10.3,1x,f10.3,1x,f10.3)
               enddo
            enddo
         endif

CM         call count_restrains   !count number of satisfied restraints
         call count_restrains(Lch,x,y,z,xrep,yrep,zrep,ica,vector,gx,
     +                           gy,gz,seq,Mcom,Kcom,Mdis,kdis,t_dist,
     +                           s_dist,dev,McomCA,KcomCA,t_combCA,
     +                           dist_CA_cut,s_combCA,
     +                           t_comb,concut2,s_comb,dist,Lch1)

ccccccccccccccccc swap replicas cccccccccccccccccccccccccccccccc
         if(icycle.eq.icycle/2*2)then
            do i=1,N_rep-1,2    !swap odd replicas
CM               call swap(i,i+1)
               call swap(i,i+1,bNSa,Lch,xrep,yrep,zrep,aT_rep,E_rep,
     +                  bNSt)
            enddo
         else
            do i=2,N_rep-1,2
CM               call swap(i,i+1) !swap even replicas
               call swap(i,i+1,bNSa,Lch,xrep,yrep,zrep,aT_rep,E_rep,
     +                   bNSt) !swap even replicas
            enddo
         endif
 1111 continue


CM      call test_neighbor        !check all the neighboring residues
      call test_neighbor(Lch,ica,goodc,vx,vy,vz)

CM      call test_overlap         !test the overlap of C_a and C_b
      call test_overlap(Lch,exc,mv,x,y,z,ex,ey,ez)

CM      energy_tot_tmp=energy_tot()

       energy_tot_tmp=energy_tot(icnto,sumcto,Lch,nop,noa,nom,seq,envir,
     +                    acorder,contt,en1,en2,en3,icnt,
     +                    sumct,codevsum,conew,didevsum,dinew,nopp,
     +                    noaa,nomm,
     +                    eprofo,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin, 
     +                hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,colim,codev,
     +                didev)


      write(20,*)'E_final=',energy_tot_tmp

      write(20,*)
      write(20,*)'<s_comb>=',s_comb/float(ncycle)
      write(20,*)'<t_comb>=',t_comb/float(ncycle)
      write(20,*)'s_comb/t_comb=',float(s_comb)/(t_comb+0.001)
      write(20,*)
      write(20,*)'<s_dist>=',s_dist/float(ncycle)
      write(20,*)'<t_dist>=',t_dist/float(ncycle)
      write(20,*)'s_dist/t_dist=',float(s_dist)/(t_dist+0.001)
      write(20,*)
      write(20,*)'<s_combCA>=',s_combCA/float(ncycle)
      write(20,*)'<t_combCA>=',t_combCA/float(ncycle)
      write(20,*)'s_combCA/t_combCA=',float(s_combCA)/(t_combCA+0.001)
      write(20,*)

      write(*,*)'<s_comb>=',s_comb/float(ncycle)
      write(*,*)'<t_comb>=',t_comb/float(ncycle)
      write(*,*)'s_comb/t_comb=',float(s_comb)/(t_comb+0.001)
      write(*,*)'<s_dist>=',s_dist/float(ncycle)
      write(*,*)'<t_dist>=',t_dist/float(ncycle)
      write(*,*)'s_dist/t_dist=',float(s_dist)/(t_dist+0.001)
      write(*,*)'<s_combCA>=',s_combCA/float(ncycle)
      write(*,*)'<t_combCA>=',t_combCA/float(ncycle)
      write(*,*)'s_combCA/t_combCA=',float(s_combCA)/(t_combCA+0.001)

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
         energ_sum(i)=energ_sum(i)/N_sum(i)
         energ_sum2(i)=energ_sum2(i)/N_sum(i)
         if(bNNt(i).gt.1)then
            cheat=(energ_sum2(i)-energ_sum(i)**2)/(aT_rep(i)**2)
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i),bNNa(i)/bNNt(i),aTs_rep(i),aTTs_rep(i)
         else
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i)
         endif
      enddo
 5006 format(I4,f7.2,f8.1,f15.3,2f12.1,f11.6,6f8.1)
      write(20,*)'E_min=',E_min


      STOP
      END
ccccccc=======================================================cccccccccc
cc          The main program ended!
ccccccc=======================================================cccccccccc










ccccc set common used parammeters
      subroutine set_common(random,Lch,Lch1,Lch2,anvec,contt,maxdis2,
     +                      ichos,vvv,mk_pair3,
     +                      t_comb,t_dist,t_combCA,s_comb,s_dist,
     +                      s_combCA,atemp1,atemp2,N_rep,mname,
     +                      izscore,er1,er3,er4,er5,eh1c,er6,er7,chuan,
     +                      ssp,i_bigbond,fract_pair1,fract_pair3,teco)

      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)
      character protein*10
      integer maxdis2(ndim)
      dimension vvv(ndim,ndim)
      character*6 mname(100)
      character type*10

      common/seed1/no

CM      common/lengths/Lch,Lch1,Lch2
CM      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
CM      common/maxdis2/maxdis2(ndim)
CM      common/arandom/  aarand,abrand,acrand,adrand
CM      common/seed1/no
CM      common/nswap/bNSa(100),bNSt(100)
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      common/commonuse1/random,ncycle
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/sw1/aT_rep(nrep),E_rep(nrep)
CM      character*6 mname
CM      character*80 line
CM      common/movename/mname(100)
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      common/distres/er4,es3c
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM      common/zscore/izscore
CM      common/excluded/vvv(ndim,ndim)
CM      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
CM      common/pair1/eh2,eh1b,eh1c
CM      common/weight/chuan
CM      character type*10
CM      common/bigbond/i_bigbond,teco
CM      common/ssp/ssp
CM      common/fractpair/fract_pair1,fract_pair3
CM      common/pair33/i_pair3,mk_pair3
CM      common/ichos/ichos
CM
ccccc set the random generator
      no=random
      print*,'NO=',no

      if(no.gt.0)no=-no

CM Why is this needed?
      firstrandom=aranzy()
CM      firstrandom=aranzy(no,ix1,ix2,ix3,r)
      print*,'firstrandom',firstrandom
                  print *,'NO=',no

ccccc read length of chain and name of protein
      read(21,*)
      read(21,*)Lch
      read(21,*)protein
      write(*,*)Lch,protein
      Lch1=Lch-1
      Lch2=Lch-2
      print*,'Lch 1 2',Lch,protein,Lch1,Lch2

ccccc float version of nvec (=312)
      anvec=nvec-0.00001
      print*,'anvec',anvec

ccccc contt is expected number of contacts   1.5*Lch
      contt=1.5*float(Lch)
      print*,'contt',contt

ccccc maxdis2(i) maximum distance of walk in i steps
      do i=1,Lch
         maxdis2(i)=20*i*i
CM      print*,'maxdis2',i,maxdis2(i)
      enddo

      write(20,*)'Target:',protein
      write(20,*)'Length:',Lch

ccccc ichos determines which "Boltzmann" contact to use in the
c     metropolis algorithm
      ichos=1

ccccc every pair should be checked
      do i=1,Lch
         do j=1,Lch
            vvv(i,j)=1
CM      print*,'vvv',i,j,vvv(i,j)
         enddo
      enddo

ccccc mk_pair=1 means use 'pair3.dat' and 'pair1.dat' (-1 for do not
c     use)
      mk_pair3=1

      t_comb=0
      t_dist=0
      t_combCA=0
      s_comb=0
      s_dist=0
      s_combCA=0

ccccc atemp1=minimal temperature, atemp2=maximal temp.
c     [80,130] is the standard:
c     small corrections for short proteins
      if(Lch.lt.80)then
         atemp1=atemp1*0.97
         atemp2=atemp2*0.91
         if(Lch.lt.55)then
            atemp1=atemp1*0.97
            atemp2=atemp2*0.91
         endif
      endif

c     small corrections for long proteins
      if(Lch.gt.130)then
         atemp1=atemp1*1.05
         atemp2=atemp2*1.1
         if(Lch.gt.165)then
            atemp1=atemp1*1.05
            atemp2=atemp2*1.1
            if(Lch.gt.200)then
               atemp1=atemp1*1.05
               atemp2=atemp2*1.1
            endif
         endif
      endif

CM      print*,'atemp1 2',atemp1,atemp2

ccccc Number of replicas (range from 40 to 80)
      if(Lch.gt.165)then
         N_rep=N_rep+10
      endif
      if(Lch.gt.240)then
         N_rep=N_rep+10
      endif
      if(Lch.gt.300)then
         N_rep=N_rep+10
      endif
      if(Lch.gt.400)then
         N_rep=N_rep+10
      endif
      if(N_rep.gt.80)N_rep=80
      
CM       print*,'N_rep',N_rep

ccccc movement name
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
      mname(19)='rot_N'
      mname(20)='rot_M'
      mname(21)='rot_C'
      mname(22)='trot_N'
      mname(23)='trot_M'
      mname(24)='trot_C'
      mname(25)='defo_N'
      mname(26)='defo_M'
      mname(27)='defo_C'

ccccc restraints weights according to zscore
c     unit=24 <--> chain.dat
      rewind(24)
      read(24,*)n_thr,type
      write(*,*)n_thr,type
      print*,'ssp',ssp

      if(type.eq.'easy')then
         izscore=1
      print*,'izscore',izscore

c     er1 for dist.dat
         er1=chuan*3.6
c     er3 for comb.dat
         er3=chuan*0.765
c     er4 for comb.dat of deviation
         er4=chuan*0.45
c     er5 for combCA.dat
         er5=chuan*2.7
c     eh1c for par.dat
         eh1c=chuan*1.8
c     er6 for distL.dat
         er6=chuan*0.45
c     er7 for RMSD
         er7=chuan*500
         if(ssp.eq.1)er7=chuan*100 !for RMSD
c     decide fragment base on distance
         i_bigbond=3
         fract_pair1=0.4
         fract_pair3=0.3
      elseif(type.eq.'medm')then
         izscore=2
         er1=chuan*4.05 
         er3=chuan*0.81 
         er4=chuan*0.405
         er5=chuan*1.08 
         eh1c=chuan*1.0 
         er6=chuan*1.0  
         er7=chuan*100  
         if(ssp.eq.1)er7=chuan*10
c     decide fragment base on +-2
         i_bigbond=1
         fract_pair1=0.7
         fract_pair3=0.3
c     the target is hard
      else              
         izscore=3
         er1=chuan*2.7  
         er3=chuan*0.4  
         er4=chuan*0.27 
         er5=chuan*0.4  
         eh1c=chuan*1.5 
         er6=chuan*0.5  
         er7=chuan*5    
c     decide fragment base on distance
         i_bigbond=3    
         fract_pair1=0.3
         fract_pair3=0.7
      endif
CM      print*,'er s',er1,er3,er4,er5,eh1c,er6,er7,i_bigbond,
CM     +        fract_pair1,fract_pair3

c     teco=1, template from teco
      if(teco.eq.1) i_bigbond=1
CM      print*,'techo i_bogbond',teco,i_bigbond

      return
      end
ccccc common parameters finished



ccccc read sequence
c     Only seq(i) is useful
      subroutine read_seq (Lch,sec,seq,sequ)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      character*3 aa(-1:20), name,sequ(ndim)
      integer sec(ndim),seq(ndim)
      integer icg(ndim)


CM      common/seqe/seq(ndim),sec(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/icgg/ icg(ndim), EH6  
CM      common/aminoacid/sequ(ndim)

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
c                -1    0     1     2     3     4     5     6
     $     'PRO','MET','ASP','ASN','LEU',
c            7     8     9    10    11
     $     'LYS','GLU','GLN','ARG',
c           12    13    14    15
     $     'HIS','PHE','TYR','TRP','CYX'/
c           16    17    18    19    20

      print*,'Lch',Lch
      do 121 i=1,Lch
c     unit=14 <--> seq.dat
         read(14,707) k,name,sec(i),tmp
CM          write(*,707) k,name,sec(i),tmp
         do j=0,19
CM            print*,'aa',j,aa(j)
            if(name.eq.aa(j)) then
               seq(i)=j
c     icg(i) is used nowhere in the program, so it's useless
               icg(i)=0
               sequ(i)=name
               if(name.eq.'ASP'.or.name.eq.'GLU') icg(i)=-1
               if(name.eq.'LYS'.or.name.eq.'ARG') icg(i)= 1	
CM       print*,'seq sequ icg',i,j,seq(i),sequ(i),icg(i)

               go to 121
            endif
         enddo
         seq(i)=0
         icg(i)=0
         sequ(i)='GLY'
 121  continue
 707  format(i5,3x,a3,2i5)
      close(14)
      return
      end
ccccc read sequence finished



ccccc read centrosymmetric potential
c     eonekd(A) controls centrosymmetric potential of C_g.
      subroutine read_centro(acrit,Lch,eonekd,eonehw)
      implicit integer(i-z)
      parameter(nvec=312)
      real eonekd(0:19),eonehw(0:19),acrit

CM      character*3 NAME
CM      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
CM      common/lengths/Lch,Lch1,Lch2
CM      common/hopp/eonehw(0:19)

CM Why the data needed?

c     read hydrophobic potential for Sg, positive for hydrophobic
c     residue
CM      data eonekd /-0.4, 1.8, -0.8, 2.5, 4.2, -0.7, 4.5, -1.6, 1.9,
CM     $     -3.5, -3.5, 3.8, -3.9, -3.5, -3.5, -4.5, -3.2, 2.8, -1.3,
CM     $     -0.9/

c     read hydrophilic potential for Sg, positive for hydrophilic
c     residue
CM      data eonehw /0.0, -0.5, 0.3, -1.0, -1.5, -0.4, -1.8, 0.0,
CM     $     -1.3, 3.0, 0.2, -1.8, 3.0, 3.0, 0.2, 3.0,-0.5, -2.5,
CM     $     -2.3, -3.4/

      eonekd(0) = -0.4
      eonekd(1) =  1.8 
      eonekd(2) = -0.8 
      eonekd(3) = 2.5 
      eonekd(4) = 4.2 
      eonekd(5) = -0.7 
      eonekd(6) = 4.5 
      eonekd(7) = -1.6 
      eonekd(8) = 1.9
      eonekd(9) = -3.5 
      eonekd(10) = -3.5 
      eonekd(11) = 3.8 
      eonekd(12) = -3.9 
      eonekd(13) = -3.5 
      eonekd(14) = -3.5 
      eonekd(15) = -4.5 
      eonekd(16) = -3.2 
      eonekd(17) = 2.8 
      eonekd(18) = -1.3
      eonekd(19) = -0.9

      eonehw(0) = 0.0
      eonehw(1) = -0.5
      eonehw(2) = 0.3 
      eonehw(3) = -1.0
      eonehw(4) = -1.5
      eonehw(5) = -0.4 
      eonehw(6) = -1.8
      eonehw(7) = 0.0
      eonehw(8) = -1.3
      eonehw(9) = 3.0 
      eonehw(10) = 0.2
      eonehw(11) = -1.8 
      eonehw(12) = 3.0 
      eonehw(13) = 3.0 
      eonehw(14) = 0.2
      eonehw(15) = 3.0
      eonehw(16) = -0.5
      eonehw(17) = -2.5 
      eonehw(18) = -2.3
      eonehw(19) = -3.4


      do i = 0,19
CM         print*,'eonekd eonehw',eonekd(i),eonehw(i)
      enddo
c     expected gyration radius: gyrat-radius~2.2*l^0.38
c     Definition of gyration-radius: acrit=sqrt(<(r-r0)^2>)
      acrit=2.2*exp(0.38*alog(float(Lch)))/0.87
      print*,'acrit',acrit
      return
      end
cicccc     read centersymmetric potential finished



ccccc read environment
c         G  A  V  L  I  S  T  C  M  P  D  N  E  Q  K  R  H  F  Y  W
c (0,0,0) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (0,0,1) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (0,0,2) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c   ...                    ...
c (4,4,3) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (4,4,4) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
ccccc 
      subroutine read_profile(envir)
      implicit integer(i-z)
      real envir(0:15,0:15,0:15,0:19,4)
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM
      print*,'In read_profile'
c     initialize envir(im,ip,ia,i,kkk) to maximum value
      do kkk=1,4
         do i=0,19
            do im=0,15
               do ip=0,15
                  do ia=0,15
                     envir(im,ip,ia,i,kkk)=2.0
CM                     print*,'envir',envir(im,ip,ia,i,kkk)
                  end do
               end do
            end do
         end do
      enddo

c     PROFILE3 environment potential, describes the probability that a
c     particular residue type makes 'im' orthogonal contacts, 'ia'
c     antiparallel contacts and 'ip' parallel contacts with the
c     surrounding amino acids. Number of contacts are taken modulo 2,
c     i.e. 0-1 contact, 2-3,...
      
c     i: residue type
      do i=0,19                 !from column
c     unit=2 <--> profile3.comm
c     we don't care first line
         read(2,*)
c     im: number of orthogonal contacts
         do im=0,8
c     ia: number of antiparallel contacts
            do ia=0,8
c     ip: number of parallel contacts
c     last index is ALLWAYS equal to 3.
               read(2,*)(envir(ia,im,ip,i,3),ip=0,8)
CM               write(*,*)(envir(ia,im,ip,i,3),ip=0,8)
            end do
            read(2,*)
         end do
         read(2,*)
      enddo
      return
      end
ccccc read profile finished


ccccc read 1-3 potential
c     interaction between 1'th CA and 3'th CA
      subroutine read_E13(csr,seq,Lch2)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      real csr(ndim,2)
      dimension seq(ndim)

      dimension csre(0:19,0:19,2)

CM      common/short2/codevsum,didevsum,csr(ndim,2)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

      print*,'In read_E13'
c     R13 potential - two bins only (helical and expanded)
c     r2<48, E=csre(i,j,1); r2>48, E=csre(i,j,2)
      do i=0,19
         do j=0,19
c     unit=5 <--> r13.comm
            read(5,*)
            read(5,*) (csre(i,j,k),k=1,2)
CM            write(*,*) (csre(i,j,k),k=1,2)
         enddo
      enddo

      do i=1,Lch2
         do k=1,2
            csr(i,k)=2.0*csre(seq(i),seq(i+2),k)
CM            print*,'csr',i,k,csr(i,k)
         enddo
      enddo
      return
      end
ccccc read E13 finished


ccccc read 1-4 potential
c     interaction between 1'th CA and 4'th CA
c     the aim is to obtain IBIN(r14), asr(i,IBIN)
      subroutine read_E14(Lch,asr,sec,seq,ibin)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      real asrr(0:19,0:19,-12:12)
      DIMENSION asrh(0:19,0:19,-12:12),asre(0:19,0:19,-12:12)
      integer ibin(-300:300),sec(ndim),seq(ndim)
      real asr(ndim,-12:12)

CM      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/forpreparemove4/asrr(0:19,0:19,-12:12)
CM      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

      print*,'read_E14'
ccccc read asrr,asrh,asre
      do i=0,19                 !asrr(Ai,Bi,dis) from 'r14.comm'
         do j=0,19
c     unit=12 <--> r14.comm
            read(12,*)
CM            write(*,*)
c     asrr(i,j,k) is the potential between Calpha of amino acid of type
c     "i" and Calpha of amino acid of type "j" separated by three amino
c     aids along the sequence, and separated by a combination of
c     distance and chirality included in "k". k takes 25 different
c     values (probably negative values indicate negative chirality).
            read(12,*) (asrr(i,j,k),k=-12,-5)
CM            write(*,*) (asrr(i,j,k),k=-12,-5)
            read(12,*) (asrr(i,j,k),k=-4,3)
CM            write(*,*) (asrr(i,j,k),k=-4,3)
c     k=4 is not allowed
            read(12,*) (asrr(i,j,k),k=5,12)
CM            write(*,*) (asrr(i,j,k),k=5,12)
c     deal with k=1 and rewrite k=3,2,1
            do k=4,1,-1
               asrr(i,j,k)=asrr(i,j,k-1)
CM               print*,'asrr',i,j,k,asrr(i,j,k)
            enddo
         enddo
      enddo

c     asrh(i,j,k) same as asrr but tailored to amino acid with
c     helical propensity
      do i=0,19
         do j=0,19
c     unit=7 r14h.comm
            read(7,*)
            read(7,*) (asrh(i,j,k),k=-12,-5)
CM            write(*,*) (asrh(i,j,k),k=-12,-5)
            read(7,*) (asrh(i,j,k),k=-4,3)
CM            write(*,*) (asrh(i,j,k),k=-4,3)
            read(7,*) (asrh(i,j,k),k=5,12)
CM            write(*,*) (asrh(i,j,k),k=5,12)
            do k=4,1,-1
               asrh(i,j,k)=asrh(i,j,k-1)
CM               print*,'asrh',i,j,k,asrh(i,j,k)
            enddo
         enddo
      enddo

c     asre(Ai,Bi,dis) same as asrr but tailored to amino acid with
c     strand propensity
      do i=0,19
         do j=0,19
c     unit=8 r14e.comm
            read(8,*)
            read(8,*) (asre(i,j,k),k=-12,-5)
CM            write(*,*) (asre(i,j,k),k=-12,-5)
            read(8,*) (asre(i,j,k),k=-4,3)
CM            write(*,*) (asre(i,j,k),k=-4,3)
            read(8,*) (asre(i,j,k),k=5,12)
CM            write(*,*) (asre(i,j,k),k=5,12)
            do k=4,1,-1
               asre(i,j,k)=asre(i,j,k-1)
CM               print*,'asre',i,j,k,asre(i,j,k)
            enddo
         enddo
      enddo
ccccc read asrr,asrh,asre finished

ccccc combine asrr,asrh, and asre into asr, which will be used
      do i=1,Lch-3
         do k=-12,12
c     initialize asr(i,k)
            asr(i,k)=asrr(seq(i+1),seq(i+2),k)
CM            print*,'asr seq',i,k,asr(i,k),seq(i+1),seq(i+2)
CM            print*,'sec',i,sec(i),sec(i+1),sec(i+2),sec(i+3)
c     impose helical propensities on the [i,i+3] segment
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
               if(sec(i).eq.2) then
                  asr(i,k)=asrh(seq(i+1),seq(i+2),k)
               endif
            endif
c     impose strand propensities on the [i,i+3] segment
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
               if(sec(i).eq.4) then
c     note 1.5 reinforces more asre() than asrr()
                  asr(i,k)=asre(seq(i+1),seq(i+2),k)*1.5
               endif
            endif
         enddo
      enddo
ccccc asr(i,ibin(r14)) finished

ccccc if i is square if a distance, then ibin(i) is that distance, but
c     in angstroms and only the integer part. This is the bin number
      do i=1,300
         kk=int((sqrt(float(i))*0.87))+1
c     if distance bigger than 12Angstroms, then assign to distance last
c     bin
         if(kk.gt.12) kk=12
         ibin(i) = kk
         ibin(-i)=-kk
CM         print*,'ibin',i,-i,ibin(i),ibin(-i)
      enddo
      ibin(0)=ibin(1)
      return
      end
ccccc read E14 finished


cccc read 1-5 potential
c     interaction between 1'th CA and 5'th CA
c     the aim is to obtain jbin(r15), bsr(i,jbin)
      subroutine read_E15(Lch,bsr,sec,seq,jbin,acops)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
CM      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/seqe/seq(ndim),sec(ndim)
      DIMENSION bsrh(0:19,0:19,16)
      dimension bsre(0:19,0:19,16)
      dimension bsrr(0:19,0:19,16)
      dimension jbin(0:500),sec(ndim),seq(ndim)
      dimension bsr(ndim,16), acops(ndim,16)

CM      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

      print *,'In read_E15'
cccc read bsrr,bsrh,bsre
      do i=0,19
         do j=0,19
c     unit=9 r15.dat
            read(9,*)
            read(9,*) (bsrr(i,j,k),k=1,8)
CM            write(*,*) (bsrr(i,j,k),k=1,8)
            read(9,*) (bsrr(i,j,k),k=9,16)
CM            write(*,*) (bsrr(i,j,k),k=9,16)
         enddo
      enddo

      do i=0,19
         do j=0,19
c     unit=10 r15h.dat
            read(10,*)
c     bsrh(i,j,k) tailored for amino acids with helical propensity
            read(10,*) (bsrh(i,j,k),k=1,8)
CM            write(*,*) (bsrh(i,j,k),k=1,8)
            read(10,*) (bsrh(i,j,k),k=9,16)
CM            write(*,*) (bsrh(i,j,k),k=9,16)
         enddo
      enddo	

      do i=0,19
         do j=0,19
c     unit=9 r15e.dat
            read(11,*)
c     bsrhei,j,k) tailored for amino acids with strand propensity
            read(11,*) (bsre(i,j,k),k=1,8)
CM            write(*,*) (bsre(i,j,k),k=1,8)
            read(11,*) (bsre(i,j,k),k=9,16)
CM            write(*,*) (bsre(i,j,k),k=9,16)
         enddo
      enddo	

c     combine bsrr, bsrh, and bsre into bsr, which will be used
      do i=1,Lch-4
         do k=1,16
c     initialize bsh(i,k)
            bsr(i,k)=bsrr(seq(i+1),seq(i+3),k)
CM            print*,'bsr seq',i,k,bsr(i,k),seq(i+1),seq(i+2)
CM            print*,'sec',i,sec(i),sec(i+1),sec(i+2),sec(i+3)
c     impose helical propensities on the [i,i+3] segment
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
               bsr(i,k)=bsrh(seq(i+1),seq(i+3),k)
            endif
c     impose strand propensities on the [i,i+3] segment
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
               bsr(i,k)=bsre(seq(i+1),seq(i+3),k)*1.5
            endif
            bsr(i,k)=bsrr(seq(i+1),seq(i+3),k)
CM            print *,'bsr',i,k,bsr(i,k)
         enddo
      enddo

ccccc if i is square if a distance, then jbin(i) is that distance, but
c     in angstroms and only the integer part. This is the bin number
      do i=0,500
         kk=int((sqrt(float(i))*0.87))+1
         if(kk.gt.16) kk=16
         jbin(i) = kk
CM          print *,'jbin',i,jbin(i)
      enddo

ccccc acops(i,jbin) to enhance the contacts between gragments for a
c     given position "i" and bin number "k", acops(i,k) is the average
c     of bsr over the [k-1,k+1] set of bins, or zero if the average is
c     negative
      do i=1,Lch-4
         acops(i,1)=( min( bsr(i,1), 0.0) )/2.0
CM         print *,'acops 1',i,acops(i,1)
         do k=2,15
            acops(i,k)=min( 0.0, bsr(i,k-1)+2.0*bsr(i,k)+bsr(i,k+1) )
CM            print *,'acops',i,k,acops(i,k)
         enddo
         acops(i,16)=(min(bsr(i,16),0.0))/2.0
CM         print *,'acops 16',i,acops(i,16)
      enddo
      return
      end
cccc  finished read_E15



ccccc read contact-pair potential
c     potential: app(Ai,Aj), apa(Ai,Aj), apm(Ai,Aj)
c     distance range: [arlp,alp], [arla,ala], [arlm,alm]
c     all the general data in 'quarsi3.comm'.
      subroutine read_quarsi3(arla,arlm,arlp,alp,ala,alm,Lch,seq,
     +                        mk_pair3,i_pair3,fract_pair3,apa,apm,app)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      character*3 NAME
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      real apabla(0:19,0:19),apablp(0:19,0:19),apablm(0:19,0:19)
      real apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      dimension apba(ndim,ndim),apbp(ndim,ndim),apbm(ndim,ndim)
      integer seq(ndim)

CM      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
CM      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
CM      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/pair1/eh2,eh1b,eh1c
CM      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
CM      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap
CM      common/paircut/ash
CM      common/zscore/izscore
CM      common/fractpair/fract_pair1,fract_pair3
CM      common/pair33/i_pair3,mk_pair3
CM      common/weight/chuan
CM

      print *,'In read_quarsi3'
c     Pairwise interactions apablp ... and cut-off parmeters
c     arlp, orientation dependent, pairwise specific, sequence
c     independent

c     unit=3 <--> quasi3.comm
      read(3,*)
      read(3,*)
      do i=0,19
c     apablp(i,j) contact potentials between residue-types "i" and "j"
c     when doing a parallel contact
         read(3,725) name, (apablp(i,j),j=0,19)
CM         write(*,725) name, (apablp(i,j),j=0,19)
      enddo

c     apablp(i,j) contact potentials between residue-types "i" and "j"
c     when doing an orthgonal contact
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) name, (apablm(i,j),j=0,19)
CM         write(*,725) name, (apablm(i,j),j=0,19)
      enddo

c     apablp(i,j) contact potentials between residue-types "i" and "j"
c     when doing an antiparallel contact
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) name, (apabla(i,j),j=0,19)
CM         write(*,725) name, (apabla(i,j),j=0,19)
      enddo

ccccc read midpoint distances in the potential range
      read(3,*)
      read(3,*)
      do i=0,19
c     max distance for parallel contacts
         read(3,725) name, (arlp(i,j),j=0,19)
CM         write(*,725) name, (arlp(i,j),j=0,19)
      enddo
c     max distance for orthogonal contacts
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) name, (arlm(i,j),j=0,19) !for perpendicular contact
CM         write(*,725) name, (arlm(i,j),j=0,19) !for perpendicular contact
      enddo
c     max distance for antiparallel contacts
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) name, (arla(i,j),j=0,19) !for antiparellel pair
CM         write(*,725) name, (arla(i,j),j=0,19) !for antiparellel pair
      enddo
 725  format(a3,1x,20f5.1)

c     set minimum and maximum distance-discontinuities in the potential
c     for a contact, in lattice units
      ash=0.17
      ash_min=1-ash
      ash_max=1+ash
      do i=0,19
         do j=0,19
c     distance ranges for parallel contacts 
            alp(i,j)=(arlp(i,j)*ash_max/0.87)**2
            arlp(i,j)=(arlp(i,j)*ash_min/0.87)**2
c     distance ranges for antiparallel contacts 
            ala(i,j)=(arla(i,j)*ash_max/0.87)**2
            arla(i,j)=(arla(i,j)*ash_min/0.87)**2
c     distance ranges for orthogonal contacts 
            alm(i,j)=(arlm(i,j)*ash_max/0.87)**2
            arlm(i,j)=(arlm(i,j)*ash_min/0.87)**2
CM            print *,'alp alm etc',alp(i,j),arlp(i,j),ala(i,j),
CM     +               arla(i,j),alm(i,j),arlm(i,j)
         enddo
      enddo
c     E=EH1/2, for r in [0,arlp];
c     E=app-es*fs, for [0,alp];
c     E=0, for r in [alp,00].
ccccc contact interaction, orientation dependent finished


ccccc incorporate contact interaction, orientation dependent from
c     threading to a PDB structural library of a multiple sequence
c     alignment of homologs to the target.
      i_pair3=-1
c     unit=25 <--> pair3.dat
      read(25,*,end=1000)
c     read first contact strength for antiparallel contacts between
c     residues at positions "i" and "j", apba(i,j)
      rewind(25)
      i_pair3=1
      Nline=100
      do i=1,Lch
         read(25,*)
c     clumsy way of reading because the Lch values are broken into
c     different lines
         do i_line=1,Nline
            line_end=min(25*i_line,Lch)
            read(25,*)(apba(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4074
         enddo
 4074    continue
         read(25,*)
      enddo

c     read now contact-strength for orthogonal contacts
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch)
            read(25,*)(apbm(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4075
         enddo
 4075    continue
         read(25,*)
      enddo

c     read now contact-strength for parallel contacts
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
            read(25,*)(apbp(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4076
         enddo
 4076    continue
         read(25,*)
      enddo
      close(25)
 1000 continue

ccccc combine the data in 'quarsi3.comm' and 'pair3.dat' to get
c     contact potential
      do i=1,lch
         ii=seq(i)
         do j=1,lch
            jj=seq(j)
            if(iabs(i-j).lt.5) then
               dd=0.0
            else
c     encourage contact of distant residues.
               dd=0.25
            endif

ccccc apa(),apm(),app() contain orientation dependent contact potential
c     initialize with target independent potentials (quasi3.comm only)
            apa(i,j)=apabla(ii,jj)-dd
            apm(i,j)=apablm(ii,jj)-dd
            app(i,j)=apablp(ii,jj)-dd
c     combine target dependent and target independent (fract_pair3=0.7)
            if(i_pair3.gt.0.and.mk_pair3.gt.0)then
               apa(i,j) = fract_pair3*apba(i,j) + 
     $              (1-fract_pair3)*apabla(ii,jj) - dd
               apm(i,j) = fract_pair3*apbm(i,j) +
     $              (1-fract_pair3)*apablm(ii,jj) - dd
               app(i,j) = fract_pair3*apbp(i,j) +
     $              (1-fract_pair3)*apablp(ii,jj) - dd
            endif
         enddo
      enddo
      return
      end
ccccc read_quarsi3 finished


ccccc read contact-pair potential, orientation independent, from files
c     par.dat and pair1.dat the sequence-dependent contact-pair data in
c     par.dat.
      subroutine read_par(Lch,mk_pair3,i_pair3,fract_pair1,chuan,eh1c,
     +                    apar)
      implicit integer(i-z)
      parameter(ndim=1500)
      dimension apar(ndim,ndim)

CM      common/lengths/Lch,Lch1,Lch2
CM      common/par/apar(ndim,ndim)
CM      common/fractpair/fract_pair1,fract_pair3
CM      common/pair33/i_pair3,mk_pair3
CM      common/weight/chuan
CM      common/pair1/eh2,eh1b,eh1c

c     for par.dat
      dimension apar1(ndim,ndim)
c     for pair1.dat
      dimension apar2(ndim,ndim)
      print *,'In read_par Lch',Lch
      Nline=1000
      do i=1,Lch
c     unit=15 <--> par.dat
         read(15,*)
         do i_line=1,Nline
            line_end=min(10*i_line,Lch)
            read(15,*)(apar1(i,j),j=(i_line-1)*10+1,line_end)
CM            write(*,*)(apar1(i,j),j=(i_line-1)*10+1,line_end)
            if(line_end.ge.Lch) go to 1
         enddo
 1       continue
      enddo

c     read from pair1.dat
c     pair1.dat exists
      if(i_pair3.gt.0)then
         Nline=100
         do i=1,Lch
c     unit=26 pair1.dat
            read(26,*)
            do i_line=1,Nline
               line_end=min(25*i_line,Lch)
               read(26,*)(apar2(i,j),j=(i_line-1)*25+1,line_end)
CM               write(*,*)(apar2(i,j),j=(i_line-1)*25+1,line_end)
               if(line_end.eq.Lch) go to 2
            enddo
 2          continue
            read(26,*)
         enddo
      endif
      
c     find maximum values of apar1(i,j) and apar2(i,j) to normalize
      apar1_max=-100
      apar2_max=-100
      do i=1,Lch
         do j=1,Lch
            if(abs(apar1(i,j)).gt.apar1_max)apar1_max=abs(apar1(i,j))
            if(abs(apar2(i,j)).gt.apar2_max)apar2_max=abs(apar2(i,j))
CM             print *,'apar1 2 max',apar1_max,apar2_max
         enddo
      enddo

c     combine apar1() and apar2() in apar()
      do i=1,Lch
         do j=1,Lch
            apar(i,j)=0
CM            print *,'mk_pair3',mk_pair3,i_pair3
            if(mk_pair3.gt.0)then
c     if we're using info from pair1.dat
               if(i_pair3.gt.0)then
c     fract_pair1=0.3, more weight to apar1() than to apar2()
CM                   print *,'apar fract_pair1',fract_pair1,apar(i,j),
CM     +                     apar2(i,j),chuan
                  apar(i,j) = apar1(i,j) * apar2_max / apar1_max *
     $                 (1-fract_pair1) + apar2(i,j)*fract_pair1
                  if(chuan.lt.0.00001)then
c     modify weight for pairwise potential between side-chains
                     eh1c=1
c     do NOT use info from par.dat
                     apar(i,j)=apar2(i,j)
                  endif
               else
c     do NOT use info from par1.dat
                  apar(i,j)=apar1(i,j)
               endif
            endif
CM            print *,'apar',apar(i,j)
         enddo
      enddo
      return
      end
ccccc pair-wise potential finished


ccccc read cut-off for contact predictions
      subroutine read_concut(izscore,concut,concut2)
      implicit integer(i-z)
      character*3 NAME
      dimension cut(0:19,0:19),cut_dev(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)

CM      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
CM      common/zscore/izscore

      print *,'In read_concut'
c     unit=1 <--> contact.comm
      rewind(1)

c     cut(i,j) are midpoint distances in the potential for contacts
c     between two side-chains for contact restrain potential
      read(1,*)
      do i=0,19
         read(1,*)name,(cut(i,j),j=0,19)
CM         write(*,*)name,(cut(i,j),j=0,19)
      enddo

c     cut_dev(i,j) are deviation for contact between two side-chains
c     for contact restrain potential
      read(1,*)
      do i=0,19
         read(1,*)name,(cut_dev(i,j),j=0,19)
CM         write(*,*)name,(cut_dev(i,j),j=0,19)
      enddo

      do i=0,19
         do j=0,19
c     izscore value denotes type of target (1:easy, 2:medium, 3:target)
c     easy target
CM            print *,'izscore',izscore
            if(izscore.eq.1)then
               concut(i,j)=7
c     medium target
            elseif(izscore.eq.2)then
               concut(i,j)=7.5
c     hard target
            else
CM               print *,'cut cut_dev',cut(i,j),cut_dev(i,j)
               concut(i,j)=cut(i,j)+cut_dev(i,j)*2.5
            endif
c     switch to lattice units
            concut(i,j)=concut(i,j)/0.87
c     precompute the square of the cutoff distance
            concut2(i,j)=concut(i,j)**2
CM            print *,'2 concut',concut2(i,j),concut(i,j)
         enddo
      enddo
      return
      end
ccccc read_concut finished


ccccc read contact restrains from 'comb.dat'
c     aim: Mcom(i)--->number of contacts on i-residue
c          Kcom(i,ii): ii's contact of i-residue is j=kres
c     contact restrain is only on C_g.
      subroutine read_contactrestrain(izscore,aweig,Mcom,Kcom,Lch,colim,
     +                                Ncom)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      dimension Mcom(ndim),Kcom(ndim,100)
      dimension aweig(4000,4000)

CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
CM      common/lengths/Lch,Lch1,Lch2
CM      common/distres/er4,es3c
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM      common/resnumber/Ncom,Ndis,accur
CM      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
CM      common/freg/aweig(4000,4000)
CM      common/zscore/izscore

      dimension r1(4000),r2(4000)

      print *,'In read_contactrestrain',Lch,izscore

c     reads Side group - side group contacts 
c     (from NMR or therading predictions or clusters)

     
c     easy target
      if(izscore.eq.1)then
ccccc cut_min is the minimun value of the confidence for a
c     contact-restrain between two SG's. cut_min has different values
c     depending on the difficulty of the target. 
         cut_min=0.2
ccccc cut0 is the "zero" for the confidence of a predicted contactd,
c     from which derive the weight of the contact
         cut0=0.6
c     medium target
      elseif(izscore.eq.2)then
         cut_min=0.1
         cut0=0.5
      else
c     hard target
         cut_min=0.1
         cut0=0.4
      endif
       print *,'cut0 cut_min',cut0,cut_min

ccccc pool all the restraints into r1(i),r2(i)
 11   continue
c     unit=16 <--> comb.dat
      read(16,*)ntmp
      i_c=0
      do i=1,ntmp
         read(16,*)i1,i2,conf
CM         write(*,*)i1,i2,conf
         if(conf.gt.1) conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.4000)then
               write(*,*)'Too many COMB restraints!!!!!!!!'
c     increase the minimun confidence for a contact restraint and scan
c     again all contacts
               cut_min=cut_min*1.1
               rewind(16)
               goto 11
            endif
c     r1(),r2() contain the positions of the contacting residues
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweig(i1,i2)=1+abs(conf-cut0)*4
            else
               aweig(i1,i2)=1-abs(conf-cut0)*2
            endif
            aweig(i2,i1)=aweig(i1,i2)
CM            print *,'aweig',aweig(i2,i1),aweig(i1,i2)
         endif
      enddo
c     Ncom is the number of read predicted restraints
      Ncom=i_c

ccccc   map r1,2(i) into Mcom(i),Kcom(i,Mcom(i))
      do i=1,Lch
c     Mcom(i) is the number of predicted contacts that SG at position
c     "i" makes with the rest of the SG's
         Mcom(i)=0
         do j=1,Ncom
            if(r1(j).eq.i)then
               Mcom(i)=Mcom(i)+1
c     Kcom(i,j) is the position of one of the SG's that contact with SG
c     "i". "j" is a running index from 1 to Mcom(i)
               Kcom(i,Mcom(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               Mcom(i)=Mcom(i)+1
               Kcom(i,Mcom(i))=r1(j)
            endif
CM            print *,'Kcom',i,Kcom(i,Mcom(i))
         enddo
      enddo

c     the larger 'colim' is, the weaker the contact restrain strength is      
      colim=1.5*Ncom
      print *,'colim',colim

ccccc output restraints
      write(20,*)'Number of restraints:',Ncom
      write(20,*)'----------- contact restraints ---------------'
      nnc=0
      do i=1,Lch
         nnc=nnc+Mcom(i)
         write(20,12)i,Mcom(i),(Kcom(i,j),j=1,Mcom(i))
CM         write(*,12)i,Mcom(i),(Kcom(i,j),j=1,Mcom(i))
 12      format(i4,'(',i2,'):',20i4)
      enddo
      write(20,*)'Number of contact=',nnc,' Lch=',Lch
      write(20,*)'fc=',float(nnc)/Lch
      write(20,*)
      return
      end
ccccc read_contactrestrain finished


ccccc read distance restrains from 'dist.dat'
      subroutine read_distantrestrain(Lch,Mdis,kdis,dist,dev,Ndis)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      parameter(maxnL=30000)

      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)

CM      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/distres/er4,es3c
CM      common/resnumber/Ncom,Ndis,accur
CM      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b

      dimension r1(maxnL),r2(maxnL),dis(maxnL),deviation(maxnL)

      print *,'In read_distant',Lch
c     unit=17 <--> dist.dat
      read(17,*)Ndis
      do i=1,Ndis
         read(17,*)r1(i),r2(i),nothing,dis(i),deviation(i)
CM         write(*,*)r1(i),r2(i),nothing,dis(i),deviation(i)
         dis(i)=dis(i)/0.87
         deviation(i)=deviation(i)/0.87
         if(deviation(i).lt.0.5)deviation(i)=0.5
      enddo

      do i=1,Lch
c     number of SG's predicted to be at certain distances from SG's "i"
         Mdis(i)=0
         do j=1,Ndis
            if(r1(j).eq.i)then
               Mdis(i)=Mdis(i)+1
c     kdis(i,j) is the position of one of the SG's that are predicted to
c     be at certain distance from SG "i". j runs from 1 to Mdis(i)
               kdis(i,Mdis(i))=r2(j)
c     predicted distance for i<->r2()
               dist(i,Mdis(i))=dis(j)
c     predicted deviation around the predicted distance
               dev(i,Mdis(i))=deviation(j)
            endif
            if(r2(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r1(j)
               dist(i,Mdis(i))=dis(j)
               dev(i,Mdis(i))=deviation(j)
            endif
CM             print *,'dev',i,j,dev(i,Mdis(i)),deviation(j)
         enddo
      enddo

c     background number for derivation
      dilim=0.5*Ndis
      print *,'dlim',dilim
      return
      end
ccccc read contact restrains finished



ccccc read distance restrains from 'distL.dat'
c     each pair (i,j) can have many different distances predictions
      subroutine read_longdistantrestrain(Lch,MdisL,kdisL,distL)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      parameter(maxnL=30000)

      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)

CM      common/lengths/Lch,Lch1,Lch2
CM      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)

      dimension r1(maxnL),r2(maxnL),dis(maxnL),deviation(maxnL)

       print *,'In read_long',Lch
c     set the long distance restrain cut-off to be 30 angstroms
      aldrco=30.0
c     unit=22 <--> distL.dat
      read(22,*) j
      NdisL=0      
      do i=1,j
         read(22,*)m,n,d
c     record restrain if only smaller than distance cut-off
         if(d.le.aldrco)then
            NdisL=NdisL+1
            r1(NdisL)=m
            r2(NdisL)=n
            dis(NdisL)=d/0.87
         endif
      enddo


      do i=1,Lch
c     MdisL(i) is the number of SG's predicted to be at long distances from SG's "i"
        MdisL(i)=0
        do j=1,NdisL
          if(r1(j).eq.i)then
            MdisL(i)=MdisL(i)+1
c     kdisL(i,j) is the position of one of the SG's that are predicted
c     to be at certain long distance from SG "i". j runs from 1 to
c     MdisL(i)
            kdisL(i,MdisL(i))=r2(j)
            distL(i,MdisL(i))=dis(j)
          endif
          if(r2(j).eq.i)then
            MdisL(i)=MdisL(i)+1
            kdisL(i,MdisL(i))=r1(j)
            distL(i,MdisL(i))=dis(j)
          endif
CM          print *,'distL',i,j,distL(i,MdisL(i)),dis(j),MdisL(i)
        enddo
      enddo
      return
      end
ccccc read contact restrains finished

ccccc read contact restrains from 'comb.dat'
c     aim: Mcom(i)--->number of contacts on i-residue
c          Kcom(i,ii): ii's contact of i-residue is j=kres
c     CAcontact restrain is only on CA.
      subroutine read_CAcontact(izscore,dist_CA_cut,aweigCA,Lch,
     +                          McomCA,KcomCA)
      implicit integer(i-z)
      parameter(ndim=1500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)

CM      common/lengths/Lch,Lch1,Lch2
CM      common/zscore/izscore
CM      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
CM      common/CAcontact1/dist_CA_cut
      dimension r1(4000),r2(4000)


       print *,'In read_CAcon',Lch,izscore
c     easy target
      if(izscore.eq.1)then
         dist_CA_cut=(6.5/0.87)**2
ccccc cut_min is the minimun value of the confidence for a
c     contact-restrain between two CA's. cut_min has different values
c     depending on the difficulty of the target. 
         cut_min=0.2
ccccc cut0 is the "zero" for the confidence of a predicted contactd,
c     from which derive the weight of the contact
         cut0=0.5
c     medium target
      elseif(izscore.eq.2)then
         dist_CA_cut=(6.0/0.87)**2
         cut_min=0.1
         cut0=0.4
c     hard target
      else
         dist_CA_cut=(6.5/0.87)**2
         cut_min=0.1
         cut0=0.3
      endif
      print *,'cut0 etc',dist_CA_cut,cut_min,cut0

ccccc pool all the restraints into r1(i),r2(i)
 11   continue
c     unit=18 <--> combCA.dat
      read(18,*)ntmp
      i_c=0
      do i=1,ntmp
         read(18,*)i1,i2,conf
CM         write(*,*)i1,i2,conf
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.4000)then
               write(*,*)'Too many COMBCA restraints!!!!!'
               cut_min=cut_min*1.1
               rewind(18)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweigCA(i1,i2)=1+abs(conf-cut0)*4
            else
               aweigCA(i1,i2)=1-abs(conf-cut0)*2
            endif
            aweigCA(i2,i1)=aweigCA(i1,i2)
CM            print *,'aweigCA',aweigCA(i2,i1),aweigCA(i1,i2)
         endif
      enddo
      NcomCA=i_c

ccccc map r1,2(i) into McomCA(i),KcomCA(i,McomCA(i))
      do i=1,Lch
c     McomCA(i) is the number of predicted contacts that CA at position
c     "i" makes with the rest of the CA's
         McomCA(i)=0
         do j=1,NcomCA
            if(r1(j).eq.i)then
               McomCA(i)=McomCA(i)+1
c     KcomCA(i,j) is the position of one of the CA's that contact with CA
c     "i". "j" is a running index from 1 to McomCA(i)
               KcomCA(i,McomCA(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               McomCA(i)=McomCA(i)+1
               KcomCA(i,McomCA(i))=r1(j)
            endif
CM             print *,'Kcom etc',i,KcomCA(i,McomCA(i)),McomCA(i)
         enddo
      enddo
      return
      end
ccccc read_CAcontact restrains finished


ccccc read solvent exposure propensities
      subroutine read_exp(nana,Lch,area)
      implicit integer(i-z)
      parameter(ndim=1500)
      dimension area(ndim)

CM      common/lengths/Lch,Lch1,Lch2
CM      common/expose/mp(20,ndim),area(ndim)
CM      common/nana1/nana

      dimension mp(20,ndim)

       print *,'In read_exp',Lch,nana
      na=nana
c     unit=27 <--> exp.dat
c     first line is unimportant
      read(27,*)
c     determine area(i), which is a weight for the propensity of the
c     residue to expose its surface (positive values), or bury it
c     (negative values)
      do i=1,Lch
         read(27,*)itmp,(mp(j,i),j=1,na)
CM         write(*,*)itmp,(mp(j,i),j=1,na)
         area(i)=0
         do j=1,na
c     mp(j,i)=-1, bury; mp(j,i)=1, expose
            if(mp(j,i).eq.0)mp(j,i)=-1
            area(i)=area(i)+mp(j,i)
CM            print *,'area',i,j,area(i)
         enddo
      enddo
      close(27)
      return
      end
ccccc read_exp finished

ccccc reset temperature
      subroutine reset_temperature(er1,er3,er4,Ncom,Lch,atemp2)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100) 
      parameter(nvec=312) 

CM      common/lengths/Lch,Lch1,Lch2
CM      common/resnumber/Ncom,Ndis,accur
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      common/distres/er4,es3c
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)

       print *,'In reset_temp',Lch,er1,er3,er4
      if(er1+er3+er4.lt.0.1) return
      
      a_rest=float(Ncom)/(1.3*Lch)
      a_rest=sqrt(a_rest)
      if(a_rest.lt.0.875)a_rest=0.875 !because of 80/70 -> rest/without_rest
      atemp2=atemp2*a_rest
      if(atemp2.gt.115) atemp2=115.0

      print *,'a_rest atemp2',a_rest,atemp2
      return
      end
ccccc reset_temperature

ccccc temperature for different replicas
      subroutine set_temperature(N_rep,aT_rep,atemp1,atemp2,aTs_rep,
     +                           aTs1,aTs2,aTTs_rep,aTTs1,aTTs2)
      implicit integer(i-z)
      parameter(nrep=100)
      dimension aT_rep(nrep),aTs_rep(nrep),aTTs_rep(nrep)

CM      common/commonuse1/random,ncycle
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      common/sw1/aT_rep(nrep),E_rep(nrep)
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      common/thrtem/aT_ann(100)
CM      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
CM      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs

       print *,'In set_temp',N_REP,atemp1,atemp2,aTs1,aTs2,aTTs1,aTTs2
ccccc for normal run
      do i=1,N_REP
c     replica temperature
         aT_rep(i)=atemp1*(atemp2/atemp1)**(float(i-1)/(N_rep-1))
c     modified temperature for square energy transformation
         aTs_rep(i)=aTs1*(aTs2/aTs1)**(float(i-1)/(N_rep-1))
c     modified temperature for arcsh energy transformation
         aTTs_rep(i)=aTTs1*(aTTs2/aTTs1)**(float(i-1)/(N_rep-1))
CM         print *,'Ats etc',aT_rep(i),aTs_rep(i),aTTs_rep(i)
      enddo
      return
      end
ccccc set_temperature finished


ccccc ehbij(i,j) - set-up secondary structure dependent strength of the
c     hyrogen bond network - stronger for helices and beta-sheets
      subroutine set_EHB(Lch,sec,ehbij)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      dimension sec(ndim),ehbij(ndim,ndim)

CM      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/seqe/seq(ndim),sec(ndim)


      print *,'In set_EHB',Lch
      do i=1,Lch
         is=sec(i)
         do j=1,Lch
            js=sec(j)
            ehbij(i,j)=1
            if(iabs(i-j).eq.3.and.is.eq.2.and.js.eq.2)then
c     hydrogen bond for helix enhanced
               ehbij(i,j)=ehbij(i,j)+0.5
            endif
c     if one of the two amino acids is predicted to be in strand, and
c     none is predicted to be in helix state
            if(is.eq.4.or.js.eq.4) then
               if(is*js.ne.8.and.iabs(i-j).gt.4)then
c     hydrogen bond for strand enhanced
                  ehbij(i,j)=ehbij(i,j)+0.5
               endif
            endif
CM            print *,'is js ehb',i,j,is,js,ehbij(i,j)
         enddo
      enddo
      return
      end
ccccc set H_bond finished


ccccc prepare possible vector v=(vx,vy,vz) satisfy |v*v|=[14,25]
c     (3.26A-4.35A) and vector(-5:5,-5:5,-5:5)
      subroutine prepare_vectors(vx,vy,vz,vector)
      implicit integer(i-z)
c     nvec value is known a posteriori. It's the maximum value of nwmax
      parameter(nvec=312)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension n(100)

      print *,'In prepare_vect'
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
c     14 corresponds to (3.25 Angstroms)^2
c     25 corresponds to (4.35 Angstroms)^2
               if(r.ge.14.and.r.le.25) then
c     nwmax is the possible number of vectors joining two CA's
                  nwmax=nwmax+1
c     vx(i),vy(i),vz(i) are the actual coordinates of one vector
                  vx(nwmax)=x
                  vy(nwmax)=y
                  vz(nwmax)=z
c     vector(x,y,z) is the index for vector with coordinates (x,y,z)
                  vector(x,y,z)=nwmax
c     n(r) is number of vectors with a particular square distance "r"
                  n(r)=n(r)+1
                  aaa=aaa+sqrt(float(r)*0.87*0.87)
               endif
CM                print *,'vx vy vz vector',x,y,z,vx(nwmax),vy(nwmax),
CM     +                 vz(nwmax),vector(x,y,z)
            enddo
         enddo
      enddo
c     aaa is the average square of the distance for all possible
c     vectors, in Angstroms, not in lattice units
      aaa=aaa/float(nwmax)
      write(*,*)'n1_all=',nwmax,'  <vr>=',aaa

      return
      end
ccccc prepare_vectors finished


c     define goodc(i,j), angle(i,j), prod(i,j)
c     prod(i,j): v(i)*v(j).
c     angle(i,j): angle of neighbor bonds, i,j--->(1:nvec)
c     goodc(i,j): ture, when angle(i,j) in [60,160]; false, otherwise.
      subroutine prepare_neighbors(m12,u21,vx,vy,vz,angle,goodc,u1,u2)
      implicit integer(i-z)
      parameter(nvec=312)

      logical goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)
      dimension angle(nvec,nvec)

CM      common/logica/goodc
CM      logical goodc(nvec,nvec)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/three/angle(nvec,nvec)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      print *,'prepare_neigh',Lch
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
c     angle between peptide bonds "i" and "j"
            angle(i,j)=acos(cosangle)*180/3.1415926
c     if angle is in [65,165];
CM            print *,'vx vy vz angle,',i,j,vx(i),vy(i),vz(i),angle(i,j)
            if(angle(i,j).gt.65.and.angle(i,j).lt.165)then
c     accept this angle as belonging to protein structures
               goodc(i,j)=.true.
               mmm=mmm+1
               ijx=vx(i)+vx(j)
               ijy=vy(i)+vy(j)
               ijz=vz(i)+vz(j)
               do k=1,nvec
c     we found a peptide bond resulting from the summ of two peptide
c     bonds
                  if(vx(k).eq.ijx.and.vy(k).eq.ijy.and.vz(k).eq.ijz)then
                     kkk=kkk+1
c     u21(i,j) is the index of peptide bond resulting from addition of
c     peptide bonds "i" and "j"
                     u21(i,j)=k
c     m12(k) is the number of ways in which addition of two peptide
c     bonds results in peptide bond "k"
                     m12(k)=m12(k)+1
c     u1(k,l) and u2(k,l) are indexes of two peptide bonds such that
c     their addition results in peptide bond "k". "l" is just a running
c     index from 1 to m12(k)
                     u1(k,m12(k))=i
                     u2(k,m12(k))=j
CM                      print *,'us m12',k,u21(i,j),m12(k),
CM     +                        u1(k,m12(k)),u2(k,m12(k))
                     if(max_m12.lt.m12(k))max_m12=m12(k)
c     since there's only one v_k such that v_k=v_i+v_j, stop the search
c     for v_k
                     goto 10
                  endif
               enddo
 10            continue
            else
c     angle between "i" and "j" is not representative of what's seen in
c     proteins
               goodc(i,j)=.false.
            endif
            nnn=nnn+1
c     go to next peptide bond "j"
         enddo
c     go to next peptide bond "i"
      enddo

CM Manisha Should delete lines below as not required

      n=0
      do i=1,nvec
         r=vx(i)**2+vy(i)**2+vz(i)**2
         if(m12(i).gt.0)n=n+1
      enddo
      print *,'n',n
      return
      end
ccccc prepare_neighbors

ccccc define hydrogen-bond, C_beta, C_group for all possible neighbors
      subroutine prepare_beta(axalf,ayalf,azalf,axbet,aybet,azbet,vx,vy,
     +                        vz,hbx,hby,hbz,cax,cay,caz,angle,gx,gy,gz)
      implicit integer(i-z)
      parameter(nvec=312)

      logical goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension angle(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension axalf(0:19),ayalf(0:19),azalf(0:19)
      dimension axbet(0:19),aybet(0:19),azbet(0:19)

CM      common/logica/goodc
CM      logical goodc(nvec,nvec)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
CM      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
CM      common/three/angle(nvec,nvec)
CM      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
CM      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
CM      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)

      print *,'In prepare_beta'
      do k=0,19
c     unit=4 <--> sidechain.comm for each amino acid, contains preferred
c     coordinates of SG when in alpha (axalf) or beta (axbet)
c     conformation
         read(4,*)
         read(4,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
CM         write(*,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
      enddo
      CLOSE(4)

ccccc define hydrogen-bond, C_beta, C_group for good (i,j)
      do 101 i=1,nvec
         avi=sqrt(float(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)))
c     (avxi,avyi,avzi) unit vector for peptide bond "i"
         avxi=vx(i)/avi
         avyi=vy(i)/avi
         avzi=vz(i)/avi
         do 102 j=1,nvec
            avj=sqrt(float(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)))
c     (avxj,avyj,avzj) unit vector for peptide bond "j"
            avxj=vx(j)/avj
            avyj=vy(j)/avj
            avzj=vz(j)/avj
c     (ax,ay,az) unit vector resulting from the normalization of the
c     superposition of (avxi,avyi,avzi) and (avxj,avyj,avzj)
            ax=avxi+avxj
            ay=avyi+avyj
            az=avzi+avzj
            aaa=sqrt(ax*ax+ay*ay+az*az)
            ax=ax/aaa
            ay=ay/aaa
            az=az/aaa
c     (bx,by,bz)==(hbx(i,j),hby(i,j),hbz(i,j)) unit vector perpendicular
c     to (avxi,avyi,avzi) and (avxj,avyj,avzj), which is direction in
c     which hydrogen bond at the central CA may happen
            bx=avyi*avzj-avzi*avyj
            by=avzi*avxj-avxi*avzj
            bz=avxi*avyj-avyi*avxj
            bbb=sqrt(bx*bx+by*by+bz*bz)
            bx=bx/bbb
            by=by/bbb
            bz=bz/bbb
            hbx(i,j)=bx
            hby(i,j)=by
            hbz(i,j)=bz
c     (cx,cy,cz)==(cax(i,j),cay(i,j),caz(i,j)) unit vector resulting
c     from the normalization of the substraction of (avxj,avyj,avzj) to
c     (avxi,avyi,avzi)
            cx=avxi-avxj
            cy=avyi-avyj
            cz=avzi-avzj
            ccc=sqrt(cx*cx+cy*cy+cz*cz)
            cx=cx/ccc   
            cy=cy/ccc   
            cz=cz/ccc   
            cax(i,j)=cx 
            cay(i,j)=cy 
            caz(i,j)=cz             
CM            print *,'cax etc',i,j,cax(i,j),cay(i,j),caz(i,j),angle(i,j)
CM            print *,'hbx etc',i,j,hbx(i,j),hby(i,j),hbz(i,j),angle(i,j)
c     side-chain coordinate from C_a
            do k=0,19
c     contiguos peptide bonds "i" and "j" are in a helical or turn-like
c     conformation
               if(angle(i,j).lt.105) then
c     (gx,gy,gz)(i,j,k) are coordinates of vector going from CA to SG
c     for residue type "k" when the CA is flanked by peptide bonds "i"
c     and "j"
                  gx(i,j,k)=(axalf(k)*ax+ayalf(k)*bx+azalf(k)*cx)/0.87
                  gy(i,j,k)=(axalf(k)*ay+ayalf(k)*by+azalf(k)*cy)/0.87
                  gz(i,j,k)=(axalf(k)*az+ayalf(k)*bz+azalf(k)*cz)/0.87
c     peptide bonds are in extended conformation
               else
                  gx(i,j,k)=(axbet(k)*ax+aybet(k)*bx+azbet(k)*cx)/0.87
                  gy(i,j,k)=(axbet(k)*ay+aybet(k)*by+azbet(k)*cy)/0.87
                  gz(i,j,k)=(axbet(k)*az+aybet(k)*bz+azbet(k)*cz)/0.87
               endif
CM               print *,'gx gy gx',i,j,k,gx(i,j,k),gy(i,j,k),gz(i,j,k)
            enddo
 102     continue
 101  continue
      return
      end
ccccc prepare_beta finished


ccccc Compute the secondary fragment biases
c     check local secondary structure from 'seq.dat'
c     if it is beta-structure in [i,i+6],  frga(i)=19.1/0.87;
c     if it is alpha-structure in [i,i+7], frgb(i)=10.5/0.87;
      subroutine prepare_frg(Lch,frga,frgb,sec)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      dimension frga(ndim),frgb(ndim)
      dimension sec(ndim)

CM      common/fr/frga(ndim),frgb(ndim)
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/lengths/Lch,Lch1,Lch2

      print *,'In prepare_frg',Lch
      do i=1,Lch
         frga(i)=0.0
         frgb(i)=0.0
      enddo

      do i=1,Lch-7
c     q accumulates values of secondary structure in an eigth residue
c     fragment
         q=0
c     we look at segments of size eigth residues (seven peptide bonds)
         do j=i,i+7            
            if(sec(j).eq.2) q=q+2
         enddo
c     8 continue alpha-residues
         if(q.eq.16)then
c     distance for 7 peptide bonds when in alpha conformation
            frga(i)=10.5/0.87 
         endif
CM         print *,'frga',i,frga(i)
      enddo

      do i=1,Lch-6
         q=0
c     we looks at segments of size six residues (five peptide bonds)
         do j=i+1,i+5
            if(sec(j).eq.4) q=q+4
         enddo
c     5 continue beta-residues
         if(q.eq.20)then        
            if(sec(i).ne.2.and.sec(i+6).ne.2)then
c     distance for 6 peptide bonds in extended conformation
               frgb(i)=19.1/0.87
            endif
         endif
CM         print *,'frgb',i,frgb(i)
      enddo
      return
      end
ccccc  prepare_frg finished


ccccc compute contact order
      subroutine get_acorder(Lch,sec,acorder)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)

      dimension sec(ndim)

CM      common/seqe/seq(ndim),sec(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/lengths/Lch,Lch1,Lch2
     
      print *,'In get_acorder',Lch
 
c     number of residues predicted in helical conformation
      n_H=0                     
c     number of residues predicted in extended conformation
      n_E=0                     
      do i=1,Lch
         if(sec(i).eq.2) n_H=n_H+1
         if(sec(i).eq.4) n_E=n_E+1
      enddo
c     structure with no secondary structure, use alpha=H/E
      if(n_H+n_E.lt.2)then      
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
c     use alpha=aH+bE
      else           
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
c     predicted contact order
      acorder=alph*Lch
      print *,'acorder',acorder
      return
      end
ccccc get_acorder finished



ccccc print out initial parameters
      subroutine write_parameter(Lch,ncycle,N_rep,phot,Ndis,Ncom,atemp2,
     +                           atemp1,exc,er1,er3,er4,er5,er6,er7,
     +                           eh1c,random,acorder,switch)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      parameter(nrep=100)      

CM      common/commonuse1/random,ncycle
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      common/lengths/Lch,Lch1,Lch2
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/distres/er4,es3c
CM      common/resnumber/Ncom,Ndis,accur
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      common/looks/exc,exc1,exc2
CM      common/pair1/eh2,eh1b,eh1c

CMM      common/seed1/no

CM                  print *,'In write_para NO=',no

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
      write(20,*)'................er3........',er3
      write(20,*)'................er4........',er4
      write(20,*)'................er5........',er5
      write(20,*)'................er6........',er6
      write(20,*)'................er7........',er7
      write(20,*)'................eh1c.......',eh1c
      write(20,*)
      write(20,*)'initial random number seed...........',random
      write(20,*)'contact order........................',acorder
      write(20,*)
      write(20,*)'the first arandom number.............',aranzy()
CM                  print *,'NO=',no
CM      +                                  aranzy(no,ix1,ix2,ix3,r)
      write(20,*) 
      write(20,*)'*******************************************'
      if(switch.gt.1)then
         write(20,*)'This is RS running of threading structure!'
      else
         write(20,*)'This is a normal simulation.'
      endif
      write(20,*)'*******************************************'
      write(20,*)
      return
      end
ccccc  write_parameter  finished

ccccc set movement percentage
c     CPU time of each movement
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
      subroutine set_move_ratio(h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,
     +                          hend,bh2,bh3s,bh3d,bh4s,bh4d,
     +                          bh5s,bh5d,bh6,bhendn)
      implicit integer(i-z)
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      common/moveratio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
CM      common/moveratio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
CM      common/moveratio3/bh7a,bh7b,bhendn,bhendc

      print *,'In set_move_ratio',
     +         h2+h3s+h3d+h4s+h4d+h5s+h5d+h6+hend

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
      print *,'bh etc',
     +   bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6,bhendn,bhendc
      return
      end
ccccc set_move_ratio finished

ccccc prepare 2-bond move, i.e. calculate v21(tx,ty,tz,i), v22(tx,ty,tz,i).
c v21(tx,ty,tz,i) --- 1th vector of i'th path from (0,0,0) to (tx,ty,tz)
c v22(tx,ty,tz,i) --- 2th vector of i'th path from (0,0,0) to (tx,ty,tz)
c it will be dangerous if number of used variable is more than 3,000,000
c i.e. the usage of memory of CUP can not beyond 90%.
      subroutine prepare_move2(Nw,goodc,vx,vy,vz,Np2,v21,v22,vector,
     +                         w21,w22)
      implicit integer(i-z)
      parameter(nvec=312)

      logical goodc(nvec,nvec)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)


      print *,'In prepare_move2'
      do i=-10,10
         do j=-10,10
            do k=-10,10
               Nw(i,j,k)=0
            enddo
         enddo
      enddo

      max=0
      maxa=0
      nnn=0
      mmm=0
      do 101 i=1,nvec
         do 102 j=1,nvec
c     angle between peptide bonds "i" and "j" in [65,165]
            if(goodc(i,j))then
c     (ijx,ijy,ijz) vector from superposition of peptide bonds "i" and
c     "j"
               ijx=vx(i)+vx(j)
               ijy=vy(i)+vy(j)
               ijz=vz(i)+vz(j)
c     Nw(ijx,ijy,ijz) number of ways one can superimpose two peptide
c     bonds and obtain vector (ijx,ijy,ijz). Note: the two peptide bonds
c     must have a "good" angle.
               Nw(ijx,ijy,ijz)=Nw(ijx,ijy,ijz)+1 !based on r
c     w21(ijx,ijy,ijz, l) and w22(ijx,ijy,ijz, l) are indexes of two
c     peptide bonds whose superposition results in vector
c     (ijx,ijy,ijz). "l" is a running index from 1 to Nw(ijx,ijy,ijz)
               w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=i
               w22(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=j
CM               print *,'w21 w22',
CM     +          w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz)),
CM     +          w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz))
               if(maxa.lt.Nw(ijx,ijy,ijz))maxa=Nw(ijx,ijy,ijz)
c     based on i,j:
               Np2(i,j)=0
               do 103 ii=1,nvec
                  jx=ijx-vx(ii)
                  jy=ijy-vy(ii)
                  jz=ijz-vz(ii)
                  jr=jx*jx+jy*jy+jz*jz
                  if(jr.ge.14.and.jr.le.25)then
c     "ii" and "jj" are two other peptide bonds whose superposition also
c     gives us (ijx,ijy,ijz)
                     jj=vector(jx,jy,jz)
                     if(goodc(ii,jj))then
c     Np2(i,j) is the number of ways one can obtain the superposition of
c     both peptide bonds
                        Np2(i,j)=Np2(i,j)+1
c     v21(i,j, l) and v22(i,j, l) are indexes of two peptide bonds whose
c     superposition gives the same vector than the superposition of "i"
c     and "j". "l" is a running index from 1 to Np2(i,j) to list all
c     possible peptide bonds. Note: the peptide bonds must make a "good"
c     angle
                        v21(i,j,Np2(i,j))=ii
                        v22(i,j,Np2(i,j))=jj
                        mmm=mmm+1 !number of total memory occupied by vx21.
                     endif
                  endif
CM            print *,'Np2 v21 v22',i,j,ii,jj,Np2(i,j),v21(i,j,Np2(i,j)),
CM     +                                v22(i,j,Np2(i,j))    
 103           continue
               if(max.lt.Np2(i,j))max=Np2(i,j)
               nnn=nnn+1        !number of possible pairs
            endif

 102     continue
 101  continue

ccccc among all 312*312=97344 pairs, 67272 pairs are legal (~70%).  all
c     Np2(i,j)>=2, i.e. there are at least one other path for any pair.
c     <Np2(i,j)>=27.
cccc  the following are the cases without goodc on prefabracated pairs:
c              N_v  N_pare N_pp   mmm
c     [14,25], 312, 97344, 67272, 1830000 (25% memory, within memory)
ccccc the following is cases without goodc limitation:
c     [14,25], 306, 93636, 3872214 (90% memory, beyond memory)
c     [14,24], 282, 79524, 2923758 (80% memory, beyond memory)
      return
      end
ccccc prepare_move2 finished





      subroutine prepare_move2_1(Nw,goodc,vx,vy,vz,Np2,v21,v22,vector,
     +                         w21,w22)
      implicit integer(i-z)
      parameter(nvec=312)

      logical goodc(nvec,nvec)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)

CM      common/logica/goodc
CM      logical goodc(nvec,nvec)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)

CM      dimension w21(-10:10,-10:10,-10:10,46)
CM      dimension w22(-10:10,-10:10,-10:10,46)

      print *,'In prepare_move2'
      do i=-10,10
         do j=-10,10
            do k=-10,10
               Nw(i,j,k)=0
            enddo
         enddo
      enddo

      max=0
      maxa=0
      nnn=0
      mmm=0
      do 101 i=1,nvec
         do 102 j=1,nvec
c     angle between peptide bonds "i" and "j" in [65,165]
            if(goodc(i,j))then
c     (ijx,ijy,ijz) vector from superposition of peptide bonds "i" and
c     "j"
               ijx=vx(i)+vx(j)
               ijy=vy(i)+vy(j)
               ijz=vz(i)+vz(j)
c     Nw(ijx,ijy,ijz) number of ways one can superimpose two peptide
c     bonds and obtain vector (ijx,ijy,ijz). Note: the two peptide bonds
c     must have a "good" angle.
               Nw(ijx,ijy,ijz)=Nw(ijx,ijy,ijz)+1 !based on r
c     w21(ijx,ijy,ijz, l) and w22(ijx,ijy,ijz, l) are indexes of two
c     peptide bonds whose superposition results in vector
c     (ijx,ijy,ijz). "l" is a running index from 1 to Nw(ijx,ijy,ijz)
               w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=i
               w22(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=j
               print *,'w21 w22',
     +          w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz)),
     +          w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz))
               if(maxa.lt.Nw(ijx,ijy,ijz))maxa=Nw(ijx,ijy,ijz)
c     based on i,j:
               Np2(i,j)=0
               do 103 ii=1,nvec
                  jx=ijx-vx(ii)
                  jy=ijy-vy(ii)
                  jz=ijz-vz(ii)
                  jr=jx*jx+jy*jy+jz*jz
                  if(jr.ge.14.and.jr.le.25)then
c     "ii" and "jj" are two other peptide bonds whose superposition also
c     gives us (ijx,ijy,ijz)
                     jj=vector(jx,jy,jz)
CM                     write(*,*),'HERE',ii,jj
CM                     write(*,*),'HERE',goodc(ii,jj)
                     if(goodc(ii,jj))then
c     Np2(i,j) is the number of ways one can obtain the superposition of
c     both peptide bonds
                        Np2(i,j)=Np2(i,j)+1
c     v21(i,j, l) and v22(i,j, l) are indexes of two peptide bonds whose
c     superposition gives the same vector than the superposition of "i"
c     and "j". "l" is a running index from 1 to Np2(i,j) to list all
c     possible peptide bonds. Note: the peptide bonds must make a "good"
c     angle
                        v21(i,j,Np2(i,j))=ii
                        v22(i,j,Np2(i,j))=jj
                        mmm=mmm+1 !number of total memory occupied by vx21.
                     endif
                  endif
          print *,'Np2 v21 v22',i,j,ii,jj,Np2(i,j),v21(i,j,Np2(i,j)),
     +                                v22(i,j,Np2(i,j))
 103           continue
               if(max.lt.Np2(i,j))max=Np2(i,j)
               nnn=nnn+1        !number of possible pairs
            endif
 102     continue
 101  continue

ccccc among all 312*312=97344 pairs, 67272 pairs are legal (~70%).  all
c     Np2(i,j)>=2, i.e. there are at least one other path for any pair.
c     <Np2(i,j)>=27.
cccc  the following are the cases without goodc on prefabracated pairs:
c              N_v  N_pare N_pp   mmm
c     [14,25], 312, 97344, 67272, 1830000 (25% memory, within memory)
ccccc the following is cases without goodc limitation:
c     [14,25], 306, 93636, 3872214 (90% memory, beyond memory)
c     [14,24], 282, 79524, 2923758 (80% memory, beyond memory)
      return
      end
ccccc prepare_move2 finished


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     from 'chain.dat':
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine template_initial2(exc,i_thr0,Lch,q,cx0,cy0,cz0,ssp,M_a,
     +                             M_i,M_f,L_a,L_cut,n_i,n_f,cx,cy,cz,
     +                             ras,nfr_i,nfr_f,mv,angle00,angle0,
     +                             d_xyz00,d_xyz0,N_rep,exrep,eyrep,
     +                             ezrep,ecxrep,ecyrep,eczrep,ebxrep,
     +                             ebyrep,ebzrep,seq,axalf,ayalf,azalf,
     +                             axbet,aybet,azbet,egxrep,egyrep,
     +                             egzrep,ex,ey,ez,ax,ay,az,x,y,z,goodc,
     +                             vx,vy,vz,xrep,yrep,zrep,sec,
     +                             icarep,vector,sequ,vvv,ex0,ey0,ez0,
     +                             i_bigbond,
     +                             nfl,Mend_N,Mend,Mend_C,nfl6,
     +                             ras6,nfl5,ras5,nfl4,ras4,nfl3,ras3,
     +                             nfl2,ras2,Lch2)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)

      logical goodc(nvec,nvec)
      dimension x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension cx(1000),cy(1000),cz(1000)
      dimension seq(ndim)
      dimension xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension icarep(ndim,nrep)
      dimension exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      dimension egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      dimension ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      dimension ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      dimension vvv(ndim,ndim)
      dimension cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)
      dimension axalf(0:19),ayalf(0:19),azalf(0:19)
      dimension axbet(0:19),aybet(0:19),azbet(0:19)
      dimension ras(ndim)
CM      real nfl
      dimension d_xyz00(ndim),angle00(ndim)
      dimension nfr_i(ndim),nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim)
CM      dimension egx(ndim),egy(ndim),egz(ndim)
      character*3 sequ
      dimension sequ(ndim)
      dimension ex0(ndim),ey0(ndim),ez0(ndim)
      dimension mv(ndim),sec(ndim)

      character sign(ndim)
      dimension i_chunk(ndim)
      dimension ras2(ndim)
      dimension ras3(ndim)
      dimension ras4(ndim)
      dimension ras5(ndim)
      dimension ras6(ndim)



CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/xyzrs/xrs(ndim,40,40),yrs(ndim,40,40),zrs(ndim,40,40)
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/rs/i_thr0
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      character c1*4,c2*2,c3*3,aaaaa*10,text,text1*22
CM      character*80 head
CM      character sign(ndim)
CM      common/initialstr/cx(1000),cy(1000),cz(1000)
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/looks/exc,exc1,exc2
CM
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      common/distres/er4,es3c
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM
CM      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
CM      common/sw3/icarep(ndim,nrep)
CM      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
CM      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
CM      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
CM      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
CM      common/excluded/vvv(ndim,ndim)
CM
CM      common/consensus1/cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)
CM
CM      common/mng/m_g(100)
CM      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
CM      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)

CM      common/chain0/ras(ndim),nfl
CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/ssp/ssp

CM      common/chainm/mv(ndim)

CM      character*3 sequ
CM      common/aminoacid/sequ(ndim)
CM      common/bigbond/i_bigbond,teco
CM      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)

      dimension n_i(1000),n_f(1000)
      dimension cx0b(1000),cy0b(1000),cz0b(1000)
      dimension cx0g(1000),cy0g(1000),cz0g(1000)
      dimension ax(ndim),ay(ndim),az(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)
      dimension M_i(1000),M_f(1000)
      dimension mf(0:ndim)
      dimension q_bk(ndim)


      print *,'In template_initial2',Lch,exc,i_thr0

******************************************************
c     set-up the template to use
******************************************************
      exc_eff=exc               !for excluded volumn
      M_consensus=0             !not use consensus
      if(i_thr0.eq.0)then       !means use consensus of top-2 template
         i_thr0=1               !if failed, using the first template
         M_consensus=1          !use consensus in the following
      endif

****************************************************
c     read q(i), cx0(i) from 'init.dat'
****************************************************
      rewind(24)
      read(24,*)n_thr           !Number of real templates in 'init.dat'
      if(i_thr0.gt.n_thr)then
         write(20,*)'without threading structure at this i_thr0!'
         write(*,*)'without threading structure at this i_thr0!'
         stop
      endif
      do k=1,i_thr0
         do i=1,Lch
            q(i)=0
         enddo
         read(24,*)Nal
         do i=1,Nal
            read(24,*)ii,a1,a2,a3
            cx0(ii)=a1
            cy0(ii)=a2
            cz0(ii)=a3
            q(ii)=1
         enddo
      enddo
 1237 format(A22,I4,A4,3F8.3)
c      if(M_consensus.eq.1)call get_consensus !get q(i),cx0(i) from consensus
c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

**********************************************************
c     using predicted secondary structures
**********************************************************
      if(ssp.eq.1)then
CM         call secondary         !redifine q(i) and cx0(i)
           call secondary(Lch,cx0,cy0,cz0,q,L_cut)
      endif
      
********************************************************
c     remove small segments, so that we have a better random walk
********************************************************
      M_a=0                     !number of aligned segments
      q(0)=0
      do i=1,Lch
         if(q(i).ne.q(i-1).and.q(i).eq.1)then
            M_a=M_a+1
            M_i(M_a)=i
         endif
         if(q(i).eq.1)M_f(M_a)=i
      enddo
      do i=1,M_a
         L_a=M_f(i)-M_i(i)+1
         if(L_a.le.L_cut)then
            do j=M_i(i),M_f(i)
               q(j)=0
            enddo
         endif
      enddo
c^^^^^^^^^^^ remove small segment finished ^^^^^^^^^^^^^

********************************************************
c     check and find GAP, only according to q(i)
********************************************************
      N_g=0                     !number of gaps
      q(0)=1
      do i=1,Lch
         if(q(i).ne.q(i-1).and.q(i).eq.0)then
            N_g=N_g+1
            n_i(N_g)=i          !initial point of the gap
         endif
         if(q(i).eq.0)n_f(N_g)=i !final point of the gap
      enddo
c^^^^^^^^^^^^^^^ check gap finished ^^^^^^^^^^^^^^^^^

      n_check_bond=0
 103  continue                  !re-walk and try to make dis(i,i+1)<7A
********************************************************
c     fill GAP: cx0(i) -> cx(i)
********************************************************
      n_walk=0                  !for number of reject by excluded volumn
 70   n_walk=n_walk+1
      if(n_walk.gt.1000)then
         exc_eff=exc_eff*0.99
      endif
      if(n_walk.gt.10000)then
         write(20,*)'unsolvable structure',pass,i,exc_eff
         write(*,*)'unsolvable structure',pass,i,exc_eff
         stop
      endif
      do i=1,Lch
         if(q(i).eq.1)then
            cx(i)=cx0(i)
            cy(i)=cy0(i)
            cz(i)=cz0(i)
         else
            cx(i)=1000000.      !for checking excluded volumn
            cy(i)=1000000.
            cz(i)=1000000.
         endif
      enddo
      do i=2,n_g
         i1=n_i(i)
         i2=n_f(i)
CM         call connect(i1,i2,pass) !fill missed cooridinates, no move others
           call connect(i1,i2,pass,cx,cy,cz,Lch,amcheck_dis) !fill missed cooridinates, no move others
         if(pass.ge.3)goto 70             !re-walk
      enddo
      if(n_g.ge.1)then
         i1=n_i(1)
         i2=n_f(1)
CM         call connect(i1,i2,pass)
         call connect(i1,i2,pass,cx,cy,cz,Lch,amcheck_dis) !fill missed cooridinates, no move others
         if(pass.ge.3)goto 70   !re-walk
      endif
*^^^^^^^^^^^^^^Fill gap done, cx(i) is continuous ^^^^^^^^^^^^^^^^^^^^^^^

**********************************************************************
c     decide moveable points nfl, ras(i), according to 3 factors:
**********************************************************************
***   1: ras => gap +- 2, smallest moveable-length is 3
      nfl=0                     !number of moveable points
      do i=1,Lch
         ras(i)=0               !residue name of i'th moveable point.
      enddo
      do i=1,n_g
         i1=n_i(i)-2
         i2=n_f(i)+2
         if(i1.lt.1)i1=1
         if(i2.gt.Lch)i2=Lch
         do ii=i1,i2
            nfl=nfl+1           !number of moveable residues
            ras(nfl)=ii         !locations
         enddo
      enddo
      if(nfl.gt.0)call sort_ras(nfl,ras)
      if(i_bigbond.eq.1)then
***   2a: ras => bigbond on  +-2 for templates:
         do i=2,Lch
            if(q(i-1).eq.1.and.q(i).eq.1)then
               bondlen=di(cx(i-1),cy(i-1),cz(i-1),cx(i),cy(i),cz(i))
               if(bondlen.gt.4.6)then
                  ik_f=i+1      !ending point
                  ik_i=i-2      !starting point of big-bond
                  if(ik_i.lt.1)ik_i=1
                  if(ik_f.gt.Lch)ik_f=Lch
                  do ik=ik_i,ik_f
                     nfl=nfl+1
                     ras(nfl)=ik
                  enddo
               endif
            endif
         enddo
         if(nfl.gt.0)call sort_ras(nfl,ras)
      else
***   2b: ras => bigbond +-1
         bond_max=0
         do i=2,Lch
            j=i-1               !starting point of big-bond
            bondlen=di(cx(i),cy(i),cz(i),cx(j),cy(j),cz(j))
            if(bondlen.gt.bond_max) bond_max=bondlen
            if(bondlen.gt.5.0)then
               do t=i,Lch
                  dh=di(cx(j),cy(j),cz(j),cx(t),cy(t),cz(t))
                  adh=dh/float(t-j)
                  if(adh.lt.3.498)goto 23 !can be walked to
               enddo
 23            ik_f=t+1         !ending point + 1
               ik_i=j-1         !starting point of big-bond - 1
               if(ik_i.lt.1)ik_i=1
               if(ik_f.gt.Lch)ik_f=Lch
               do ik=ik_i,ik_f
                  nfl=nfl+1
                  ras(nfl)=ik
               enddo
            endif
         enddo
         if(nfl.gt.0)call sort_ras(nfl,ras)
      endif
***   3: ras0 => smallpiece:
      do i=1,Lch
         mf(i)=0                !frozen
      enddo
      do i=1,nfl
         mf(ras(i))=1           !moveable points
      enddo
      nfr=0                     !number of frozen fragments
      mf(0)=1
      do i=1,Lch
         if(mf(i).ne.mf(i-1).and.mf(i).eq.0)then
            nfr=nfr+1
            nfr_i(nfr)=i        !starting point of nfr'th fragment
         endif
         if(mf(i).eq.0)nfr_f(nfr)=i !ending point of nfr'th fragment
      enddo
c     'nfr' frozen fragments in [nfr_i(i),nfr_f(i)], convert small piece ---->
      do i=1,nfr
         l_fr=nfr_f(i)-nfr_i(i)+1 !length of the frozen fragment
         if(l_fr.le.L_cut)then
            do j=nfr_i(i),nfr_f(i)
               nfl=nfl+1
               ras(nfl)=j
            enddo
         endif
      enddo
      if(nfl.gt.0)call sort_ras(nfl,ras)
*^^^^^^^^^^^^^^^^^^^^ nfl, ras(i) finished ^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
c     decide nfl2, ras2(i), mv(i):
*****************************************************************
      if(nfl.lt.1)then
         write(20,*)'without moveable points in this template!'
         write(*,*)'without moveable points in this template!'
      endif
CM      call move_point           !decide movement set from ras(i)
      call move_point(Lch,mv,nfl,ras,Mend_N,Mend,Mend_C,
     +                      nfl6,ras6,nfl5,ras5,nfl4,ras4,nfl3,
     +                      ras3,nfl2,ras2)
*^^^^^^^^^^ decision of movement set finished ^^^^^^^^^^^^^^^^^^^

****************************************************************************
c     decide nfr_i(i), nfr_f(i), i=1,...,nfr, according to nfl,ras(nfl):
****************************************************************************
      do i=1,Lch
         mf(i)=0
         sign(i)="f"            !frozen
      enddo
      do i=1,nfl
         mf(ras(i))=1
         sign(ras(i))="m"       !moveable points
      enddo
      mf(0)=1
      nfr=0                     !number of frozen fragments
      do i=1,100
         nfr_i(i)=0
         nfr_f(i)=0
      enddo
      do i=1,Lch
         if(mf(i).ne.mf(i-1).and.mf(i).eq.0)then
            nfr=nfr+1
            nfr_i(nfr)=i         !starting point of nfr'th fragment
         endif
         if(mf(i).eq.0)nfr_f(nfr)=i !ending point of nfr'th fragment
      enddo
c^^^^^^^^^^^^^^^^^^^ nfr, nfr_i(i), nfr_f(i) finished ^^^^^^^^^^^^^^^^^^

****************************************************************
c     redefine mv(i) for frozen fragments (for excluded volumn)
****************************************************************
      do i=1,nfr
         do j=nfr_i(i),nfr_f(i)
            mv(j)=-i
         enddo
      enddo
c^^^^^^^^^^^^^^^^^ mv(i) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****************************************************************
c     translation distance d_xyz00(i), rotation angle00(i)
      do i=1,nfr
         siz=nfr_f(i)-nfr_i(i)+1
         angle00(i)=angle0*10/float(siz)
         d_xyz00(i)=d_xyz0*10/float(siz)
         if(siz.gt.40)d_xyz00(i)=d_xyz00(i)/2.0
      enddo
c^^^ d_xyz00(i), angle00(i) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****************************************************************
c     Record frozen parts:
c     there are 4 arraies needed to record and transfer:
c     common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)    CA
c     common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) SG
c     common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) cc
c     common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) bb
****************************************************************
      cx0(0)=cx(1)+(cx(2)-cx(3))
      cy0(0)=cy(1)+(cy(2)-cy(3))
      cz0(0)=cz(1)+(cz(2)-cz(3))
      cx0(Lch+1)=cx(Lch)+(cx(Lch1)-cx(Lch2))
      cy0(Lch+1)=cy(Lch)+(cy(Lch1)-cy(Lch2))
      cz0(Lch+1)=cz(Lch)+(cz(Lch1)-cz(Lch2))
      do itemp=1,N_rep
         do i=1,Lch
            if(mv(i).lt.0)then  !frozen points
c     Ca:
               exrep(i,itemp)=cx0(i)/0.87
               eyrep(i,itemp)=cy0(i)/0.87
               ezrep(i,itemp)=cz0(i)/0.87
c     uniform vector:
               amx=cx0(i)-cx0(i-1)
               amy=cy0(i)-cy0(i-1)
               amz=cz0(i)-cz0(i-1)
               aaa=sqrt(amx**2+amy**2+amz**2)
               amx=amx/aaa      !ax(i)-ax(i-1)
               amy=amy/aaa
               amz=amz/aaa
               apx=cx0(i+1)-cx0(i)
               apy=cy0(i+1)-cy0(i)
               apz=cz0(i+1)-cz0(i)
               aaa=sqrt(apx**2+apy**2+apz**2)
               apx=apx/aaa      !ax(i+1)-ax(i)
               apy=apy/aaa
               apz=apz/aaa
               ang=acos(-(amx*apx+amy*apy+amz*apz))*180/3.1415926 !
               aaax=amx+apx
               aaay=amy+apy
               aaaz=amz+apz
               aaa=sqrt(aaax**2+aaay**2+aaaz**2)
               aaax=aaax/aaa
               aaay=aaay/aaa
               aaaz=aaaz/aaa
c     cc:
               ccx=amx-apx
               ccy=amy-apy
               ccz=amz-apz
               aaa=sqrt(ccx**2+ccy**2+ccz**2)
               ccx=ccx/aaa
               ccy=ccy/aaa
               ccz=ccz/aaa
               ecxrep(i,itemp)=ccx !cc
               ecyrep(i,itemp)=ccy
               eczrep(i,itemp)=ccz
c     Hb:
               bbx=amy*apz-amz*apy
               bby=amz*apx-amx*apz
               bbz=amx*apy-amy*apx
               aaa=sqrt(bbx**2+bby**2+bbz**2)
               bbx=bbx/aaa
               bby=bby/aaa
               bbz=bbz/aaa
               ebxrep(i,itemp)=bbx !Hb
               ebyrep(i,itemp)=bby
               ebzrep(i,itemp)=bbz
c     SG:
               k=seq(i)
               if(ang.lt.105)then !alpha
                  dx=(axalf(k)*aaax+ayalf(k)*bbx+azalf(k)*ccx)/0.87
                  dy=(axalf(k)*aaay+ayalf(k)*bby+azalf(k)*ccy)/0.87
                  dz=(axalf(k)*aaaz+ayalf(k)*bbz+azalf(k)*ccz)/0.87
               else
                  dx=(axbet(k)*aaax+aybet(k)*bbx+azbet(k)*ccx)/0.87
                  dy=(axbet(k)*aaay+aybet(k)*bby+azbet(k)*ccy)/0.87
                  dz=(axbet(k)*aaaz+aybet(k)*bbz+azbet(k)*ccz)/0.87
               endif
               egxrep(i,itemp)=exrep(i,itemp)+dx
               egyrep(i,itemp)=eyrep(i,itemp)+dy
               egzrep(i,itemp)=ezrep(i,itemp)+dz
            endif
         enddo
      enddo
c     current ex(i) for check of excluded volumn:
      do i=1,Lch
         if(mv(i).lt.0)then     !frozen points
            ex(i)=cx0(i)/0.87
            ey(i)=cy0(i)/0.87
            ez(i)=cz0(i)/0.87
         endif
      enddo
c^^^^^^^^^^^^^^^ record of frozen-fragment are set^^^^^^^^^^^^^^^^^^^^^^

*************************************************************************
c     project chain onto lattices, decide (x,y,z):
*************************************************************************
***   ax(i)=cx(i)/0.87, for the decision of (x,y,z):
      do i=1,Lch
         ax(i)=cx(i)/0.87       !C_alpha scaled by 0.87
         ay(i)=cy(i)/0.87
         az(i)=cz(i)/0.87
      enddo
c     (ax,ay,az) into (x,y,z)----------------------->
      x(1)=nint(ax(1))
      y(1)=nint(ay(1))
      z(1)=nint(az(1))
      xm=x(1)
      ym=y(1)
      zm=z(1)
      do 101 i=2,Lch
         dis_min=1000000.       !minimum distance between c_alpha and lattice
         j_ch=0                 !chosen vector
         do 100 j=1,nvec
            if(i.gt.2)then      !check good neighbor
               if(.not.goodc(jm,j))goto 100
            endif
            x_tmp=xm+vx(j)
            y_tmp=ym+vy(j)
            z_tmp=zm+vz(j)
c     check excluded volumn---->
            if(mv(i).gt.0)then  !on-lattice
               do m=1,Lch       !!!!! off-lattice part.
                  if(abs(i-m).ge.3.and.mv(m).lt.0)then
                     disaa=(x_tmp-ex(m))**2+(y_tmp-ey(m))**2+
     $                    (z_tmp-ez(m))**2
                     if(disaa.lt.exc_eff) goto 100
                  endif
               enddo
               do m=1,i-3       !!!!! on-lattice part.
                  if(mv(m).gt.0)then
                     disaa=(x_tmp-x(m))**2+(y_tmp-y(m))**2+
     $                    (z_tmp-z(m))**2
                     if(disaa.lt.exc_eff) goto 100
                  endif
               enddo
            endif
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^
            dis=(ax(i)-x_tmp)**2+(ay(i)-y_tmp)**2+(az(i)-z_tmp)**2
            if(dis.lt.dis_min)then
               j_ch=j
               x_ch=x_tmp
               y_ch=y_tmp
               z_ch=z_tmp
               dis_min=dis
            endif
 100     continue
         if(j_ch.lt.1)goto 70   !refill the gaps
         x(i)=x_ch              !get (x(i),y(i),z(i)) here
         y(i)=y_ch
         z(i)=z_ch
         jm=j_ch
         xm=x(i)
         ym=y(i)
         zm=z(i)
 101  continue
c^^^^^^^^^^^^^^^^^^^project lattice done ^^^^^^^^^^^^^^^^^^^^^^^

***********************************************************
*     check the bond-length
***********************************************************
      do i=1,Lch1
         r2=
     $        (nint(aax(i,mv,x,ex))-nint(aax(i+1,mv,x,ex)))**2+
     $        (nint(aay(i,mv,y,ey))-nint(aay(i+1,mv,y,ey)))**2+
     $        (nint(aaz(i,mv,z,ez))-nint(aaz(i+1,mv,z,ez)))**2
         if(r2.gt.65)then       !dis>7A
            n_check_bond=n_check_bond+1
            if(n_check_bond.lt.100)goto 103
         endif
c         write(*,*)i,r2,sqrt(float(r2))*.87
      enddo
cccc  

*************************************************************************
c     record replcas of (x,y,z), icarep:
*************************************************************************
      do itemp=1,N_rep
         do i=1,Lch
            xrep(i,itemp)=x(i)  !x(i)
            yrep(i,itemp)=y(i)
            zrep(i,itemp)=z(i)
            if(i.lt.Lch)then
               wx=x(i+1)-x(i)   !ica(i)
               wy=y(i+1)-y(i)
               wz=z(i+1)-z(i)
               icarep(i,itemp)=vector(wx,wy,wz)
            endif
         enddo
         icarep(Lch,itemp)=icarep(Lch2,itemp) !just for unity
      enddo

CM       call get_vvv !get vvv(i,j). vvv(i,j)>0 check; vvv(i,j)<0 not check
      call get_vvv(Lch,vvv,mv,exc) !get vvv(i,j). vvv(i,j)>0 check; vvv(i,j)<0 not check

***********************************************************************
c      decide the frequence of bulk movements: f=1/n_fra
***********************************************************************
      n_fra=nint(8*sqrt(float(Lch)/100.0)) !move bulk once in n_fra local move
cccccc n_fra decrease with number of frozen residues
      nnn=0
      do i=1,Lch
        if(mv(i).lt.0)nnn=nnn+1
      enddo
      n_fra=nint(n_fra*(float(Lch)/(2.0*float(nnn)+0.001)))
cccccc n_fra decrease with number of frozen fragments:
      np=1
      if(nfr.eq.1)np=2
      n_fra=nint(n_fra*(2.8/(float(nfr)+0.001))**np)
      if(n_fra.lt.1) n_fra=1 !nfr big and nnn big
      if(n_fra.gt.20) n_fra=20 !nfr small and nnn small
c^^^^^^^^^^^ n_fra decided ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

********************************************************************
*     output the threading results:
********************************************************************
      write(20,*)'Lch=',Lch
      write(20,*)"#align=",Nal,'   #unalign=',Lch-Nal,"n_g=",n_g
      write(20,*)'#frozen=',Lch-nfl,'    #moveable=',nfl
      write(20,*)
      write(20,*)'number of frozen pieces: nfr=',nfr
      do i=1,nfr
         siz=nfr_f(i)-nfr_i(i)+1
         write(20,42)i,siz,nfr_i(i),nfr_f(i),d_xyz00(i),angle00(i)*57.3
      enddo
 42   format(i5,i5,' [',i3,',',i3,']',f8.3,f8.3)
      write(20,*)'---number of local movement for each bulk move:',n_fra
      write(20,*)
      write(20,*)'#chain fix/move  SEC   SEQ'
      do i=1,Lch
         write(20,41)i,sign(i),sec(i),sequ(i)
      enddo
 41   format(i7,a5,i8,a8)
      do i=1,Lch
         do j=1,Lch
            if(vvv(i,j).eq.-2)then
               write(20,*)i,j,'  not be checked for excluded volumn'
            endif
         enddo
      enddo
      write(20,*)'bond_max=',bond_max

*******for calculation of RMSD of fragment------------------------>
      do i=1,Lch
         i_chunk(i)=-1
      enddo
      do i=1,nfr
         do j=nfr_i(i),nfr_f(i)
            i_chunk(j)=i
            ex0(j)=ex(j)
            ey0(j)=ey(j)
            ez0(j)=ez(j)
         enddo
      enddo

c^^^^^^^^^^^^^^^^^^^^^^^^^ template_initial2 finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cSSSccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine secondary(Lch,cx0,cy0,cz0,q,L_cut)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)

      dimension cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)


CM      common/lengths/Lch,Lch1,Lch2
CM      common/consensus1/cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)
CM      common/aminoacid/sequ(ndim)
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)



      character*3 sequ
      dimension sec_high(0:ndim)
      dimension nse_i(200),nse_f(200),nse_type(200)

      dimension cx00(0:ndim),cy00(0:ndim),cz00(0:ndim)
      dimension ax(ndim),ay(ndim),az(ndim)

ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

ccc   read SSP --------------------------------->
      open(unit=23,file='seq.dat',status='old') !ssp from native
      do i=1,Lch
         read(23,707)k,name,sec_high(i)
c     write(*,*)i,k,name,sec_high(i)
      enddo
      close(23)
 707  format(i5,3x,a3,2i5)

ccc   find continuous SSP pieces -------------------->
      nse=0                     !number of fragments
      do i=1,200
         nse_i(i)=0
         nse_f(i)=0
      enddo
      sec_high(0)=1
      do i=1,Lch
         if(sec_high(i).ge.2.and.sec_high(i).ne.sec_high(i-1))then
            nse=nse+1
            nse_i(nse)=i        !initial location of nse'th ssp
            nse_type(nse)=sec_high(i) !type of ssp
         endif
         if(sec_high(i).ne.1)nse_f(nse)=i !final location of nse'th ssp
      enddo

ccc   remove small secondary structures -------------->
      nse_old=nse
      nse=0
      do i=1,nse_old
         len=nse_f(i)-nse_i(i)+1
         if(len.ge.L_cut)then
            nse=nse+1
            nse_i(nse)=nse_i(i)
            nse_f(nse)=nse_f(i)
            nse_type(nse)=nse_type(i)
         endif
      enddo

c      do i=1,nse
c         write(*,*)i,nse_i(i),nse_f(i),nse_f(i)-nse_i(i)+1,nse_type(i)
c      enddo

ccc   generate new (cx0,cy0,cz0) --------------------------->
      do i=1,Lch
         if(q(i).eq.1)then
            cx00(i)=cx0(i)
            cy00(i)=cy0(i)
            cz00(i)=cz0(i)
         endif
      enddo
      do 10 i=1,nse             !nse of SSP
         len=nse_f(i)-nse_i(i)+1 !length of SSP

ccc   generate standard fragments --------------------->
         if(nse_type(i).eq.2)then
            call alpha_helix(ax,ay,az,len) !alpha-helix
         else
            call beta_sheet(ax,ay,az,len) !beta-sheet
         endif
         
ccc   rotation matrix based on initial template ------->
*1: aligned
         n=0
         do j=nse_i(i),nse_f(i)
            if(q(j).eq.1)then
               n=n+1
               k=j-nse_i(i)+1
               r_1(1,n)=ax(k)
               r_1(2,n)=ay(k)
               r_1(3,n)=az(k)
               r_2(1,n)=cx00(j)
               r_2(2,n)=cy00(j)
               r_2(3,n)=cz00(j)
            endif
         enddo
*2: gapped at middle
         if(n.lt.2)then         !it is a gap at ssp region
            r_1(1,1)=ax(1)
            r_1(2,1)=ay(1)
            r_1(3,1)=az(1)
            r_1(1,2)=ax(nse_f(i)-nse_i(i)+1)
            r_1(2,2)=ay(nse_f(i)-nse_i(i)+1)
            r_1(3,2)=az(nse_f(i)-nse_i(i)+1)
            n=0
            j=nse_i(i)
 11         j=j-1
            if(q(j).ne.1.and.j.gt.1)goto 11
            if(q(j).eq.1)then
               n=n+1
               ax1=cx00(j)
               ay1=cy00(j)
               az1=cz00(j)
               j0=j
               mm=1
            endif
            j=nse_f(i)
 12         j=j+1
            if(q(j).ne.1.and.j.lt.Lch)goto 12
            if(q(j).eq.1)then
               n=n+1
               ax2=cx00(j)
               ay2=cy00(j)
               az2=cz00(j)
               mm=2
             endif
             if(n.eq.2)then
               aaa=sqrt(((ax2-ax1)**2)+((ay2-ay1)**2)+((az2-az1)**2))
               bbb=sqrt(((ax(len)-ax(1))**2)+((ay(len)-ay(1))**2)
     &                  +((az(len)-az(1))**2))
               ccc=abs(nse_i(i)-j0)
CM               ccc=abs(nse_i(i)-j0)*1
               r_2(1,1)=ax1+(ax2-ax1)/aaa*ccc
               r_2(2,1)=ay1+(ay2-ay1)/aaa*ccc
               r_2(3,1)=az1+(az2-az1)/aaa*ccc
               r_2(1,2)=ax1+(ax2-ax1)/aaa*(ccc+bbb)
               r_2(2,2)=ay1+(ay2-ay1)/aaa*(ccc+bbb)
               r_2(3,2)=az1+(az2-az1)/aaa*(ccc+bbb)
             endif
         endif
*3: gap at terminal
         if(n.lt.2)then         !n=1, the ssp gap is at terminal
            n=2
            if(mm.eq.1)then     !C-terminal
*3a: C-terminal
               m=0
               j=nse_i(i)
 13            j=j-1
               if(q(j).eq.1)m=m+1
               if(m.eq.1.and.q(j).eq.1)then
                  j0=j
                  ax2=cx00(j)
                  ay2=cy00(j)
                  az2=cz00(j)
               endif
               if(m.lt.5)goto 13 !the minimum alignment length is 5
               if(m.eq.5.and.q(j).eq.1)then
                  ax1=cx00(j)
                  ay1=cy00(j)
                  az1=cz00(j)
               endif
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(nse_i(i)-j0)*2
               r_2(1,1)=ax2+(ax2-ax1)/aaa*ccc
               r_2(2,1)=ay2+(ay2-ay1)/aaa*ccc
               r_2(3,1)=az2+(az2-az1)/aaa*ccc
               r_2(1,2)=ax2+(ax2-ax1)/aaa*(ccc+bbb)
               r_2(2,2)=ay2+(ay2-ay1)/aaa*(ccc+bbb)
               r_2(3,2)=az2+(az2-az1)/aaa*(ccc+bbb)
            else
*3a: N-terminal
               m=0
               j=nse_f(i)
 14            j=j+1
               if(q(j).eq.1)m=m+1
               if(m.eq.1.and.q(j).eq.1)then
                  j0=j
                  ax1=cx00(j)
                  ay1=cy00(j)
                  az1=cz00(j)
               endif
               if(m.lt.5)goto 14 !the minimum alignment length is 5
               if(m.eq.5.and.q(j).eq.1)then
                  ax2=cx00(j)
                  ay2=cy00(j)
                  az2=cz00(j)
               endif
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(j0-nse_f(i))*2
               r_2(1,2)=ax1+(ax1-ax2)/aaa*ccc
               r_2(2,2)=ay1+(ay1-ay2)/aaa*ccc
               r_2(3,2)=az1+(az1-az2)/aaa*ccc
               r_2(1,1)=ax1+(ax1-ax2)/aaa*(ccc+bbb)
               r_2(2,1)=ay1+(ay1-ay2)/aaa*(ccc+bbb)
               r_2(3,1)=az1+(az1-az2)/aaa*(ccc+bbb)
            endif
         endif
c         write(*,*)'n=',n,'  i=',i
c         do j=1,n
c            write(*,*)j,r_2(1,j),r_2(2,j),r_2(3,j)
c         enddo
         call u3b(w,r_1,r_2,n,1,rms,u,t,ier) !u rotate r_1 to r_2
ccc   rotate ssp onto initial template ------------>
         do j=nse_i(i),nse_f(i)
            k=j-nse_i(i)+1
            cx0(j)=t(1)+u(1,1)*ax(k)+u(1,2)*ay(k)+u(1,3)*az(k)
            cy0(j)=t(2)+u(2,1)*ax(k)+u(2,2)*ay(k)+u(2,3)*az(k)
            cz0(j)=t(3)+u(3,1)*ax(k)+u(3,2)*ay(k)+u(3,3)*az(k)
         enddo
 10   continue

ccc   re-define q(i) ------------------------------>
      do i=1,Lch
         q(i)=0
      enddo
      do i=1,nse
         do j=nse_i(i),nse_f(i)
            q(j)=1
         enddo
      enddo

      return
      end


ccccc generate a random walk of i2-i1 peptide bonds
c     contiguous peptide bonds not require to make a good angle. This
c     peptide bonds are off-lattice and in Angstrom units
      subroutine connect(i1,i2,pass,cx,cy,cz,Lch,amcheck_dis)
      implicit integer(i-z)
      parameter(nvec=312)


      dimension cx(1000),cy(1000),cz(1000)


CM      common/lengths/Lch,Lch1,Lch2
CM      common/initialstr/cx(1000),cy(1000),cz(1000)
CM      common/mcheck_dis/amcheck_dis

      print *,'In connect i1',i1,i2
      bdis0=3.5
c     minimum CA-CA distance to avoid steric problem
      amcheck_dis=3.1            
      pass=1
c     i1 is the N-terminal, begin from here
      if(i1.eq.1)then           

ccccc go residue by residue from i1 up to i2
         do i=i2,1,-1
c     number of times we try a peptide bond vector
            n_check=0
c     get a 3.8Angstrom peptide bond with random orientation
 10         call get_bond(axx,ayy,azz,3.8)
c     tentative coordinates
            cx(i)=cx(i+1)+axx
            cy(i)=cy(i+1)+ayy
            cz(i)=cz(i+1)+azz
            n_check=n_check+1
            if(n_check.gt.100000)then
c     we tried too many times, indicate failure in building the
c     unaligned fragment and exit subroutine
               pass=3
               return
            endif
c     check the tentative coordinates. If steric problem (.eq.3), try
c     other peptide bond
          if(mcheck(i,cx(i),cy(i),cz(i),cx,cy,cz,Lch,amcheck_dis).eq.3) 
     +      goto 10
c     go to next unaligned residue
         enddo

c     i2 is the C-terminal
      elseif(i2.eq.Lch)then     
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
            if(mcheck(i,cx(i),cy(i),cz(i),cx,cy,cz,Lch,amcheck_dis)
     +         .eq.3) goto 11
         enddo
c     i1 is not N-terminal and i2 is not C-terminal
      else                    
c     distance between i1-1 and i2+1
CM         print *,'cx cy cz',cx(i2+1),cy(i2+1),cz(i2+1),
CM     +            cx(i1-1),cy(i1-1),cz(i1-1)
         diss=di(cx(i2+1),cy(i2+1),cz(i2+1),cx(i1-1),cy(i1-1),cz(i1-1))
c     average distance per residue from going from i1-1 to i2+1
         adis=diss/float((i2+1)-(i1-1))
CM         print *,'adis',diss,adis
c     minimum peptide bond distance to try the single-line construction
c     (see below)
         adis0=3.5
         x_check=0
c     number of times we try the single-line construction
 15      x_check=x_check+1

ccccc single-line construction: connect unaligned residues along a
c     straight line if distance
c     between i1-1 and i2+1 is large enough
         if(adis.ge.adis0)then
            dex=3.5*(cx(i1-1)-cx(i2+1))/diss
            dey=3.5*(cy(i1-1)-cy(i2+1))/diss
            dez=3.5*(cz(i1-1)-cz(i2+1))/diss
            do j=i2,i1,-1
               cx(j)=cx(j+1)+dex
               cy(j)=cy(j+1)+dey
               cz(j)=cz(j+1)+dez
            enddo
c     finish building the random walk
            return
         endif

ccccc random walk from i1 to i2
         m_check=0             
c     number of times we try to build a random walk for the unaligned
c     residues
 13      m_check=m_check+1
         do 14 i=i1,i2
            n_check=0          
 12         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i-1)+axx
            cy(i)=cy(i-1)+ayy
            cz(i)=cz(i-1)+azz
            n_check=n_check+1
c     if we tried peptide bonds n*1000 times, then...
            if(int(n_check/1000)*1000.eq.n_check)then 
c     increase the maximum peptide bond distance
               bdis0=bdis0*1.03
               if(bdis0.ge.3.8)bdis0=3.8
c     minimum CA-CA distance to avoid steric problem
               amcheck_dis=2.8               
            endif
            if(n_check.eq.4000)then
c     we tried this peptide bond too many times, which may mean other
c     previously built peptide bonds are interferring. Thus trash them
c     and start afresh from i1
               goto 13
            endif
c     if we tried more than 2000 times to build a random walk, then...
            if(m_check.ge.2000)then
c     if we tried the single-line construction less than six times,
c     then...
               if(x_check.le.6)then
c     try the single-line construction with a smaller peptide bond
c     distance paramenter
                  adis0=adis0*0.995
                  goto 15
               endif
               pass=5
               return
            endif
c     check for the excluded volume of CA's
            if(mcheck(i,cx(i),cy(i),cz(i),cx,cy,cz,Lch,amcheck_dis)
     +         .eq.3) goto 12
            aaa=float(i2+1-i)
c     average peptide bond distance from i to i2+1
            bdis=di(cx(i),cy(i),cz(i),cx(i2+1),cy(i2+1),cz(i2+1))/aaa
c     we require that the remaining residues are not too extended
            if(i.lt.i2.and.bdis.ge.bdis0) goto 12
            if(i.eq.i2)then
c     for the last residue, we require peptide bond joining i2 and i2+1
c     has a distance in [3.4, 4.2]
               if(bdis.gt.4.2.or.bdis.lt.3.4) goto 12 !last step
            endif
c     reset bdis0
            bdis0=3.5
c     reset amcheck_dis
            amcheck_dis=3.1
 14      continue
      endif
CM      print *,'Out connect amcheck_dis',amcheck_dis
      return
      end
ccccc connect finished


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sort ras, recalculate nfl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sort_ras(nfl,ras)
      implicit integer(i-z)
      parameter(ndim=1500)
CM      parameter(nrep=100)
CM      parameter(nvec=312)

      dimension ras(ndim)

CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain0/ras(ndim),nfl

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


ccccc decide movement point from ras(i)
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
      subroutine move_point(Lch,mv,nfl,ras,Mend_N,Mend,Mend_C,
     +                      nfl6,ras6,nfl5,ras5,nfl4,ras4,nfl3,
     +                      ras3,nfl2,ras2)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)

      dimension ras(ndim)
      dimension ras2(ndim)
      dimension ras3(ndim)
      dimension ras4(ndim)
      dimension ras5(ndim)
      dimension ras6(ndim)
      dimension mv(ndim)

CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain0/ras(ndim),nfl
CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/chainm/mv(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

      print *,'In move_point',Lch,nfl

c     mv(i).eq.1 if position "i" is moveable. We extend over Lch because
c     of the algorithm employed to find nfl2, nfl3,...
      do i=1,Lch+10
         mv(i)=-1               
      enddo

c     go through list of moveable points and set mv(i)
      do i=1,nfl
         mv(ras(i))=1
      enddo

ccccc beginning from the N-terminal, find number of moveable positions
c     (not to exceed Mend)
      Mend_N=0            
      k=0
c     go through every amino acid position
      do i=1,Lch
         k=k+1
c     if we reached the maximum number of moveable positions (Mend) or
c     we reached a frozen fragment (mv(i).eq.-1), then stop going
c     through every amino acid position
         if(k.gt.Mend.or.mv(i).lt.0) goto 111
         Mend_N=Mend_N+1
      enddo
c     go to next position
 111  continue

cccc  beginning from the C-terminal, find number of moveable positions
c     (not to exceed Mend)
      Mend_C=0            
      k=0
      do i=Lch,1,-1
         k=k+1
         if(k.gt.Mend.or.mv(i).lt.0) goto 222
         Mend_C=Mend_C+1
      enddo
 222  continue

ccccc find 6-bond-move position
c     nfl6 total number of 6-bond-move points
      nfl6=0                    
      do i=1,Lch               
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then
                     if(mv(i+4).gt.0)then
                        if(mv(i+5).gt.0)then
                           if(mv(i+6).gt.0)then
                              nfl6=nfl6+1
c     ras6(j) denotes amino acid position of the beginning of a moveable
c     chunk of six amino acids
                              ras6(nfl6)=i
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

ccccc find 5-bond-move position
      nfl5=0                   
      do i=1,Lch               
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then
                     if(mv(i+4).gt.0)then
                        if(mv(i+5).gt.0)then 
                           nfl5=nfl5+1
                           ras5(nfl5)=i
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

ccccc find 4-bond-move position
      nfl4=0
      do i=1,Lch
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then
                     if(mv(i+4).gt.0)then
                        nfl4=nfl4+1
                        ras4(nfl4)=i
                     endif
                  endif
               endif
            endif
         endif
      enddo

ccccc find 3-bond-move position
      nfl3=0                   
      do i=1,Lch               
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  if(mv(i+3).gt.0)then
                     nfl3=nfl3+1
                     ras3(nfl3)=i
                  endif
               endif
            endif
         endif
      enddo

ccccc find 2-bond-move position
      nfl2=0
      do i=1,Lch
         if(mv(i).gt.0)then
            if(mv(i+1).gt.0)then
               if(mv(i+2).gt.0)then
                  nfl2=nfl2+1
                  ras2(nfl2)=i
               endif
            endif
         endif
      enddo

      print *,'Mend_N,Mend,Mend_C',Mend_N,Mend,Mend_C
      print *,'nfl ras etc',nfl6,ras6(nfl6),nfl5,ras5(nfl5),
     +         nfl4,ras4(nfl4),nfl3,ras3(nfl3),nfl2,ras2(nfl2)
      return
      end
ccccc move_point finished



ccccc decide pairs for which excluded-volume should not be checked:
c     1, pairs in same frozen fragments;
c     2, pairs in different frozen fragments but initially violated.
      subroutine get_vvv(Lch,vvv,mv,exc)
      implicit integer (i-z)
      parameter(ndim=1500)

      dimension vvv(ndim,ndim)
      dimension mv(ndim)

CM      common/excluded/vvv(ndim,ndim)
CM      common/chainm/mv(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/looks/exc,exc1,exc2

      print *,'In get_vvv',Lch,exc

      do i=1,Lch
         do j=1,Lch
ccccc initialize to 1 means check for excluded-volume violation for pair
c     (i,j)
            vvv(i,j)=1
c     if residues are frozen, then..
            if(mv(i).lt.0.and.mv(j).lt.0)then
c     if residues belong to the same fragment, then..
               if(mv(i).eq.mv(j))then
c     do not check (i,j) for exclude-volume violation
                  vvv(i,j)=-1   
               else
                  dist2=(aax(i,mv,x,ex)-aax(j,mv,x,ex))**2+
     +                  (aay(i,mv,y,ey)-aay(j,mv,y,ey))**2
     $                 +(aaz(i,mv,z,ez)-aaz(j,mv,z,ez))**2
                  if(dist2.lt.exc)then
c     (i,j) were already in steric hindrance in the template, thus don't
c     check for excluded-volume violation
                     vvv(i,j)=-2
                  endif
               endif
            endif
CM            print *,'vvv',i,j,vvv(i,j)
         enddo
      enddo
      return
      end
ccccc get_vvv finished


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     generate standard alpha-helix:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine alpha_helix(x,y,z,n)
      dimension x(1000),y(1000),z(1000)

      rad=2.3                   !redius of helix
      do i=1,n
         angle=100*3.1415926/180*(i-1) !100 degree per residues
         x(i)=rad*cos(angle)
         y(i)=rad*sin(angle)
         z(i)=1.5*(i-1)            !increase 1.5 per residue
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     generate standard beta-sheet:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine beta_sheet(x,y,z,n)
      dimension x(1000),y(1000),z(1000)

      dist2=6.70098             !distance between i and i+2
      high=1.792820             !distance between i+1 and middle of (i,i+2)
      do i=1,n
         x(i)=dist2/2*(i-1)
         if(int(i/2)*2.eq.i)then
            y(i)=0
         else
            y(i)=high
         endif
         z(i)=0
      enddo

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


ccccc Calculate eigenvalue and normalized eigenvectors
c     A{3,3}-matrix E{3}-eigenvalue T{3,3}-eigenvector
c     E(1)=<x^2>; E(1)=<y^2>; E(1)=<z^2> in the rotated system
      subroutine eigenvalue(A,E,T)
      real A(3,3),E(3),T(3,3)

CM      common/eigen/A(3,3),E(3),T(3,3)

      pi=3.1415926
      p1=-1
      p2=A(1,1)+A(2,2)+A(3,3) 
      p3=-A(1,1)*A(2,2)-A(1,1)*A(3,3)-A(2,2)*A(3,3)+
     &     A(2,3)*A(2,3)+A(1,3)*A(1,3)+A(1,2)*A(1,2)
      p4=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+
     &     A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)-
     &     A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)
      p5=(-(1.0/3)*(p2/p1)**2+p3/p1)/3 
      ap5=sqrt(-p5)
      p6=((2.0/27)*(p2/p1)**3-(1.0/3)*(p2*p3/p1**2)+p4/p1)/2
      p7=acos(-p6/sqrt(-p5**3))
      p8=2*ap5*cos(p7/3.0)
      p9=-2*ap5*cos((p7+pi)/3.0)
      p10=-2*ap5*cos((p7-pi)/3.0) 
      p11=p2/(3*p1)
c     Eigenvalues
      E(1)=p8-p11  
      E(2)=p9-p11  
      E(3)=p10-p11
      
ccccc normalized eigenvectors
      do i=1,3
         fnorm1=A(2,1)*A(1,2)-(A(1,1)-E(i))*(A(2,2)-E(i))
         x=((A(2,2)-E(i))*A(1,3)-A(1,2)*A(2,3))/fnorm1
         y=((A(1,1)-E(i))*A(2,3)-A(2,1)*A(1,3))/fnorm1
         T(i,3)=1/sqrt(x*x+y*y+1)
         T(i,1)=x*T(i,3)
         T(i,2)=y*T(i,3)
      enddo
      return
      end
ccccc eigenvalue finished



cccc produce a random vector of length=3.8A
c     dW=d(phi)*d(cos(athita)), where d() is the differential operator
      subroutine get_bond(axx,ayy,azz,al)
      implicit integer(i-z)
c     thita angle in random, [0,pi]
      athita=acos(1.-2.*aranzy())
CMM                  print *,'NO=',no
c     phi angle in random, [0,2pi]      
      aphi=2.*3.1415926*aranzy()
CMM                  print *,'NO=',no
      axx=al*sin(athita)*cos(aphi)
      ayy=al*sin(athita)*sin(aphi)
      azz=al*cos(athita)
CM      print *,'In get_bond',axx,ayy,azz
      return
      end
ccccc get_bond finished


ccccc check weak excluded volumn for cx(). Note that (cx(i),cy(i),cz(i))
ccccc are stupid (cx(i)=10000,cy(i)=10000,cz(i)=10000) for residues
ccccc belonging to unaligned fragments not yet built. Thus these
ccccc residues will not pose a steric problem
      function mcheck(i0,bx0,by0,bz0,cx,cy,cz,
     +                Lch,amcheck_dis)
      implicit integer(i-z)
      parameter(nvec=312)

      dimension cx(1000),cy(1000),cz(1000)


CM      common/lengths/Lch,Lch1,Lch2
CM      common/initialstr/cx(1000),cy(1000),cz(1000)
CM      common/mcheck_dis/amcheck_dis

      mcheck=1
      do i=1,Lch
         if(i0.ne.i)then
            dis=di(bx0,by0,bz0,cx(i),cy(i),cz(i))
            if(dis.le.amcheck_dis)then
c     steric problem
               mcheck=3
               return
            endif
         endif
      enddo

CM       print *,'In mcheck',mcheck,amcheck_dis

      return
      end
ccccc mcheck finished




ccccc set current (x,y,z), ica, T, from itemp's replica
      subroutine set_current(x,y,z,xrep,yrep,zrep,ica,vector,Lch,itemp,
     +                       Lch1)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)

      dimension xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension vector(-5:5,-5:5,-5:5)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)



CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/sw1/aT_rep(nrep),E_rep(nrep)
CM      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
CM      common/sw3/icarep(ndim,nrep)
CM      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
CM      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
CM      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/chain0/ras(ndim),nfl
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim)
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)

      print *,'In set_current',Lch,Lch1

ccccc put the itemp'th replica conformation as current conformation
      do i=1,Lch
         x(i)=xrep(i,itemp)
         y(i)=yrep(i,itemp)
         z(i)=zrep(i,itemp)
      enddo

ccccc ica(is) stores peptide bond indexes
      do i=1,Lch1
         j=i+1
         wx=x(j)-x(i)
         wy=y(j)-y(i)
         wz=z(j)-z(i)
         ica(i)=vector(wx,wy,wz)
CM         print *,'ica',i,vector(wx,wy,wz)
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
CM      print *,'ica',ica(0),ica(Lch)
      return
      end
ccccc set_current finished


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set current (x,y,z), ica, T, from itemp's replica
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_current_RS(Lch,x,y,z,xrep,yrep,zrep,ica,icarep,
     +                          ex,ey,ez,exrep,eyrep,ezrep,egx,egy,
     +                          egz,egxrep,egyrep,egzrep,ecx,ecy,ecz,
     +                          ecxrep,ecyrep,eczrep,ebx,eby,ebz,
     +                          ebxrep,ebyrep,ebzrep,
     +                          itemp,mv,Lch2)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)

      dimension icarep(ndim,nrep)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension mv(ndim)
      dimension xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      dimension egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      dimension ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      dimension ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)

CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

CM      common/sw1/aT_rep(nrep),E_rep(nrep)
CM      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
CM      common/sw3/icarep(ndim,nrep)
CM      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
CM      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) !SG
CM      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) !cc
CM      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) !hb

CM      common/chain0/ras(ndim),nfl
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)

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
ccccc set_current_RS finished


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  2-bond movement: prefabricated
c  residues of 'm, m1=m+1,m2=m+2' are involved,
c  but positions of m1 are changed.
c  old ---> new ---> old ----> new
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move2(nfl2,ras2,Np2,ica,v21,goodc,v22,x,y,z,
     +                      vx,vy,vz,Lch,nop,nopp,nom,nomm,noa,noaa,
     +                      afsn,afs,seq,envir,dord,en1,eprofo,
     +                      bNa,bNNa,itemp,icnto,icnt,sumcto,sumct,
     +                      codevsum,codev,didevsum,didev,energ,bNt,
     +                      bNNt,jbin,es2,er1,er3,er4,es3,es3a,es3b,
     +                      es3c,es4,es5,es6,er5,er6,er7,conew,colim,
     +                      coold,dinew,diold,eprofn,ehb1,ehb1a,ehb1b,
     +                      ehb1c,ehb2,ehb3,ehb4,ehb5a,enb5b,mv,gx,gy,
     +                      gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +                      egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +                      apar,alp,arlp,app,ala,arla,apa,alm,arlm,apm,
     +                      sec,concut2,acut_cc,acut_bb,acut_vv,bcut_cc,
     +                      bcut_bb,bcut_vv,en2,en3,contt,Cr2a,Cr2b,eh1,
     +                      eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +                      acut_hh,bcut_hh,EHBIJ,hh,ee,area,mdis,n_dis,
     +                      kdis,dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                      McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                      frga,frgb,ESHORT2,ESHORT3,ESHORT4,ESHORT4a,
     +                      ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c,ESHORT6,
     +                      ESHORT7,ESHORT8,ESHORT9,ESHORT10,ESHORT11,
     +                      cex,cey,cez,acrit,distL,aweig,dist_CA_cut,
     +                      dis,switch,AA,exc,vvv,atemp,ichos,aTs) 

      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1500)
      parameter(nvec=312)

      dimension ras2(ndim)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      logical look,goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension afs(ndim),afsn(ndim)
      real envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim),sec(ndim)
      dimension bNa(nrep),bNt(nrep),bNNa(nrep),bNNt(nrep)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension bsr(ndim,16), acops(ndim,16)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension ehbij(ndim,ndim)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2),asr(ndim,-12:12)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension AA(3,3),EE(3),HH(3,3)
      dimension vvv(ndim,ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
CM      real eprofn
      integer id


CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      double precision bNa,bNt,bNNa,bNNt
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)


CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6

      print *,'In move2',nfl2,itemp
     
c     choose the position of (m) ----------->
      i=int(aranzy()*nfl2)+1  ![1,nfl2]
CM      print *,'i',i
      m=ras2(i)
      m1=m+1
      m2=m+2

CM      print *,'i m',i,m,m1,m2

cccc  all the pairs have at least one another pair, so nc>=2
      nc=Np2(ica(m),ica(m1))    !number of pathes from m to m2

c     choose p'th new path from (m) to (m2) -------------->
 201  p=int(aranzy()*nc)+1    ![1,nc]
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

CM      if(look(m,m2))then        ! check excluded volumn for passage of [m,m2]
      if(LOOK(m,m2,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,m2,1)+ESHORT(m,m2,1) !icnt,nop are repeated.

CM      do igx=1,nvec
CM      do jgx=1,nvec
CM      do kgx=0,19
CM         print *,'gx=',igx,jgx,kgx,gx(igx,jgx,kgx)
CM      enddo
CM      enddo
CM      enddo

CM      CALL EXIT(-1)
      Enew=ehb(m,m2,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m2,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin) !icnt,nop are repeated

        print *,'move2 Enew=',Enew
         do kkk=m,m2
            afsn(kkk)=afs(kkk)  !afs(i) is useful in ESHORT
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
CM         Eold=EHB(m,m2,-1)+ESHORT(m,m2,-1) !repeated part of icnt, nop removed

         Eold=ehb(m,m2,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m2,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)
!repeated part of icnt, nop removed

         print *,'move2 Eold=',Eold

c     calculate eprofn while dord was calculated when call EHB(m,m2,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)          
            ia=noa(pp)          
            ip=nop(pp)          !nopp(old)+1(new)-1(old)=new
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

       print *,'eprofn=',eprofn

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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

      print *,'Out Move2'
c ^^^^^^^^^^^^^^^^^ move2 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3-bond movement: single-move
c  residues of 'i,i+1,i+2,i+3' are involved,
c  but positions bond-vectors of i+1, i+2 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move3s(nfl3,ras3,Np2,ica,v21,goodc,v22,x,y,z,vx,
     +                  vy,vz,nop,nopp,nom,nomm,noa,noaa,afsn,
     +                  afs,seq,envir,eprofo,bNa,bNNa,
     +                  icnto,icnt,sumcto,sumct,Lch,eprofn,
     +                  codevsum,codev,didevsum,didev,bNt,bNNt,
     +                  itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 

      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1500)
      parameter(nvec=312)

      dimension ras3(ndim)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension afs(ndim),afsn(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension seq(ndim),sec(ndim)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension ehbij(ndim,ndim)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension AA(3,3),EE(3),HH(3,3)


CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy()*nfl3)+1
      m=ras3(i)
      m1=m+1
      m2=m+2
      m3=m+3

ccccccccccccc 1th 2-bond movement cccccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccc 2th 2-bond movement cccccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy()*nc)+1    ![1,nc]
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

CM      if(look(m,m3))then        ! check excluded volumn for passage of [i,m]
      if(LOOK(m,m3,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,m3,1)+ESHORT(m,m3,1) !use ifs

        Enew=ehb(m,m3,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m3,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)
!use ifs

       print *,'In move3s enew=',Enew

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
CM         Eold=EHB(m,m3,-1)+ESHORT(m,m3,-1)

         Eold=ehb(m,m3,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m3,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move3s Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center       !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move3d(nfl3,ras3,Np2,ica,v21,goodc,v22,x,y,z,vx,vy,
     +                  vz,Lch,nop,nopp,nom,nomm,noa,noaa,afsn,afs,
     +                   seq,envir,energ,bNa,bNNa,icnto,icnt,sumct,
     +                   sumcto,eprofo,eprofn,codevsum,codev,didevsum,
     +                   didev,bNt,bNNt,itemp,
     +                 jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 

      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1500)
      parameter(nvec=312)

      dimension ras3(ndim)
      logical look,goodc(nvec,nvec)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension afs(ndim),afsn(ndim)
      dimension seq(ndim),sec(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension ehbij(ndim,ndim)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension AA(3,3),EE(3),HH(3,3)


CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy()*nfl3)+1
      m=ras3(i)
      m1=m+1
      m2=m+2
      m3=m+3

ccccccccccccc temporal 1th 2-bond movement cccccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy()*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccc temporal 2th 2-bond movement cccccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy()*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202

ccccccccccccc 1th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 203  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 203
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,t2))goto 203
ccccccccccccc 2th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,t2)
 204  p=int(aranzy()*nc)+1    ![1,nc]
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

CM      if(look(m,m3))then       ! check excluded volumn for passage of [i,m]
      if(LOOK(m,m3,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,m3,1)+ESHORT(m,m3,1) !use ifs

         Enew=ehb(m,m3,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m2,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move3d Enew=',Enew

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
CM         Eold=EHB(m,m3,-1)+ESHORT(m,m3,-1)

         Eold=ehb(m,m3,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m3,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move3d Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center       !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move4s(nfl4,ras4,Np2,ica,v21,v22,goodc,x,y,z,vx,vy,
     +                  vz,Lch,nop,nopp,nom,nomm,noa,noaa,afsn,
     +                  afs,eprofn,seq,envir,energ,bNa,bNNa,icnto,
     +                  icnt,sumcto,sumct,eprofo,codevsum,codev,
     +                  didevsum,didev,bNt,bNNt,itemp,     
     +                 jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 

      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1500)
      parameter(nvec=312)

      dimension ras4(ndim)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension afs(ndim),afsn(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim),sec(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension AA(3,3),EE(3),HH(3,3)
      dimension ehbij(ndim,ndim)




CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy()*nfl4)+1
      m=ras4(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy()*nc)+1    ![1,nc]
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

CM      if(look(m,m4))then       ! check excluded volumn for passage of [i,m]
      if(LOOK(m,m4,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,m4,1)+ESHORT(m,m4,1) !use ifs

         Enew=ehb(m,m4,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m4,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move4s Enew=',Enew

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
CM         Eold=EHB(m,m4,-1)+ESHORT(m,m4,-1)

         Eold=ehb(m,m4,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m4,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move4s Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center       !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move4d(nfl4,ras4,Np2,ica,v21,v22,goodc,x,y,z,vx,
     +                  vy,vz,Lch,nop,nopp,nom,nomm,noa,noaa,afsn,afs,
     +                  eprofn,seq,envir,bNa,bNNa,icnto,icnt,
     +                  sumcto,sumct,eprofo,codevsum,codev,didevsum,
     +                  didev,bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 

      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1500)
      parameter(nvec=312)


      dimension ras4(ndim)
      logical look, goodc(nvec,nvec)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(nrep),bNt(nrep),bNNa(nrep),bNNt(nrep)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      dimension JBIN(0:500)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)




CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy()*nfl4)+1
      m=ras4(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4

ccccccccccccccc 1th temperor 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy()*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy()*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy()*nc)+1    ![1,nc]
      tt2=v21(t2,ica(m3),p)
      if(.not.goodc(tt1,tt2))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203

ccccccccccccc first 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 204  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 204
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,tt2))goto 204
ccccccccccccc second 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,tt2)
 205  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,tt2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 205
      ttt2=v22(ttt1,tt2,p)
      if(.not.goodc(ttt2,t3))goto 205
ccccccccccccc third 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt2,t3)
 206  p=int(aranzy()*nc)+1    ![1,nc]
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

CM      if(look(m,m4))then       ! check excluded volumn for passage of [i,m]
      if(LOOK(m,m4,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,m4,1)+ESHORT(m,m4,1) !use ifs

         Enew=ehb(m,m4,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m4,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move4d Enew=',Enew

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
CM         Eold=EHB(m,m4,-1)+ESHORT(m,m4,-1)

         Eold=ehb(m,m4,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m4,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move4d Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center       !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move5s(nfl5,ras5,Np2,ica,v21,goodc,v22,x,y,z,vx,vy,
     +                  vz,Lch,nop,nopp,nom,nomm,noa,noaa,afsn,afs,
     +                  eprofn,envir,bNa,bNNa,icnto,icnt,sumcto,
     +                  sumct,eprofo,codevsum,codev,didevsum,didev,
     +                  bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,seq,AA,atemp,ichos,aTs) 
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1500)
      parameter(nvec=312)

      logical look, goodc(nvec,nvec)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      dimension ras5(ndim)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension ehbij(ndim,ndim)
      dimension AA(3,3),EE(3),HH(3,3)






CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM
c     choose the position of (m) ----------->
      i=int(aranzy()*nfl5)+1
      m=ras5(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4
      m5=m+5

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m2)=v21(t2,ica(m3),p)
      if(.not.goodc(nn(m1),nn(m2)))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203
ccccccccccccccc 4th 2-bond movement cccccccccccccccc
      nc=Np2(t3,ica(m4))
 204  p=int(aranzy()*nc)+1    ![1,nc]
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

CM      if(look(m,m5))then       ! check excluded volumn for passage of [i,m]
      if(LOOK(m,m5,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv))then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,m5,1)+ESHORT(m,m5,1) !use ifs

         Enew=ehb(m,m5,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m5,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move5s Enew=',Enew

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
CM         Eold=EHB(m,m5,-1)+ESHORT(m,m5,-1)

         Eold=ehb(m,m5,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m5,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move5s Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center       !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move5d(nfl5,ras5,Np2,ica,v21,goodc,v22,x,y,z,vx,vy,vz,
     +                  Lch,nop,nopp,nom,nomm,noa,noaa,afsn,afs,eprofn,
     +                  seq,envir,bNa,bNNa,icnto,icnt,sumcto,sumct,
     +                  eprofo,codevsum,codev,didevsum,didev,bNt,bNNt,
     +                  itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1500)
      parameter(nvec=312)


      dimension ras5(ndim)
      logical look, goodc(nvec,nvec)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      dimension NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension ehbij(ndim,ndim)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension AA(3,3),EE(3),HH(3,3)



CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6

c     choose the position of (m) ----------->
      i=int(aranzy()*nfl5)+1
      m=ras5(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4
      m5=m+5

ccccccccccccccc temperor 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy()*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc temperor 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy()*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc temperor 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy()*nc)+1    ![1,nc]
      tt2=v21(t2,ica(m3),p)
      if(.not.goodc(tt1,tt2))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203
ccccccccccccccc temperor 4th 2-bond movement cccccccccccccccc
      nc=Np2(t3,ica(m4))
 204  p=int(aranzy()*nc)+1    ![1,nc]
      tt3=v21(t3,ica(m4),p)
      if(.not.goodc(tt2,tt3))goto 204
      t4=v22(t3,ica(m4),p)
      if(.not.goodc(t4,ica(m5)))goto 204

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(t,tt1)
 205  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 205
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,tt2))goto 205
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(ttt1,tt2)
 206  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,tt2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 206
      ttt2=v22(ttt1,tt2,p)
      if(.not.goodc(ttt2,tt3))goto 206
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(ttt2,tt3)
 207  p=int(aranzy()*nc)+1    ![1,nc]
      nn(m2)=v21(ttt2,tt3,p)
      if(.not.goodc(nn(m1),nn(m2)))goto 207
      ttt3=v22(ttt2,tt3,p)
      if(.not.goodc(ttt3,t4))goto 207
ccccccccccccccc 4th 2-bond movement cccccccccccccccc
      nc=Np2(ttt3,t4)
 208  p=int(aranzy()*nc)+1    ![1,nc]
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

CM      if(look(m,m5))then       ! check excluded volumn for passage of [i,m]
      if(LOOK(m,m5,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv))then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,m5,1)+ESHORT(m,m5,1) !use ifs

         Enew=ehb(m,m5,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m5,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move5d Enew=',Enew

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
CM         Eold=EHB(m,m5,-1)+ESHORT(m,m5,-1)

         Eold=ehb(m,m5,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,m5,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move5d Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center       !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move6(nfl6,ras6,vx,vy,vz,ica,vector,goodc,x,y,z,Lch,
     +                 nop,nopp,nom,nomm,noa,noaa,afsn,afs,seq,eprofn,
     +                 envir,eprofo,bNa,bNNa,icnto,icnt,sumcto,sumct,
     +                 codevsum,codev,didevsum,didev,bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)


      logical look,goodc(nvec,nvec)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension NOPP(ndim),NOMM(ndim),NOAA(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(nrep),bNt(nrep),bNNa(nrep),bNNt(nrep)
      dimension ras6(ndim)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension ehbij(ndim,ndim)
      dimension AA(3,3),EE(3),HH(3,3)


CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6

c      nob=6+int(aranzy()*6.999) !number of moved bonds, [6,12]
      nob=6

      print *,'In move6',nfl6
c     choose the position (m) to be moved
      i=int(aranzy()*nfl6)+1 ![1,nfl6]
      m=ras6(i)
      m1=m+1
      mnob=m+nob
      mnob1=mnob-1
      mnob2=mnob-2

 201  mx=int(aranzy()*2.999999)-1 ![-1,0,1]
      my=int(aranzy()*2.999999)-1 ![-1,0,1]
      mz=int(aranzy()*2.999999)-1 ![-1,0,1]
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

CM      if(look(m,mnob))then      ! check excluded volumn for passage of [i,m]
      if(LOOK(m,mnob,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv))then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,mnob,1)+ESHORT(m,mnob,1) !use ifs

         Enew=ehb(m,mnob,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,mnob,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move6 Enew=',Enew

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
CM         Eold=EHB(m,mnob,-1)+ESHORT(m,mnob,-1)

         Eold=ehb(m,mnob,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,mnob,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move6 Eold=',Eold
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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move_n_end(Mend_N,x,y,z,ica,anvec,goodc,vx,vy,vz,Lch,
     +                      nop,nopp,nom,nomm,noa,noaa,afsn,afs,seq,
     +                      eprofn,envir,bNa,bNNa,icnto,icnt,sumcto,
     +                      sumct,eprofo,codevsum,codev,didevsum,didev,
     +                      bNt,bNNt,itemp,
     +                 energ,jbin, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,
     +                coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +                hh,ee,area,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                switch,AA,atemp,ichos,aTs) 
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)

      logical look, goodc(nvec,nvec)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension ehbij(ndim,ndim)
      dimension AA(3,3),EE(3),HH(3,3)


CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      double precision bNa,bNt,bNNa,bNNt
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

      m=int(aranzy()*(Mend_N-1))+1 !m=[1,Mend_N-1], [1,m] will be moved.

c     back-up old path---------->
      do i=1,m
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path and move to new path -------------->
      do i=m,1,-1
 111     iv=int(aranzy()*anvec)+1
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

CM      if(look(1,m+1))then       ! check excluded volumn for passage of [2,m+1]
      if(LOOK(1,m+1,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(1,m+1,1)+ESHORT(1,m+1,1) !use ifs

         Enew=ehb(1,m+1,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,m+1,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move_n_end Enew=',Enew

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
CM         Eold=EHB(1,m+1,-1)+ESHORT(1,m+1,-1) !cut off old nop by istat=-1

         Eold=ehb(1,m+1,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,m+1,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)
!cut off old nop by istat=-1

       print *,'In move_n_end Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
      subroutine move_c_end(Mend_C,Lch,x,y,z,ica,anvec,goodc,vx,
     +             vy,vz,nop,nopp,nom,nomm,noa,noaa,afsn,afs,
     +             eprofn,seq,envir,bNa,bNNa,icnto,icnt,
     +             sumcto,sumct,eprofo,codevsum,codev,didevsum,
     +             didev,bNt,bNNt,itemp,
     +             energ,jbin, 
     +             es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +             er5,er6,er7,conew,colim,
     +             coold,dinew,diold,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,
     +             gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,
     +             hh,ee,area,mdis,n_dis,kdis,
     +             dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +             McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +             frga,frgb, 
     +             ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +             ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +             ESHORT10,ESHORT11,cex,
     +             cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +             switch,AA,Lch2,exc,vvv,atemp,ichos,aTs) 
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)

      logical look, goodc(nvec,nvec)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim),sec(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2)
      real asr(ndim,-12:12)
      dimension bsr(ndim,16), acops(ndim,16)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension ehbij(ndim,ndim)
      dimension vvv(ndim,ndim)


CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      double precision bNa,bNt,bNNa,bNNt
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

      print *,'In move_c_end',Lch,Mend_C,anvec,Lch2

      m=int(aranzy()*(Mend_C-1))+1 !m=[1,Mend_C-1], m points will be moved
      mm=Lch-m+1                !(mm,Lch] will be relocated
c     mm-1 should be moveable residue, ica(m-2) exists.

      print *,'mm=',mm
c     backup old path----------->
      do i=mm,Lch
         ox(i)=x(i)             !memory of old conformation
         oy(i)=y(i)             !memory of old conformation
         oz(i)=z(i)             !memory of old conformation
         oo(i-1)=ica(i-1)       !memory of old conformation
      enddo
c     prepare new path------------->

      do i=mm,Lch
CM         print *,'NO=',no
 111        ff=aranzy()
         print *,'FF=',ff
            iv=int(ff*anvec)+1
CM 111     iv=int(aranzy()*anvec)+1
         print *,'in do iv=',i,iv,ica(i-2),goodc(ica(i-2),iv)
CM         print *,'NO=',no
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

CM      if(look(mm-1,Lch))then    ! check excluded volumn for passage of [i,m]
      if(LOOK(mm-1,Lch,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(mm-1,Lch,1)+ESHORT(mm-1,Lch,1) !use ifs

         Enew=ehb(mm-1,Lch,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(mm-1,Lch,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

       print *,'In move_c_end Enew=',Enew

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
CM         Eold=EHB(mm-1,Lch,-1)+ESHORT(mm-1,Lch,-1) !nop go to new

         Eold=ehb(mm-1,Lch,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(mm-1,Lch,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)


       print *,'In move_c_end Eold=',Eold

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

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
CM            call get_center       !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
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
c     translation (d_xyz0) + rotation (angle0)
c     Move the off-lattice N-fragment.
c     This movement include an off-lattice rotation of fragment
c     and a two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trot_N(nfr_f,d_xyz00,angle00,ex,ey,ez,x,y,z,anvec,ica,
     +                  vx,vy,vz,goodc,Nw,W21,W22,egx,egy,egz,ecx,ecy,
     +                  ecz,ebx,eby,ebz,Lch,nop,nopp,nom,nomm,noa,noaa,
     +                  afsn,afs,aprofn,seq,envir,id,bNa,bNNa,icnto,
     +                  icnt,sumcto,sumct,eprofo,codevsum,codev,
     +                  didevsum,didev,bNt,bNNt,itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)

      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)
      dimension d_xyz00(ndim),angle00(ndim)
      real A(3,3),E(3),T(3,3)





CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)

CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      n=nfr_f(ifr)+2            !ending point of movement
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy()*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      ax=ex(n2)
      ay=ey(n2)
      az=ez(n2)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check n-terminal of fragment--------------->
      ax1=ax0+ax+(ex(n2)-ax)*a11+(ey(n2)-ay)*a12+(ez(n2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(n2)-ax)*a21+(ey(n2)-ay)*a22+(ez(n2)-az)*a23
      az1=az0+az+(ex(n2)-ax)*a31+(ey(n2)-ay)*a32+(ez(n2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=1,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=1,n2
        ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
        ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
        ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
        egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
        egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
        egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
        ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
        ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
        ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
        ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
        eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
        ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
      enddo
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=1,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(1,n))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(1,n,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(1,n,1)+ESHORT(1,n,1) !use ifs

         Enew=ehb(1,n,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,n,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=1,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(1,n,-1)+ESHORT(1,n,-1)

         Eold=ehb(1,n,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,n,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)


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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(22)=bNa(22)+1
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
            do kkk=1,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=1,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(22)=bNt(22)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ trot_N finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0) + rotation (angle0)
c     Move C-terminal of off-lattice fragment
c     This movement include an off-lattice rotation of fragment
c     and 1 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trot_C(nfr_i,d_xyz00,angle00,ex,ey,ez,x,y,z,anvec,
     +                  goodc,vx,vy,vz,Nw,W21,W22,ica,Lch,
     +                  egx,egy,egz,ecx,ecy,ecz,ebx,eby,ebz,nop,nopp,
     +                  nom,nomm,noa,noaa,afsn,afs,eprofn,seq,envir,
     +                  id,bNa,bNNa,icnto,icnt,sumcto,sumct,eprofo,
     +                  codevsum,codev,didevsum,didev,bNt,bNNt,itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)

      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension d_xyz00(ndim),angle00(ndim)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_i(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      real A(3,3),E(3),T(3,3)



CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)

CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)

      m=nfr_i(ifr)-2            !starting point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy()*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      ax=ex(m2)
      ay=ey(m2)
      az=ez(m2)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202

c     backup old path ------------>
 203  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      do i=m2,Lch
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      do i=m2,Lch
        ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
        ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
        ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
        egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
        egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
        egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
        ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
        ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
        ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
        ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
        eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
        ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      do i=m2,Lch
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(m,Lch))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(m,Lch,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,Lch,1)+ESHORT(m,Lch,1) !use ifs

         Enew=ehb(m,Lch,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,Lch,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=m,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(m,Lch,-1)+ESHORT(m,Lch,-1)

         Eold=ehb(m,Lch,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,Lch,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(24)=bNa(24)+1
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
            do kkk=m,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            do i=m2,Lch
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(24)=bNt(24)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ trot_C finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Translation (d_xyz0) + rotation (angle).
c     This movement include an off-lattice rotation of fragment
c     and 2 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3, n3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trot_M(nfr_i,d_xyz00,angle00,ex,ey,ez,x,y,z,anvec,
     +                  goodc,vx,vy,vz,Nw,W21,W22,ica,egx,egy,egz,ecx,
     +                  ecy,ecz,ebx,eby,ebz,nop,nopp,nom,nomm,noa,
     +                  noaa,afsn,afs,eprofn,Lch,seq,envir,id,bNa,bNNa,
     +                  icnto,icnt,sumcto,sumct,eprofo,codevsum,codev,
     +                  didevsum,didev,bNt,bNNt,itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)
      parameter(apai=3.1415926)

      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      dimension NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension d_xyz00(ndim),angle00(ndim)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)

      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_i(ndim),nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      real A(3,3),E(3),T(3,3)




CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)

CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)

      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5                  !maximum translation coordinates
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy()*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      n_rot=m2+int((n2-m2+1)*aranzy()) ![m2,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202
***   check n-terminal of fragment--------------->
 203  ax1=ax0+ax+(ex(n2)-ax)*a11+(ey(n2)-ay)*a12+(ez(n2)-az)*a13
      ay1=ay0+ay+(ex(n2)-ax)*a21+(ey(n2)-ay)*a22+(ez(n2)-az)*a23
      az1=az0+az+(ex(n2)-ax)*a31+(ey(n2)-ay)*a32+(ez(n2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=m2,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=m2,n2
        ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
        ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
        ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
        egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
        egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
        egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
        ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
        ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
        ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
        ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
        eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
        ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=m2,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(m,n))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(m,n,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs

         Enew=ehb(m,n,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,n,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

         Eold=ehb(m,n,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,n,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(23)=bNa(23)+1
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
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=m2,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(23)=bNt(23)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ trot_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0)
c     Move the off-lattice N-fragment.
c     This movement include an off-lattice rotation of fragment
c     and a two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tran_N(nfr_f,d_xyz00,ex,ey,ez,x,y,z,anvec,goodc,vx,vy,
     +                  vz,Nw,W21,W22,ica,egx,egy,egz,ecx,ecy,ecz,ebx,
     +                  eby,ebz,nop,nopp,nom,nomm,noa,noaa,afsn,afs,
     +                  Lch,seq,eprofn,envir,id,bNa,bNNa,icnto,icnt,
     +                  sumcto,sumct,eprofo,codevsum,codev,didevsum,
     +                  didev,bNt,bNNt,itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)


      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      dimension d_xyz00(ndim),angle00(ndim)
      real A(3,3),E(3),T(3,3)



CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)

CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      n=nfr_f(ifr)+2            !ending point of movement
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
***   check n-terminal of fragment--------------->
      ax1=ax0+ex(n2)
      ay1=ay0+ey(n2)
      az1=az0+ez(n2)
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=1,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=1,n2
         ex_n(i)=ax0+ex(i)
         ey_n(i)=ay0+ey(i)
         ez_n(i)=az0+ez(i)
         egx_n(i)=ax0+egx(i)
         egy_n(i)=ay0+egy(i)
         egz_n(i)=az0+egz(i)
      enddo
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=1,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(1,n))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(1,n,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(1,n,1)+ESHORT(1,n,1) !use ifs

         Enew=ehb(1,n,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,n,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=1,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(1,n,-1)+ESHORT(1,n,-1)

           Eold=ehb(1,n,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,n,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(16)=bNa(16)+1
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
            do kkk=1,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=1,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(16)=bNt(16)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ tran_N finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0)
c     Move C-terminal of off-lattice fragment
c     This movement include an off-lattice rotation of fragment
c     and 1 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tran_C(nfr_i,ifr,d_xyz00,ex,ey,ez,x,y,z,anvec,ica,
     +                  goodc,vx,vy,vz,Nw,egx,egy,egz,ecx,ecy,ecz,ebx,
     +                  eby,ebz,Lch,nop,nopp,nom,nomm,noa,noaa,afsn,afs,
     +                  seq,eprofn,envir,id,bNa,bNNa,icnto,icnt,sumcto,
     +                  sumct,eprofo,codevsum,codev,didevsum,didev,
     +                  bNt,bNNt,itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)


      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension d_xyz00(ndim),angle00(ndim)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_i(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      real A(3,3),E(3),T(3,3)




CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)

CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)

      m=nfr_i(ifr)-2            !starting point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
***   check m-terminal of fragment--------------->
      ax1=ax0+ex(m2)
      ay1=ay0+ey(m2)
      az1=az0+ez(m2)
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202

c     backup old path ------------>
 203  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      do i=m2,Lch
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      do i=m2,Lch
         ex_n(i)=ax0+ex(i)
         ey_n(i)=ay0+ey(i)
         ez_n(i)=az0+ez(i)
         egx_n(i)=ax0+egx(i)
         egy_n(i)=ay0+egy(i)
         egz_n(i)=az0+egz(i)
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      do i=m2,Lch
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(m,Lch))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(m,Lch,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,Lch,1)+ESHORT(m,Lch,1) !use ifs

         Enew=ehb(m,Lch,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,Lch,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=m,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(m,Lch,-1)+ESHORT(m,Lch,-1)

         Eold=ehb(m,Lch,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,Lch,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(18)=bNa(18)+1
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
            do kkk=m,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            do i=m2,Lch
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(18)=bNt(18)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ tran_C finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0)
c     This movement include an off-lattice rotation of fragment
c     and 2 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3, n3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tran_M(nfr_i,nfr_f,ifr,d_xyz00,ex,ey,ez,x,y,z,anvec,
     +                  goodc,ica,vx,vy,vz,Nw,W21,W22,egx,egy,egz,ecx,
     +                  ecy,ecz,ebx,eby,ebz,Lch,nop,nopp,nom,nomm,noa,
     +                  noaa,afsn,afs,eprofn,seq,envir,id,bNa,bNNa,
     +                  icnto,icnt,sumcto,sumct,eprofo,codevsum,codev,
     +                  didevsum,didev,bNt,bNNt,itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)
      parameter(apai=3.1415926)

      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension NOPP(ndim),NOMM(ndim),NOAA(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension d_xyz00(ndim),angle00(ndim)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_i(ndim),nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      real A(3,3),E(3),T(3,3)





CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)

CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)

      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** translate the fragment **************************
***   
c     d_xyz0=0.5                  !maximum translation coordinates
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
***   check m-terminal of fragment--------------->
      ax1=ax0+ex(m2)
      ay1=ay0+ey(m2)
      az1=az0+ez(m2)
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202
***   check n-terminal of fragment--------------->
 203  ax1=ax0+ex(n2)
      ay1=ay0+ey(n2)
      az1=az0+ez(n2)
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=m2,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=m2,n2
         ex_n(i)=ax0+ex(i)
         ey_n(i)=ay0+ey(i)
         ez_n(i)=az0+ez(i)
         egx_n(i)=ax0+egx(i)
         egy_n(i)=ay0+egy(i)
         egz_n(i)=az0+egz(i)
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=m2,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(m,n))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(m,n,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs

         Enew=ehb(m,n,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,n,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

         Eold=ehb(m,n,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,n,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(17)=bNa(17)+1
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
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=m2,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(17)=bNt(17)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ tran_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0) + deformation (angle0)
c     Move the off-lattice N-fragment.
c     This movement include an off-lattice rotation of fragment
c     and a two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defo_N(nfr_f,ifr,d_xyz00,angle00,ex,ey,ez,x,y,z,anvec,
     +                  goodc,ica,vx,vy,vz,Nw,W21,W22,egx,egy,egz,ecx,
     +                  ecy,ecz,ebx,eby,ebz,Lch,nop,nopp,nom,nomm,noa,
     +                  noaa,afsn,afs,seq,eprofn,envir,id,bNa,bNNa,
     +                  armsd_sum,N_rmsd,icnto,icnt,sumcto,sumct,
     +                  eprofo,codevsum,codev,didevsum,didev,bNt,bNNt,
     +                  itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)

      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      dimension d_xyz00(ndim),angle00(ndim)
      dimension armsd_sum(100),N_rmsd(100)
      real A(3,3),E(3),T(3,3)




CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)
CM
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
CM      common/defoangle/defo_angle
CM      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)

      n=nfr_f(ifr)+2            !ending point of movement
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
ccc   d_xyz0=0.5               !range of translation movement
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
ccc   rotation axis (awx,awy,awz)->N-terminal; (bwx,bwy,bwz)->C-terminal:
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      bwx=asin_theta*cos(aphi)
      bwy=asin_theta*sin(aphi)
      bwz=acos_theta
ccc   rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy()*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
      bng=(aranzy()*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
ccc   rotation points:
      n_rot=1+int(n2*aranzy()) ![1,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
***   rotation matrix:
      acos=cos(ang)             !for N-terminal of the fragments--->
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
      acos=cos(ang)             !for C-terminal of the fragments--->
      asin=sin(ang)
      b11=bwx*bwx*(1-acos)+acos
      b12=bwx*bwy*(1-acos)-bwz*asin
      b13=bwx*bwz*(1-acos)+bwy*asin
      b21=bwx*bwy*(1-acos)+bwz*asin
      b22=bwy*bwy*(1-acos)+acos
      b23=bwy*bwz*(1-acos)-bwx*asin
      b31=bwx*bwz*(1-acos)-bwy*asin
      b32=bwy*bwz*(1-acos)+bwx*asin
      b33=bwz*bwz*(1-acos)+acos
***   check n-terminal of fragment--------------->
      ax1=ax0+ax+(ex(n2)-ax)*b11+(ey(n2)-ay)*b12+(ez(n2)-az)*b13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(n2)-ax)*b21+(ey(n2)-ay)*b22+(ez(n2)-az)*b23
      az1=az0+az+(ex(n2)-ax)*b31+(ey(n2)-ay)*b32+(ez(n2)-az)*b33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=1,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=1,n_rot
         ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
         ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
         egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
         egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
         ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
         ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
         ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
         ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
         eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
         ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
      enddo
      do i=n_rot+1,n2
         ex_n(i)=ax0+ax+(ex(i)-ax)*b11+(ey(i)-ay)*b12+(ez(i)-az)*b13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*b21+(ey(i)-ay)*b22+(ez(i)-az)*b23
         ez_n(i)=az0+az+(ex(i)-ax)*b31+(ey(i)-ay)*b32+(ez(i)-az)*b33
         egx_n(i)=ax0+ax+(egx(i)-ax)*b11+(egy(i)-ay)*b12+(egz(i)-az)*b13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*b21+(egy(i)-ay)*b22+(egz(i)-az)*b23
         egz_n(i)=az0+az+(egx(i)-ax)*b31+(egy(i)-ay)*b32+(egz(i)-az)*b33
         ecx_n(i)=ecx(i)*b11+ecy(i)*b12+ecz(i)*b13 !cc
         ecy_n(i)=ecx(i)*b21+ecy(i)*b22+ecz(i)*b23
         ecz_n(i)=ecx(i)*b31+ecy(i)*b32+ecz(i)*b33
         ebx_n(i)=ebx(i)*b11+eby(i)*b12+ebz(i)*b13 !bb
         eby_n(i)=ebx(i)*b21+eby(i)*b22+ebz(i)*b23
         ebz_n(i)=ebx(i)*b31+eby(i)*b32+ebz(i)*b33
      enddo
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=1,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(1,n))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(1,n,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(1,n,1)+ESHORT(1,n,1) !use ifs

         Enew=ehb(1,n,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,n,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)


         do kkk=1,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(1,n,-1)+ESHORT(1,n,-1)

           Eold=ehb(1,n,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,n,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(25)=bNa(25)+1
            bNNa(itemp)=bNNa(itemp)+1
            armsd_sum(itemp)=armsd_sum(itemp)+armsd
            N_rmsd(itemp)=N_rmsd(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=1,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=1,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(25)=bNt(25)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ defo_N finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0) + deformation (angle0)
c     Move C-terminal of off-lattice fragment
c     This movement include an off-lattice rotation of fragment
c     and 1 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defo_C(nfr_i,ifr,d_xyz00,angle00,ex,ey,ez,x,y,z,anvec,
     +                  goodc,ica,vx,vy,vz,Nw,W21,W22,Lch,egx,egy,egz,
     +                  ecx,ecy,ecz,ebx,eby,ebz,nop,nopp,nom,nomm,noa,
     +                  noaa,afsn,afs,seq,eprofn,envir,id,bNa,bNNa,
     +                  armsd_sum,N_rmsd,icnto,icnt,sumcto,sumct,eprofo,
     +                  codevsum,codev,didevsum,didev,bNt,bNNt,itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)


      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension NOPP(ndim),NOMM(ndim),NOAA(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension d_xyz00(ndim),angle00(ndim)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_i(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      dimension armsd_sum(100),N_rmsd(100)
      real A(3,3),E(3),T(3,3)



CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
CM
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)
CM
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)
CM      common/defoangle/defo_angle
CM      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)

      m=nfr_i(ifr)-2            !starting point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
************** rotate the fragment **************************
***   
ccc   d_xyz0=0.5               !range of translation movement
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
ccc   rotation axis (awx,awy,awz)->N-terminal; (bwx,bwy,bwz)->C-terminal:
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      bwx=asin_theta*cos(aphi)
      bwy=asin_theta*sin(aphi)
      bwz=acos_theta
ccc   rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy()*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
      bng=(aranzy()*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
ccc   rotation points:
      n_rot=m2+int((Lch-m2+1)*aranzy()) ![m2,Lch]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
      acos=cos(ang)
      asin=sin(ang)
      b11=bwx*bwx*(1-acos)+acos
      b12=bwx*bwy*(1-acos)-bwz*asin
      b13=bwx*bwz*(1-acos)+bwy*asin
      b21=bwx*bwy*(1-acos)+bwz*asin
      b22=bwy*bwy*(1-acos)+acos
      b23=bwy*bwz*(1-acos)-bwx*asin
      b31=bwx*bwz*(1-acos)-bwy*asin
      b32=bwy*bwz*(1-acos)+bwx*asin
      b33=bwz*bwz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202

c     backup old path ------------>
 203  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      do i=m2,Lch
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      do i=m2,n_rot
         ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
         ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
         egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
         egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
         ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
         ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
         ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
         ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
         eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
         ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
      enddo
      do i=n_rot+1,Lch
         ex_n(i)=ax0+ax+(ex(i)-ax)*b11+(ey(i)-ay)*b12+(ez(i)-az)*b13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*b21+(ey(i)-ay)*b22+(ez(i)-az)*b23
         ez_n(i)=az0+az+(ex(i)-ax)*b31+(ey(i)-ay)*b32+(ez(i)-az)*b33
         egx_n(i)=ax0+ax+(egx(i)-ax)*b11+(egy(i)-ay)*b12+(egz(i)-az)*b13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*b21+(egy(i)-ay)*b22+(egz(i)-az)*b23
         egz_n(i)=az0+az+(egx(i)-ax)*b31+(egy(i)-ay)*b32+(egz(i)-az)*b33
         ecx_n(i)=ecx(i)*b11+ecy(i)*b12+ecz(i)*b13 !cc
         ecy_n(i)=ecx(i)*b21+ecy(i)*b22+ecz(i)*b23
         ecz_n(i)=ecx(i)*b31+ecy(i)*b32+ecz(i)*b33
         ebx_n(i)=ebx(i)*b11+eby(i)*b12+ebz(i)*b13 !bb
         eby_n(i)=ebx(i)*b21+eby(i)*b22+ebz(i)*b23
         ebz_n(i)=ebx(i)*b31+eby(i)*b32+ebz(i)*b33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      do i=m2,Lch
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM      if(look(m,Lch))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(m,Lch,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,Lch,1)+ESHORT(m,Lch,1) !use ifs

         Enew=ehb(m,Lch,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,Lch,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=m,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(m,Lch,-1)+ESHORT(m,Lch,-1)

         Eold=ehb(m,Lch,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,Lch,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(27)=bNa(27)+1
            bNNa(itemp)=bNNa(itemp)+1
            armsd_sum(itemp)=armsd_sum(itemp)+armsd
            N_rmsd(itemp)=N_rmsd(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            do i=m2,Lch
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(27)=bNt(27)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ rot_C finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Translation (d_xyz0) + deformation (angle).
c     This movement include an off-lattice rotation of fragment
c     and 2 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3, n3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defo_M(nfr_i,ifr,nfr_f,d_xyz00,angle00,ex,ey,ez,x,y,z,
     +                  anvec,
     +                  goodc,ica,vx,vy,vz,Nw,W21,W22,Lch,nop,nopp,nom,
     +                  nomm,noa,noaa,afsn,afs,seq,eprofn,envir,id,bNa,
     +                  bNNa,armsd_sum,N_rmsd,icnto,icnt,sumcto,sumct,
     +                  eprofo,codevsum,codev,didevsum,didev,bNt,bNNt,
     +                  itemp)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=312)
      parameter(apai=3.1415926)

      logical look, goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension afs(ndim),afsn(ndim)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim)
      dimension bNa(100),bNt(100),bNNa(100),bNNt(100)
      dimension d_xyz00(ndim),angle00(ndim)
      dimension w21(-10:10,-10:10,-10:10,46)
      dimension w22(-10:10,-10:10,-10:10,46)
      dimension Nw(-10:10,-10:10,-10:10)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension nfr_i(ndim),nfr_f(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !CA
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      dimension armsd_sum(100),N_rmsd(100)
      real A(3,3),E(3),T(3,3)




CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/temperature/itemp,atemp
CM      common/lengths/Lch,Lch1,Lch2
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
CM      common/naccept1/bN5a(100),bN5t(100)
CM      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
CM      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb

CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)
CM
CM      common/mov2/v21(nvec,nvec,46),v22(nvec,nvec,46),Np2(nvec,nvec)
CM      common/w2a/w21(-10:10,-10:10,-10:10,46)
CM      common/w2b/w22(-10:10,-10:10,-10:10,46)
CM      common/w2c/Nw(-10:10,-10:10,-10:10)
CM
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/chainm/mv(ndim)
CM      common/defoangle/defo_angle
CM      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)

      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
ccc   d_xyz0=0.5                  !maximum translation coordinates
      ax0=(aranzy()*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy()*2-1)*d_xyz00(ifr)
      az0=(aranzy()*2-1)*d_xyz00(ifr)
ccc   rotation axis (awx,awy,awz)->N-terminal; (bwx,bwy,bwz)->C-terminal:
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
      acos_theta=1-2.*aranzy()
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy() ![0,2*pai]
      bwx=asin_theta*cos(aphi)
      bwy=asin_theta*sin(aphi)
      bwz=acos_theta
ccc   rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy()*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
      bng=(aranzy()*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
ccc   rotation points:
      n_rot=m2+int((n2-m2+1)*aranzy()) ![m2,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
ccc   rotation matrix:
      acos=cos(ang)             !for N-terminal of the fragments--->
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
      acos=cos(bng)             !for N-terminal of the fragments--->
      asin=sin(bng)
      b11=bwx*bwx*(1-acos)+acos
      b12=bwx*bwy*(1-acos)-bwz*asin
      b13=bwx*bwz*(1-acos)+bwy*asin
      b21=bwx*bwy*(1-acos)+bwz*asin
      b22=bwy*bwy*(1-acos)+acos
      b23=bwy*bwz*(1-acos)-bwx*asin
      b31=bwx*bwz*(1-acos)-bwy*asin
      b32=bwy*bwz*(1-acos)+bwx*asin
      b33=bwz*bwz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202
***   check n-terminal of fragment--------------->
 203  ax1=ax0+ax+(ex(n2)-ax)*b11+(ey(n2)-ay)*b12+(ez(n2)-az)*b13
      ay1=ay0+ay+(ex(n2)-ax)*b21+(ey(n2)-ay)*b22+(ez(n2)-az)*b23
      az1=az0+az+(ex(n2)-ax)*b31+(ey(n2)-ay)*b32+(ez(n2)-az)*b33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.25)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy()*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy()*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy()*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=m2,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=m2,n_rot
         ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
         ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
         egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
         egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
         ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
         ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
         ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
         ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
         eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
         ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
      enddo
      do i=n_rot+1,n2
         ex_n(i)=ax0+ax+(ex(i)-ax)*b11+(ey(i)-ay)*b12+(ez(i)-az)*b13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*b21+(ey(i)-ay)*b22+(ez(i)-az)*b23
         ez_n(i)=az0+az+(ex(i)-ax)*b31+(ey(i)-ay)*b32+(ez(i)-az)*b33
         egx_n(i)=ax0+ax+(egx(i)-ax)*b11+(egy(i)-ay)*b12+(egz(i)-az)*b13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*b21+(egy(i)-ay)*b22+(egz(i)-az)*b23
         egz_n(i)=az0+az+(egx(i)-ax)*b31+(egy(i)-ay)*b32+(egz(i)-az)*b33
         ecx_n(i)=ecx(i)*b11+ecy(i)*b12+ecz(i)*b13 !cc
         ecy_n(i)=ecx(i)*b21+ecy(i)*b22+ecz(i)*b23
         ecz_n(i)=ecx(i)*b31+ecy(i)*b32+ecz(i)*b33
         ebx_n(i)=ebx(i)*b11+eby(i)*b12+ebz(i)*b13 !bb
         eby_n(i)=ebx(i)*b21+eby(i)*b22+ebz(i)*b23
         ebz_n(i)=ebx(i)*b31+eby(i)*b32+ebz(i)*b33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=m2,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
      enddo

CM       if(look(m,n))then         ! check excluded volumn for passage of [m,n]
      if(LOOK(m,n,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)) then
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
CM         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs

         Enew=ehb(m,n,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,n,1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
CM         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

         Eold=ehb(m,n,-1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(m,n,-1,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

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

CM         call metro(dE,id) !id=1, accepted; id=3, rejected
         call metro(dE,id,atemp,ichos,aTs) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(26)=bNa(26)+1
            bNNa(itemp)=bNNa(itemp)+1
            armsd_sum(itemp)=armsd_sum(itemp)+armsd
            N_rmsd(itemp)=N_rmsd(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=m2,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
            enddo
CM            call get_center     !update (amx,amy,amz)
        call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(26)=bNt(26)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ defo_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccc Try to swap i1-replica and i2-replica ccccccccccccccccc
c     swap the replicas with combined chains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine swap_RS(i1,i2)
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1500)	!maximum length of chain-length
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

CM Added metro_swap function here instaed of calling

CM      w_12=wi(j)wj(i)/wi(i)wj(j) ------------->
      aaa=(1/aT_rep(i)-1/aT_rep(j))*(E_rep(i)-E_rep(j))
      aweight=exp(aaa)
  
      if(aranzy().le.aweight)then
        metro_swap=1           !swap accepted
      else
        metro_swap=3           !swap rejected
      endif

CMM      if(metro_swap(i1,i2).eq.1)then
      if(metro_swap.eq.1)then
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


ccccc record the centroid
c     record EE(3),HH(3)
      subroutine get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)

      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension seq(ndim)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension mv(ndim)
      real AA(3,3),EE(3),HH(3,3)

CM      real A(3,3),E(3),T(3,3)


CM      common/lengths/Lch,Lch1,Lch2
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      common/chain0/ras(ndim),nfl
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim)
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)
CM      common/chainm/mv(ndim)
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)      
CM      common/center/cex,cey,cez
CM      common/eigen/AA(3,3),EE(3),HH(3,3)

      print *,'In get_center',Lch,switch

cccc  find center of mass of CA's
      cex=0
      cey=0
      cez=0
c     all residues are moveable
      if(switch.eq.1)then
         do i=1,Lch
            cex=cex+x(i)
            cey=cey+y(i)
            cez=cez+z(i)
         enddo
c     there are some frozen fragments
      else                     
         do i=1,Lch
            if(mv(i).gt.0)then
               cex=cex+x(i)
               cey=cey+y(i)
               cez=cez+z(i)
            else
               cex=cex+ex(i)
               cey=cey+ey(i)
               cez=cez+ez(i)
            endif
         enddo
      endif

      cex=cex/float(Lch)
      cey=cey/float(Lch)
      cez=cez/float(Lch)      

cccc  find tensor of inertia components of SG's
      AA(1,1)=0
      AA(1,2)=0
      AA(1,3)=0
      AA(2,2)=0
      AA(2,3)=0
      AA(3,3)=0

      do i=1,Lch
c     (agxi,agyi,agzi) are coordinates of SG, assuming center of mass of
c     CA's in origin
         if(mv(i).gt.0)then
            agxi=x(i)+gx(ica(i-1),ica(i),seq(i))-cex
            agyi=y(i)+gy(ica(i-1),ica(i),seq(i))-cey
            agzi=z(i)+gz(ica(i-1),ica(i),seq(i))-cez
         else
            agxi=egx(i)-cex
            agyi=egy(i)-cey
            agzi=egz(i)-cez
         endif
         AA(1,1)=AA(1,1)+agxi*agxi
         AA(1,2)=AA(1,2)+agxi*agyi
         AA(1,3)=AA(1,3)+agxi*agzi
         AA(2,2)=AA(2,2)+agyi*agyi
         AA(2,3)=AA(2,3)+agyi*agzi
         AA(3,3)=AA(3,3)+agzi*agzi
      enddo
      AA(1,1)=AA(1,1)/float(Lch)
      AA(1,2)=AA(1,2)/float(Lch)
      AA(1,3)=AA(1,3)/float(Lch)
      AA(2,2)=AA(2,2)/float(Lch)
      AA(2,3)=AA(2,3)/float(Lch)
      AA(3,3)=AA(3,3)/float(Lch)
      AA(2,1)=AA(1,2)
      AA(3,1)=AA(1,3)
      AA(3,2)=AA(2,3)

CM      print *,'AAs',AA(1,1),AA(1,2),AA(1,3),AA(2,2),AA(2,3),
CM     +              AA(3,3),AA(2,1),AA(3,1),AA(3,2)

c     obtain ellipcities and principal axis for the inertia tensor
      call eigenvalue(AA,EE,HH)
      return
      end
ccccc get_center finished


cccccccccccccccc Metropolis: id=1, accepted; id=3, rejected ccccccccccccccccc
      subroutine metro(dE,id,atemp,ichos,aTs)

      print *,'In metro dE=',dE,id,atemp,ichos,aTs
      id=1
      if(dE.gt.0)then
         if(aranzy().gt.weight12(dE,atemp,ichos,aTs))then
            id=3
         endif
      endif
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



ccccc read initial (x,y,z) from chain.dat
      subroutine read_initial2(Lch,q,cx0,cy0,cz0,cx,cy,cz,x,y,z,vx,vy,
     +                         vz,xrep,yrep,zrep,N_rep,nfl,ras,nfr,
     +                         L_cut,exc,amcheck_dis,goodc,mv,Mend_N,
     +                         Mend,Mend_C,nfl6,ras6,nfl5,ras5,nfl4,
     +                         ras4,nfl3,ras3,nfl2,ras2)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)


      logical look, goodc(nvec,nvec)
      dimension x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension cx(1000),cy(1000),cz(1000)
      dimension xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)
      dimension ras(ndim),ras2(ndim),ras3(ndim),ras4(ndim),ras5(ndim)
      dimension ras6(ndim)
      dimension mv(ndim)

      dimension ax(ndim),ay(ndim),az(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)
      dimension M_i(1000),M_f(1000)
      dimension n_i(1000),n_f(1000)

CM      common/logica/goodc
CM      logical look, goodc(nvec,nvec)
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/commonuse2/atemp1,atemp2,N_rep,phot
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      character c1*4,c2*2,c3*3,aaaaa*10,text,text1*22
CM      character*80 head
CM      character sign(ndim)
CM      common/initialstr/cx(1000),cy(1000),cz(1000)
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/looks/exc,exc1,exc2
CM      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
CM      dimension ax(ndim),ay(ndim),az(ndim)
CM      dimension x0(ndim),y0(ndim),z0(ndim)
CM      dimension M_i(1000),M_f(1000)
CM      common/chain0/ras(ndim),nfl
CM      common/gapmove2/ras2(ndim),nfl2
CM      common/gapmove3/ras3(ndim),nfl3
CM      common/gapmove4/ras4(ndim),nfl4
CM      common/gapmove5/ras5(ndim),nfl5
CM      common/gapmove6/ras6(ndim),nfl6
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim)
CM      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
CM      common/chainm/mv(ndim)
CM      dimension q_bk(ndim)
CM      character*3 sequ
CM      common/aminoacid/sequ(ndim)
CM      common/bigbond/i_bigbond,teco
CM      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)
CM      common/consensus1/cx0(0:1000),cy0(0:1000),cz0(0:1000),q(0:1000)
CM      dimension n_i(1000),n_f(1000)

      print *,'In read_initial2',Lch,L_cut

c     replica index
      i_rep=0                   
c     unit=24 <--> chain.dat
 104  rewind(24)
c     Number of templates in 'chain.dat'
      read(24,*)n_thr           

ccccc go through every template
      do 102 i_init=1,n_thr
c     q(i)=1 if the template contains coordinates for amino acid at
c     position "i", that is, amino acid "i" is aligned
         do i=1,Lch
            q(i)=0
         enddo
c     L_ali number of aligned residues (L_ali<=Lch)
         read(24,*)L_ali
         do i=1,L_ali
c     ii:amino acid number, (a1,a2,a3) coordinates
            read(24,*)ii,a1,a2,a3
c     store coordinates in temporary array
            cx0(ii)=a1
            cy0(ii)=a2
            cz0(ii)=a3
            q(ii)=1
         enddo
CM  1237    format(A22,I4,A4,3F8.3)

ccccc remove small segments, so that we have a better random walk
         M_a=0                  
c     0 does not correspond to any amino acid position. It's just a
c     boundary condition
         q(0)=0
ccccc find all aligned residues signaling the begginning of an contiguos
c     fragment of aligned residues
         do i=1,Lch
c     if the current position is aligned and the previous position is
c     unaligned, then...
            if(q(i).ne.q(i-1).and.q(i).eq.1)then
c     M_a: number of aligned segments
               M_a=M_a+1
c     M_i(l) initial amino acid position of a contiguous segment of
c     aligned residue. "l" is a running index from 1 to M_a to list all
c     aligned fragments
               M_i(M_a)=i
            endif
c     M_f(l) final amino acid position of a contiguous segment of
c     aligned residue. "l" is a running index from 1 to M_a to list all
c     aligned fragments
            if(q(i).eq.1)M_f(M_a)=i
CM            print *,'M_i M_F',i,M_i(M_a),M_f(M_a)
         enddo

ccccc go through all initial positions of the aligned fragments and set
c     them as unaligned (q(i)=0) if they are too small (L_cut=5)
         
         do i=1,M_a
            L_a=M_f(i)-M_i(i)+1
            if(L_a.le.L_cut)then
               do j=M_i(i),M_f(i)
                  q(j)=0
               enddo
            endif
         enddo

ccccc find number of unaligned fragments (N_g) and their initial (
c     n_i(l) ) and final ( n_f(l) ) positions
         N_g=0                  
         q(0)=1
         do i=1,Lch
            if(q(i).ne.q(i-1).and.q(i).eq.0)then
               N_g=N_g+1
               n_i(N_g)=i       
            endif
            if(q(i).eq.0)n_f(N_g)=i
         enddo

CM      print*,'n_g',n_g

c     number of times we try to build the all the unaligned fragments
         n_walk=0               !for number of reject by excluded volumn
c     exc=23 (corresponds to 4.2Angstroms, squared)
         exc_eff=exc            !for excluded volumn
 70      n_walk=n_walk+1
CM         print*,'NWALK',n_walk
         if(n_walk.gt.1000)then
c     reduce exclude volume parameter to increase success
            exc_eff=exc_eff*0.99
         endif
         if(n_walk.gt.10000)then
            write(20,*)'unsolvable structure',i
            write(*,*)'unsolvable structure',i
c     exit program with non-zero (error) code            
            call exit(1)
         endif

ccccc fill cx(i),cu(i),cz(i)
         do i=1,Lch
c     q(i)=1 means aligned residue, thus store its coordinates
            if(q(i).eq.1)then
               cx(i)=cx0(i)
               cy(i)=cy0(i)
               cz(i)=cz0(i)
            else
c     put some stupid initial coordinates for unaligned residue
               cx(i)=1000000.   
               cy(i)=1000000.
               cz(i)=1000000.
            endif
CM            print *,'cx cy cz',i,cx(i),cy(i),cz(i)
         enddo

ccccc find a random walk for each unaligned fragment
         do i=2,n_g
c     (i1,i2) initial and final amino acid positions of the unaligned
c     fragment
            i1=n_i(i)
            i2=n_f(i)
CM            print *,'In n_g loop',i,i1,i2
c     assign coordinates to unaligned residues
CM            call connect(i1,i2,pass)
           call connect(i1,i2,pass,cx,cy,cz,Lch,amcheck_dis) !fill missed cooridinates, no move others
CM      print*,'After Connect 1 pass',pass,amcheck_dis
c     unsuccessful buildinf of the unaligned fragment. Maybe other
c     unaligned fragments already constructed make it impossible. Thus,
c     trash all unaligned fragments already built and begin afresh
            if(pass.ge.3)goto 70
         enddo

         
         if(n_g.ge.1)then
            i1=n_i(1)
            i2=n_f(1)
CM            call connect(i1,i2,pass)
            call connect(i1,i2,pass,cx,cy,cz,Lch,amcheck_dis) !fill missed cooridinates, no move others
CM      print*,'After Connect 2 pass',pass,amcheck_dis
            if(pass.ge.3)goto 70 !re-walk
         endif


ccccc project chain onto lattices, decide (x,y,z):

c     rescale to lattice units and store in (ax(i),ay(i),az(i))
         do i=1,Lch
            ax(i)=cx(i)/0.87  
            ay(i)=cy(i)/0.87
            az(i)=cz(i)/0.87
         enddo

c     (ax,ay,az) into (x,y,z)
c     first residue is easy. Just take the integer part
         x(1)=nint(ax(1))
         y(1)=nint(ay(1))
         z(1)=nint(az(1))
         xm=x(1)
         ym=y(1)
         zm=z(1)
c     continue on next Lch-1 residues
         do 101 i=2,Lch
c     minimum distance between c_alpha and lattice
            dis_min=1000000.    
            j_ch=0              

cccc  find a good peptide vector to join next residue
            do 100 j=1,nvec
c     check peptide bond "j" makes angle in [65,165] with previous
c     peptide bond "jm"
               if(i.gt.2)then
                  if(.not.goodc(jm,j))goto 100
               endif
c     (xm,ym,zm) lattice-coordinates of previous residue
               x_tmp=xm+vx(j)
               y_tmp=ym+vy(j)
               z_tmp=zm+vz(j)
c     check excluded volume with previous [1,i-3] residues
               do m=1,i-3  
                  disaa=(x_tmp-x(m))**2+(y_tmp-y(m))**2+(z_tmp-z(m))**2
CM                  print *,'disaa',disaa
                  if(disaa.lt.exc_eff) goto 100
               enddo
c     check distance to off-lattice coordinates
               dis=(ax(i)-x_tmp)**2+(ay(i)-y_tmp)**2+(az(i)-z_tmp)**2
CM               print *,'dis',dis
c     if this is the vector found to be closest to the off-lattice
c     coordinates, then...
               if(dis.lt.dis_min)then
                  j_ch=j
                  x_ch=x_tmp
                  y_ch=y_tmp
                  z_ch=z_tmp
                  dis_min=dis
               endif
 100        continue

c     if we could not find a peptide bond, then trash all the built
c     unalined residues and start afressh
CM      print*,'j_ch',j_ch
            if(j_ch.lt.1)goto 70
c     store lattice coordinates
            x(i)=x_ch           
            y(i)=y_ch
            z(i)=z_ch
c     update peptide bond index for previous residue
            jm=j_ch
c     update coordinates of previous residue
            xm=x(i)
            ym=y(i)
            zm=z(i)
CM            print *,'xm ym zm',i,xm,ym,zm
c     go to next residue to project on-lattice
 101     continue


ccccc record the initial conformation of current replica
         i_rep=i_rep+1
         do i=1,Lch
            xrep(i,i_rep)=x(i)
            yrep(i,i_rep)=y(i)
            zrep(i,i_rep)=z(i)
CM            print *,'xrep etc',i,i_rep,xrep(i,i_rep),
CM     +               yrep(i,i_rep),zrep(i,i_rep)
         enddo
c     stop reading templates from chain.dat if we reached required
c     number of replicas
         if(i_rep.ge.N_rep)goto 105
c     continue to next template
 102  continue
c     read again chain.dat, and start with the first template, if the
c     number of templates stored in chain.dat is less than the required
c     number of replicas
      if(i_rep.lt.N_rep)goto 104

ccccc prepare movement for normal movement
c     all residues are flexible
 105  nfl=Lch
      do i=1,nfl
c     amino acid position of i-th moveable residue. Since all residues
c     are flexible, the ras(i) allways coincides with "i"
         ras(i)=i
      enddo
c     set up amino acids positions that can be selected as starting
c     points for the different MC moves
CM      call move_point           
      call move_point(Lch,mv,nfl,ras,Mend_N,Mend,Mend_C,
     +                      nfl6,ras6,nfl5,ras5,nfl4,ras4,nfl3,
     +                      ras3,nfl2,ras2)
c     number of frozen fragments.
      nfr=0                     
      return
      end
ccccc read_initial2 finished

ccccccccccccccc Try to swap i1-replica and i2-replica ccccccccccccccccc
      subroutine swap(i1,i2,bNSa,Lch,xrep,yrep,zrep,aT_rep,E_rep,bNSt)
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1500)	!maximum length of chain-length
      parameter(nvec=312)

      dimension bNSa(100),bNSt(100)
      dimension aT_rep(nrep),E_rep(nrep)
      dimension xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)



CM      common/lengths/Lch,Lch1,Lch2
CM      common/nswap/bNSa(100),bNSt(100)
CM      common/sw1/aT_rep(nrep),E_rep(nrep)
CM      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
CM      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
CM      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
CM      common/chain0/ras(ndim),nfl
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      
CM Added metro_swap function here instaed of calling

CM      w_12=wi(j)wj(i)/wi(i)wj(j) ------------->
      aaa=(1/aT_rep(i1)-1/aT_rep(i2))*(E_rep(i1)-E_rep(i2))
      aweight=exp(aaa)
  
      if(aranzy().le.aweight)then
        metro_swap=1           !swap accepted
      else
        metro_swap=3           !swap rejected
      endif

CMM      if(metro_swap(i1,i2).eq.1)then
      if(metro_swap.eq.1)then
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



ccccc count the number of restraints satisfied
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine count_restrains(Lch,x,y,z,xrep,yrep,zrep,ica,vector,gx,
     +                           gy,gz,seq,Mcom,Kcom,Mdis,kdis,t_dist,
     +                           s_dist,dev,McomCA,KcomCA,t_combCA,
     +                           dist_CA_cut,s_combCA,
     +                           t_comb,concut2,s_comb,dist,Lch1)
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1500)	!maximum length of chain-length
      parameter(nvec=312)


      dimension xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension vector(-5:5,-5:5,-5:5)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension seq(ndim)
      dimension concut2(0:19,0:19)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension Mcom(ndim),Kcom(ndim,100)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)



CM      common/lengths/Lch,Lch1,Lch2
CM      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
CM      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
CM      common/CAcontact1/dist_CA_cut
CM      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA

      dimension ax(ndim),ay(ndim),az(ndim) !SG
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim) !CA


      do i=1,3
         do j=1,Lch
            x(j)=xrep(j,i)
            y(j)=yrep(j,i)
            z(j)=zrep(j,i)
         enddo
         do j=1,Lch1
            wx=x(j+1)-x(j)
            wy=y(j+1)-y(j)
            wz=z(j+1)-z(j)
            ica(j)=vector(wx,wy,wz) !identify order number of each bond-vector
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         do j=1,Lch
            ax(j)=x(j)+gx(ica(j-1),ica(j),seq(j))
            ay(j)=y(j)+gy(ica(j-1),ica(j),seq(j))
            az(j)=z(j)+gz(ica(j-1),ica(j),seq(j))
         enddo

         do j=1,Lch
c     contact restraints--------->
            do k=1,Mcom(j)
               m=Kcom(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_comb=t_comb+1
                  dis=(ax(j)-ax(m))**2+(ay(j)-ay(m))**2
     $                 +(az(j)-az(m))**2
c                  write(*,*)j,m,dis,concut2(seq(j),seq(m)),'comb'
                  if(dis.le.concut2(seq(j),seq(m)))then
c                     write(*,*)'***'
                     s_comb=s_comb+1
                  endif
               endif
            enddo
c     dist restraints--------->
            do k=1,Mdis(j)
               m=kdis(j,k)      !j-m restraints
               if(m.gt.j) then  !to avoid repeat
                  t_dist=t_dist+1
                  dis2=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  dis=sqrt(dis2)
                  err=abs(dis-dist(j,k)) !dist: predicted dis
c                  write(*,*)j,m,err,dev(j,k),'dist'
                  if(err.lt.dev(j,k)) then !dev: deviation for arca
c                     write(*,*)'---'
                     s_dist=s_dist+1
                  endif
               endif
            enddo
c     CAcontact restraints--------->
            do k=1,McomCA(j)
               m=KcomCA(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_combCA=t_combCA+1
                  dis=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
c                  write(*,*)j,m,dis,dist_CA_cut,'combCA'
                  if(dis.le.dist_CA_cut)then
c                     write(*,*)'==='
                     s_combCA=s_combCA+1
                  endif
               endif
            enddo
         enddo
      enddo

c^^^^^^^^^^^^^^^^^ restraints count finished ^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccc 1, translate (x,y,z) to calculate (amx,amy,amz)
c     2, calculate initial energy, 
c     3, prepare initial parameters:
       subroutine initial_move(icnto,sumcto,Lch,nop,noa,nom,eprofo,seq,
     +                        envir,icnt,sumct,nopp,noaa,nomm,
     +                        codevsum,conew,didevsum,dinew,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2,
     +                hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb,vvv,exc,acut_cc,acut_bb,acut_vv,bcut_cc,
     +                bcut_bb,bcut_vv,dord,en2,en3,contt,
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,switch,
     +             Cr2a,Cr2b,eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,AA,
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,colim,codev,
     +                didev,jbin)


      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)


      dimension envir(0:15,0:15,0:15,0:19,4)
      dimension nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim),nhbnn(ndim)
      dimension seq(ndim),sec(ndim)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension ehbij(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      dimension arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      dimension ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension jbin(0:500),acops(ndim,16)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension mv(ndim)
      dimension apar(ndim,ndim)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)

      dimension IBIN(-300:300),asr(ndim,-12:12)
      dimension afs(ndim)
      dimension frga(ndim),frgb(ndim)
      dimension csr(ndim,2)
      dimension Mcom(ndim),Kcom(ndim,100)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      dimension aweig(4000,4000)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension area(ndim)
      dimension bgx(ndim),bgy(ndim),bgz(ndim)
      dimension AA(3,3),EE(3),HH(3,3)


CM      common/lengths/Lch,Lch1,Lch2
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/backup1/nopp(ndim),nomm(ndim),noaa(ndim),nhbnn(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) 
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim)
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)
CM      common/chainm/mv(ndim)

      dimension axx(1000),ayy(1000),azz(1000)

        print *,'In initial_move',Lch,exc

ccccc update vvv(i,j), ie, which position pairs (i,j) are to be checked
c     for excluded volume violations
CM      call get_vvv
      call get_vvv(Lch,vvv,mv,exc)
ccccc compute center of mass (cex,cey,cez) of CA's.
c     compute ellipcities E(i) and principial axis of SG's T(3,3)
CM      call get_center
      call get_center(switch,x,y,z,Lch,ex,ey,ez,AA,gx,gy,gz,ica,
     +                      seq,egx,egy,egz,mv,
     +                      cex,cey,cez,EE,HH)

ccccc prepare all the initial parameters of this replica
c     need when istat>0
      icnto=0                   
c     need when istat>0
      sumcto=0                  
      do i=1,Lch
c     number of parallel contacts that amino acid at position "i" makes
c     with rest of amino acids
         nop(i)=0               
c     number of antiparallel contacts
         noa(i)=0
c     number of orthogonal contacts
         nom(i)=0
      enddo

CM      energ=ehb(1,Lch,1)+eshort(1,Lch,10) !initialize all below parameters

      energ=ehb(1,Lch,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,Lch,10,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)
!initialize all below parameters

      print *,'energ',energ

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
c     eprofo stores the summ of the environment potential
      eprofo=0.0
      do k=1,Lch
         is=seq(k)
c     number of antiparallel contact-pairs
         ia=noa(k)            
c     number of parallel contact-pairs
         ip=nop(k)              
c     number of orgonal contact-apirs
         im=nom(k)              
         eprofo=eprofo+envir(ia,im,ip,is,3)
CM         print *,'eprofo=',k,eprofo
      enddo
c     backup of contact-order
      icnto=icnt                
      sumcto=sumct              !backup of contact-order
      do i=1,Lch
         nopp(i)=nop(i)         !backup of number of contact-pair
         noaa(i)=noa(i)         !backup of number of contact-pair
         nomm(i)=nom(i)         !backup of number of contact-pair
CM         print *,'nopp noaa nonn=',i,nopp(i),noaa(i),nomm(i)
      enddo
      codevsum=conew            !backup  panelity for total deviation of comb
      didevsum=dinew            !backup  panelity for total deviation of dist

      print *,'End of initial_move=',icnto,icnt,sumcto,sumct,
     +         codevsum,conew,didevsum,dinew
      return
      end
ccccc subroutine initial_move finished



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test whether all the neighboring backbones are good neighbors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine test_neighbor(Lch,ica,goodc,vx,vy,vz)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nrep=100)
      parameter(nvec=312)

      dimension ica(0:ndim)
      logical goodc(nvec,nvec)
      dimension vx(nvec),vy(nvec),vz(nvec)


CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/logica/goodc
CM      logical goodc(nvec,nvec)
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test for overlaps of Ca:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine test_overlap(Lch,exc,mv,x,y,z,ex,ey,ez)
      implicit integer(i-z)

      parameter(ndim=1500)

      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension mv(ndim)
      dimension x(ndim),y(ndim),z(ndim)

CM      common/lengths/Lch,Lch1,Lch2
CM      common/looks/exc,exc1,exc2

      do i=1,Lch
         do j=i+3,Lch
            dis2=(aax(i,mv,x,ex)-aax(j,mv,x,ex))**2+
     +           (aay(i,mv,y,ey)-aay(j,mv,y,ey))**2+
     $           (aaz(i,mv,z,ez)-aaz(j,mv,z,ez))**2
            if(dis2.lt.exc)then
               write(*,*)i,j,dis2,'  Ca-Ca overlap'
            endif
         enddo
      enddo

      return
      end


















CMMM FUNCTIONS


****************************************************
*     no=-1 (less than 0)                          *
*     rrr=ranzy(no)                                *
****************************************************
      function aranzy()
c     random numbers uniformly distributed between 0 and 1.  (code by
c     Y. Zhang, ITP, Acdemia Sincica, 1999.6 can be replaced by any
c     other suitable generator for random numbers)

      common/seed1/no ! will overwrite whatever value we pass with global no
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









CMM ****************************************************
CMM *     no=-1 (less than 0)                          *
CMM *     rrr=ranzy(no)                                *
CMM ****************************************************
CMM       function aranzy()
CMM CM      function aranzy(no,ix1,ix2,ix3,r)
CMM c     random numbers uniformly distributed between 0 and 1.  (code by
CMM c     Y. Zhang, ITP, Acdemia Sincica, 1999.6 can be replaced by any
CMM c     other suitable generator for random numbers)
CMM       common/seed1/no ! will overwrite whatever value we pass with global no
CMM       common/asdfa/ix1,ix2,ix3,rm1,rm2,r(99)
CMM CM      double precision r(99)
CMM 
CMM       data ia1,ic1,m1/1279,351762,1664557/
CMM       data ia2,ic2,m2/2011,221592,1048583/
CMM       data ia3,ic3,m3/15551,6150,29101/

CMM       if(no.ge.0) go to 2
CMM       ix1=mod(-no,m1)
CMM       ix1=mod(ia1*ix1+ic1,m1)
CMM       ix2=mod(ix1,m2)
CMM       ix1=mod(ia1*ix1+ic1,m1)
CMM       ix3=mod(ix1,m3)
CMM       rm1=1./float(m1)
CMM       rm2=1./float(m2)
CMM       do 1 j=1,99
CMM       ix1=mod(ia1*ix1+ic1,m1)
CMM       ix2=mod(ia2*ix2+ic2,m2)
CMM  1    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
CMM  2    ix1=mod(ia1*ix1+ic1,m1)
CMM       ix2=mod(ia2*ix2+ic2,m2)
CMM       ix3=mod(ia3*ix3+ic3,m3)
CMM       j=1+(99*ix3)/m3
CMM       aranzy=r(j)
CMM       r(j)=(float(ix1)+float(ix2)*rm2)*rm1
CMM       no=ix1
CMM       return
CMM       end



cccc return distance between two vectors
      function di(x1,y1,z1,x2,y2,z2)
      di=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
CM      print *,'In function di',x1,y1,z1,x2,y2,z2,di
      return
      end
cccc  di finished


ccccc returns x-coordinate of CA at position "i", whether on-lattice or
c     off-lattice
      function aax(i,mv,x,ex)
      implicit integer(i-z)
      parameter(ndim=1500)

      dimension mv(ndim),x(ndim),ex(ndim)

CM      common/chainm/mv(ndim)
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)

CM      print *,'In aax i=',i,mv(i)
      if(mv(i).gt.0)then
         aax=x(i)
      else
         aax=ex(i)
      endif
      return
      end
ccccc aax finished



ccccc returns y-coordinate of CA at position "i", whether on-lattice or
c     off-lattice
      function aay(i,mv,y,ey)
      implicit integer(i-z)
      parameter(ndim=1500)

      dimension mv(ndim),y(ndim),ey(ndim)

CM      common/chainm/mv(ndim)
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aay=y(i)
      else
         aay=ey(i)
      endif
      return
      end
ccccc aay finished



ccccc returns z-coordinate of CA at position "i", whether on-lattice or
c     off-lattice
      function aaz(i,mv,z,ez)
      implicit integer(i-z)
      parameter(ndim=1500)

      dimension mv(ndim),z(ndim),ez(ndim)

CM      common/chainm/mv(ndim)
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aaz=z(i)
      else
         aaz=ez(i)
      endif
      return
      end
ccccc aaz finished


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate energy from the beginning
c     parameters will be modified after running this subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energy_tot(icnto,sumcto,Lch,nop,noa,nom,seq,envir,
     +                    acorder,contt,en1,en2,en3,icnt,
     +                    sumct,codevsum,conew,didevsum,dinew,nopp,
     +                    noaa,nomm,
     +                    eprofo,
     +             ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             apar,alp,arlp,app,ala,arla,apa,
     +             alm,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin, 
     +                hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,ibin,
     +                frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,colim,codev,
     +                didev)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)
      parameter(nrep=100)

      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension vx(nvec),vy(nvec),vz(nvec)
      dimension nop(ndim),nom(ndim),noa(ndim)
      dimension nopp(ndim),nomm(ndim),noaa(ndim)
      dimension afs(ndim),afsn(ndim)
      real envir(0:15,0:15,0:15,0:19,4)
      dimension seq(ndim),sec(ndim)
      dimension bNa(nrep),bNt(nrep),bNNa(nrep),bNNt(nrep)
      dimension JBIN(0:500)
      dimension mv(ndim)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension bsr(ndim,16), acops(ndim,16)
      dimension apar(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      real arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      real ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension ehbij(ndim,ndim)
      dimension area(ndim)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      real csr(ndim,2),asr(ndim,-12:12)
      integer ibin(-300:300)
      dimension aweig(4000,4000)
      dimension AA(3,3),EE(3),HH(3,3)
      dimension vvv(ndim,ndim)










CM      common/lengths/Lch,Lch1,Lch2
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
CM      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
CM      common/backup2/eprofo,eprofn,energ
CM      common/seqe/seq(ndim),sec(ndim)
CM      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
CM      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
CM      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
CM      common/pair1/eh2,eh1b,eh1c
CM      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM      common/otherenergy/E_cord,E_cnum

CM      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,EHB5,EHB6
CM      common/ehbenergy1/EHB5a,EHB5b

      ICNTO=0
      SUMCTO=0
      DO i=1,Lch
         NOP(i)=0
         NOA(i)=0
         NOM(i)=0
      enddo

CM      energy_tot=EHB(1,Lch,1)+ESHORT(1,Lch,10)

      energ_tot=ehb(1,Lch,1,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2, 
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b, 
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin) +
     +             ESHORT(1,Lch,10,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb, 
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch, 
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)

      eprofo=0.0
      DO k=1,Lch
         is=seq(k)
         ia=noa(k)
         ip=nop(k)
         im=nom(k)
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c check whether passage of [jj,kk] overlap with other parts of chain.
c look =.ture., without overlap
c look =.false., with overlap
c only C_alpha is checked.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION LOOK(jj,kk,switch,Lch,exc,x,y,z,mv,ex,ey,ez,vvv)
      IMPLICIT INTEGER(I-Z)
      LOGICAL LOOK
      parameter(ndim=1500)
      parameter(nvec=312)

      dimension x(ndim),y(ndim),z(ndim)
      dimension mv(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension vvv(ndim,ndim)


CM      common/seqe/seq(ndim),sec(ndim) 	 	
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/looks/exc,exc1,exc2
CM      common/chainm/mv(ndim)
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      common/excluded/vvv(ndim,ndim)

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
            if(mv(k).gt.0)then  !off-lattice
               axk=x(k)
               ayk=y(k)
               azk=z(k)
            else                !on-lattice
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


ccccc calculate ehb for the contigous fragment [jjjj,kkkk]
c     1, pairwise interaction;
c     2, H-bond;
c     3, Electric;
c     4, contact number;
c     5, contact order.
      FUNCTION EHB(jjjj,kkkk,istat,ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,
     +             ehb4,ehb5a,enb5b,icnt,icnto,sumct,sumcto,seq,mv,x,y,
     +             z,ica,gx,gy,gz,HBX,HBY,HBZ,CAX,CAY,CAZ,ex,ey,ez,egx,
     +             egy,egz,ebx,eby,ebz,ecx,ecy,ecz,acops,
     +             Lch,apar,alp,nop,arlp,app,ala,noa,arla,apa,
     +             alm,nom,arlm,apm,sec,concut2,
     +             acut_cc,acut_bb,acut_vv,bcut_cc,bcut_bb,bcut_vv,
     +             dord,en2,en3,contt,Cr2a,Cr2b,
     +             eh1,eh2,eh3,eh4,eh1a,eh1b,eh1c,eh5a,eh5b,
     +             acut_hh,bcut_hh,EHBIJ,jbin)
      implicit integer(i-z)
      parameter(ndim=1500)
      parameter(nvec=312)

      dimension seq(ndim),sec(ndim)
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension ehbij(ndim,ndim)
      dimension apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      dimension arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      dimension ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      dimension jbin(0:500)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)	
      dimension mv(ndim)
      dimension apar(ndim,ndim)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension egx(ndim),egy(ndim),egz(ndim)
      dimension ecx(ndim),ecy(ndim),ecz(ndim)
      dimension ebx(ndim),eby(ndim),ebz(ndim)
      dimension concut2(0:19,0:19)
      dimension acops(ndim,16)







CM      common/seqe/seq(ndim),sec(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
CM      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
CM      common/energy/eh5,es3,es3a,es3b,eh1,eh1a,eh4,ehbij(ndim,ndim)
CM      common/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
CM      common/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
CM      common/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
CM      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
CM      common/short1/ jbin(0:500),bsr(ndim,16),acops(ndim,16)
CM      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
CM      common/icgg/ icg(ndim),eh6  		
CM      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3 
CM      common/envir/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)	
CM      common/pair1/eh2,eh1b,eh1c
CM      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
CM      common/ehbenergy/ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,ehb4,ehb5,ehb6
CM      common/distres/er4,es3c
CM      common/chainm/mv(ndim)
CM      common/hba/eh5a,cr2a,acut_bb,acut_cc,acut_vv,acut_hh
CM      common/hbb/eh5b,cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
CM      common/ehbenergy1/ehb5a,ehb5b
CM      common/par/apar(ndim,ndim)
CM      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim)
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)

c     +1/r of Ca-SG
      ehb1=0       
c     +1/r for non-parallel contact of Ca-Ca
      ehb1a=0                   
c     excluded volume of SG-SG
      ehb1b=0                 
c     pair-wise potential of SG-SG
      ehb1c=0                   
c     quarsi3 for SG-SG
      ehb2=0                    
c     enhance good-piece contacts
      ehb3=0                    
c     -1/r for parallel contact of Ca-Ca
      ehb4=0                    
c     H-bond energy for alpha-type
      ehb5a=0                    
c     H-bond energy for beta-type
      ehb5b=0                    

      print *,'In ehb',istat,jjjj,kkkk,icnto,sumcto,Lch
      print *,'In ehb acut bcut',acut_cc,acut_bb,acut_vv,
     +                           bcut_cc,bcut_bb,bcut_vv
      print *,'Initial EHB',eh1,eh1a,eh1b,eh1c,eh2,eh3,eh4,eh5a,eh5b,
     +         ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,ehb4,ehb5a,ehb5b

CM         do ijk = 0,500
CM            print *,'EHB begin jbin=',ijk,jbin(ijk)
CM         enddo

c     update contact order and
      if(istat.gt.0)then
         icnt=icnto
         sumct=sumcto
      endif

c     coupling of secondary structure with pairwise interactions
c     included - thus expanded range
      jj=jjjj-1
      if(jj.lt.1)jj=1
      kk=kkkk+1
      if(kk.gt.Lch)kk=Lch

CM      print *,'kk kk=',jj,kk
      
CM      do igx=1,nvec
CM      do jgx=1,nvec
CM      do kgx=0,19
CM         print *,'gx=',igx,jgx,kgx,gx(igx,jgx,kgx)
CM      enddo
CM      enddo
CM      enddo

CM      CALL EXIT(-1)

ccccc go through every residue of the extended fragment [jj,kk]
      do 1002 k=jj,kk
c     residue type at position "k"
         kseq=seq(k)
CM         print *,'kseq=',k,kseq,mv(k)
ccccc load CA and SG coordinates, also directions for hydrogen bond and
c     for orthogonal direction
         if(mv(k).gt.0)then   
c     (axk,ayk,azk) CA coordinates
            axk=x(k)            
            ayk=y(k)
            azk=z(k)
c     index of previous peptide bond. Remember we defined ica(0)=ica(2)
c     and ica(Lch)=ica(Lch-2)
            nv1=ica(k-1)
c     index of next peptide bond
            nv2=ica(k)
c     (agxk,agyk,agzk) SG coordinates
CM            print *,'HERE EHB', nv1,nv2,kseq
            agxk=axk+gx(nv1,nv2,kseq) 
            agyk=ayk+gy(nv1,nv2,kseq)
            agzk=azk+gz(nv1,nv2,kseq)
c     (bxk,byk,bzk) direction where H-bond can be made
            bxk=HBX(nv1,nv2)    
            byk=HBY(nv1,nv2)
            bzk=HBZ(nv1,nv2)
c     (cxk,cyk,czk) orthogonal direction to nv1 and nv2
            cxk=CAX(nv1,nv2)    
            cyk=CAY(nv1,nv2)
            czk=CAZ(nv1,nv2)
c     the residue is frozen
         else
            axk=ex(k)          
            ayk=ey(k)
            azk=ez(k)
            agxk=egx(k)   
            agyk=egy(k)
            agzk=egz(k)
            bxk=ebx(k)   
            byk=eby(k)
            bzk=ebz(k)
            cxk=ecx(k)   
            cyk=ecy(k)
            czk=ecz(k)
         endif

CM      print *,'if mv k gt 0=',axk,ayk,azk,agxk,agyk,agzk,
CM     +            bxk,byk,bzk,cxk,cyk,czk

         km2=k-2
         kp2=k+2
CM         print *,'km2 kp2=',km2,kp2

         if(km2.ge.1.and.kp2.le.Lch)then
CM            print *,'In km2 ge 1 kp2 le Lch'
c     calculate distance from km2 to kp2 (in lattice units)
            xxx=nint((aax(km2,mv,x,ex)-aax(kp2,mv,x,ex))**2+
     $           (aay(km2,mv,y,ex)-aay(kp2,mv,x,ey))**2+
     +           (aaz(km2,mv,z,ez)-aaz(kp2,mv,z,ez))**2)
CM            if(mv(km2).gt.0)then
CM              aax1=x(km2)
CM              aay1=y(km2)
CM              aaz1=z(km2)
CM            else
CM              aax1=ex(km2)
CM              aay1=ey(km2)
CM              aaz1=ez(km2)
CM            endif
CM            if(mv(kp2).gt.0)then
CM              aax2=x(kp2)
CM              aay2=y(kp2)
CM              aaz2=z(kp2)
CM            else
CM              aax2=ex(kp2)
CM              aay2=ey(kp2)
CM              aaz2=ez(kp2)
CM            endif
CM            xxx=nint((aax1-aax2)**2+(aay1-aay2)**2+(aaz1-aaz2)**2)
CM            print *,'xxx',km2,kp2,xxx
            if(xxx.gt.500)xxx=500
c     acops(i,k) is the average of bsr over the
c     [jbin(xxx)-1,jbin(xxx)+1] set of bins, or zero if the average is
c     negative. bsr(km2,l) is the r15 potential at residue "km2" when
c     [km2,km2+5] residues have a distance corresponding to bin "l"
            print *,'jbin=',xxx,jbin(xxx)
            ek5=acops(km2,jbin(xxx))
         else
            ek5=0
         endif

CM         print *,'ek5=',ek5

ccccc find interactions of amino acid "k" with rest of the protein
         do 1001 i=1,Lch
            iend=max(k+1,kk)
CM            print *,'iend in 1001 loop',i,iend

c     avoid counting twice the (i,k) interaction
CM            if(i.ge.k-1.and.i.le.iend)goto 1001 
            if(i.ge.k-1.and.i.le.iend) then
CM               print *,'if i ge k-1 i le iend go to 1001'
               goto 1001
            endif
c     load now coordinates of CA at "i"
            if(mv(i).gt.0)then
               axi=x(i)
               ayi=y(i)
               azi=z(i)
            else
               axi=ex(i)
               ayi=ey(i)
               azi=ez(i)
            endif
c     (axki,ayki,azki) relative vector going from CA of "i" to CA of "k"
            axki=(axk-axi)      
            ayki=(ayk-ayi)
            azki=(azk-azi)
c     square of relative vector (axki,ayki,azki)
            cr2=axki**2+ayki**2+azki**2+0.000001
CM            print *,'a etc cr2=',cr2,axki,ayki,azki
c     120 corresponds to 9.53Angstroms squared
            IF(cr2.lt.120)then
c     sequence separation
               idist=iabs(i-k)
c     residue type at position "i"
               iseq=seq(i)
CM               print *,'if cr2 lt 120=',idist,iseq,mv(i)
c     load peptide bonds flanking "i", coordinates of its SG, H-bond
c     direction, and orthogonal direction to the peptide bonds
c     (bisector)
               if(mv(i).gt.0)then
                  nv1=ica(i-1)
                  nv2=ica(i)
                  agxi=axi+gx(nv1,nv2,iseq)
                  agyi=ayi+gy(nv1,nv2,iseq)
                  agzi=azi+gz(nv1,nv2,iseq)
                  bxi=hbx(nv1,nv2)
                  byi=hby(nv1,nv2)
                  bzi=hbz(nv1,nv2)
                  cxi=cax(nv1,nv2)
                  cyi=cay(nv1,nv2)
                  czi=caz(nv1,nv2)	
               else
                  agxi=egx(i) 
                  agyi=egy(i)
                  agzi=egz(i)
                  bxi=ebx(i)  
                  byi=eby(i)
                  bzi=ebz(i)
                  cxi=ecx(i)  
                  cyi=ecy(i)
                  czi=ecz(i)
               endif
CM               print *,'if mv i gt 0 in 1001 loop =',agxi,agyi,agzi,
CM     +                 bxi,byi,bzi,cxi,cyi,czi

ccccc 1/r excluded for Ca-SG pair

c     if kseq is not GLY, then..
               if(kseq.gt.0) then
c     square distance from SG of "k" to CA of "i"
                  aks=(agxk-axi)**2+(agyk-ayi)**2+(agzk-azi)**2+.0001
CM                  print *,'in kseq gt 0 aks =',aks,ehb1
c     36 corresponds to 5.2Angstroms squared
                  if(aks.lt.36) then 
c     13 corresponds to 3.17Angstroms squared
                     if(aks.lt.13.0) aks=13.0
c     ehb1 is negative in the range [26,36], corresponding to
c     [4.4,5.2]Angstroms. ehb1 is zero above 36 ( or 5.2A)
                     ehb1=ehb1+((13.0/aks)-0.5)
                  endif
CM                  print *,'kseq',aks,ehb1
               endif

               if(iseq.gt.0) then
c     square distance from SG of "i" to CA of "k"
                  aks=(axk-agxi)**2+(ayk-agyi)**2+(azk-agzi)**2+.0001 !Ca_kSG_i
CM                  print *,'in iseq gt 0 aks =',aks,ehb1
                  if(aks.lt.36.0) then
                     if(aks.lt.13.0) aks=13.0
                     ehb1=ehb1+((13.0/aks)-0.5)
                  endif
CM                  print *,'iseq',aks,ehb1
               endif

ccccc pair-wise potential of SG-SG
               Gr2=(agxi-agxk)**2+(agyi-agyk)**2+(agzi-agzk)**2 !SG_i,SG_k
               if(Gr2.lt.concut2(iseq,kseq))then
                  EHB1c=EHB1c+apar(i,k)
               endif
CM               print *,'Gr2',i,Gr2,EHB1c

c     quarsi3 for SG-SG pair, 1/r for Ca-Ca ----------------->
               cc=cxi*cxk+cyi*cyk+czi*czk !c_i*c_k
CM               print *,'cc',i,cc
               IF(cc.gt.0.5)THEN !c_i//c_k
CM                  print *,'cc gt 0.5 Cr2=',Cr2
                  if(Cr2.lt.60)EHB4=EHB4-(30/(max(30.,Cr2))-0.5) !6.74A
CM                  print *,'if Cr2 lt 60 EHB4=',EHB4
CM                     print *,'Gr2 lt alp',Gr2,alp(iseq,kseq)
                  if(Gr2.lt.alp(iseq,kseq))then
CM                      print *,'HERE 1 Gr2 lt alp',i,idist,ek5,
CM     +                           EHB3
                     EHB3=EHB3-ek5*
     +                    ei5(i,idist,Lch,mv,x,y,z,ex,ey,ez,acops,jbin)
CM                     print *,'EHB3=',EHB3
                     NOP(k)=NOP(k)+ISTAT
                     NOP(i)=NOP(i)+ISTAT
                     ICNT=ICNT+istat*idist
                     SUMCT=SUMCT+ISTAT
CM                     print *,'arlp=',arlp(iseq,kseq),EHB2,EHB1b
                     if(Gr2.gt.arlp(iseq,kseq))THEN
                        EHB2=EHB2+app(i,k) !quarsi3
                     else
                        EHB1b=EHB1b+1 !excluded volume
                     endif
                  endif
CM                  print *,'cc 1 nop=',i,cc,EHB4,EHB3,NOP(k),NOP(i),icnt,
CM     +                  SUMCT,EHB2,EHB1b
               ELSE
                  if(Cr2.lt.33)EHB1a=EHB1a+(16.0/Cr2-0.5) !5A
CM                  print *,'if Cr2 lt 33 EHB1a=',EHB1a
                  IF(cc.lt.-0.5) THEN !antiparallel pair-wise of (i,k)
CM                     print *,'cc lt -0.5 Gr2=',Gr2,ala(iseq,kseq)
                     if(Gr2.lt.ala(iseq,kseq))then
CM                         print *,'HERE 2 Gr2 lt ala',i,idist,ek5,
CM     +                           EHB3
CM                        EHB3=EHB3-ek5*ei5(i,idist)
                        EHB3=EHB3-ek5*
     +                    ei5(i,idist,Lch,mv,x,y,z,ex,ey,ez,acops,jbin)
CM                        print *,'EHB3=',EHB3
                        NOA(k)=NOA(k)+ISTAT
                        NOA(i)=NOA(i)+ISTAT
                        ICNT=ICNT+istat*idist
                        SUMCT=SUMCT+ISTAT
CM                        print *,'arla=',arla(iseq,kseq),EHB2,EHB1b
                        if(Gr2.gt.arla(iseq,kseq))THEN
                           EHB2=EHB2+apa(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
CM                  print *,'cc 2 noa=',i,cc,EHB1a,EHB3,NOA(k),
CM     +                     NOA(i),icnt,SUMCT,EHB2,EHB1b
                  ELSE          !neither parallel nor antiparallel
CM                     print *,'Gr2 lt alm',Gr2,alm(iseq,kseq)
                     if(Gr2.lt.alm(iseq,kseq)) then
CM                        print *,'HERE Gr2 lt alm',i,idist,ek5,
CM     +                          EHB3
CM                        EHB3=EHB3-ek5*ei5(i,idist)
                         EHB3=EHB3-ek5*
     +                    ei5(i,idist,Lch,mv,x,y,z,ex,ey,ez,acops,jbin)
CM                        print *,'EHB3=',EHB3
                        NOM(k)=NOM(k)+ISTAT
                        NOM(i)=NOM(i)+ISTAT
                        ICNT=ICNT+istat*idist !distance of pairs
                        SUMCT=SUMCT+ISTAT !number of contact pairs
CM                        print *,'Gr2 gt arlm',Gr2,arlm(iseq,kseq),EHB2,
     +                           EHB1b
                        if(Gr2.gt.arlm(iseq,kseq))THEN
                           EHB2=EHB2+apm(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
CM             print *,'cc 2 nom=',i,cc,EHB1a,EHB3,NOM(k),
CM     +                   NOM(i),icnt,SUMCT,EHB2,EHB1b
                  ENDIF
CM                  print *,'end if Cr2 lt 120'
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
CM                 print *,'idist eq 3=',sec(k),sec(i)
***   alpha-helix:
                  if(sec(k).ne.4.and.sec(i).ne.4)then
CM                    print *,'Cr2 lt Cr2a',Cr2,Cr2a,cc,acut_cc
                  if(Cr2.lt.Cr2a.and.cc.gt.acut_cc)then
CM                    print *,'Cr2 lt Cr2a',bxi,bxk,byi,byk,bzi,bzk
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
CM                  print *,'bb=',bb,acut_bb
                  if(bb.gt.acut_bb)then
CM                    print *,'bb gt acut_bb',i,k,mv(i-1),mv(k-1)
                  av11=avv(aax(i-1,mv,x,ex),aay(i-1,mv,y,ey),
     +                     aaz(i-1,mv,z,ez),axi,ayi,azi,
     $                    aax(k-1,mv,x,ex),aay(k-1,mv,y,ey),
     +                    aaz(k-1,mv,z,ez),axk,ayk,azk)
CM                    if(mv(i-1).gt.0)then
CM                      aax1=x(i-1)
CM                      aay1=y(i-1)
CM                      aaz1=z(i-1)
CM                    else
CM                      aax1=ex(i-1)
CM                      aay1=ey(i-1)
CM                      aaz1=ez(i-1)
CM                    endif
CM                    if(mv(k-1).gt.0)then
CM                      aax2=x(k-1)
CM                      aay2=y(k-1)
CM                      aaz2=z(k-1)
CM                    else
CM                      aax2=ex(k-1)
CM                      aay2=ey(k-1)
CM                      aaz2=ez(k-1)
CM                    endif
CM                    av11=avv(aax1,aay1,aaz1,axi,ayi,azi,
CM     +                       aax2,aay2,aaz2,axk,ayk,azk)
CM                    print *,'av11=',i,k,av11,acut_vv,mv(i+1),mv(k+1)
                  if(av11.gt.acut_vv)then
CM                  av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
CM     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                    if(mv(i+1).gt.0)then
                      aax1=x(i+1)
                      aay1=y(i+1)
                      aaz1=z(i+1)
                    else
                      aax1=ex(i+1)
                      aay1=ey(i+1)
                      aaz1=ez(i+1)
                    endif
                    if(mv(k+1).gt.0)then
                      aax2=x(k+1)
                      aay2=y(k+1)
                      aaz2=z(k+1)
                    else
                      aax2=ex(k+1)
                      aay2=ey(k+1)
                      aaz2=ez(k+1)
                    endif
                    av22=avv(axi,ayi,azi,aax1,aay1,aaz1,
     $                    axk,ayk,azk,aax2,aay2,aaz2)
CM                    print *,'av22=',av22,acut_vv,EHB5a
                  if(av22.gt.acut_vv)then
***   ->
                     fact=(1-abs(cc-0.4))*(1-abs(bb-0.815))
                     EHB5a=EHB5a+energyHBa(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk,
     $                    acut_hh,EHBIJ)
                  endif
                  endif
                  endif
                  endif
                  endif
CM                  print *,'if idist 3=',i,idist,Cr2a,acut_cc,bb,acut_bb,
CM     +                     av11,acut_vv,av22,EHB5a
               elseif(idist.gt.4.and.idist.lt.20)then !!!!!!!!!!!!!!!!!!!!!!!
CM                 print *,'idist gt 4 lt 20=',sec(k),sec(i)
***   antiparallel-sheet
                  if(sec(k).ne.2.and.sec(i).ne.2)then
CM                  print *,'Before if Cr2 lt Cr2b',Cr2,Cr2b,cc,bcut_cc,
CM     +                     bxi,bxk,byi,byk,bzi,bzk
CM                    print *,'Cr2 lt Cr2b=',Cr2,Cr2b,cc,bcut_cc
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
CM                  print *,'Cr2 lt Cr2b',Cr2b,bcut_cc,
CM     +                     bxi,bxk,byi,byk,bzi,bzk
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
CM                  print *,'BB=',bb,-bcut_bb
                  if(bb.lt.-bcut_bb)then !antiparallel
CM                    print *,'bb lt -buct_bb=',i,k,mv(i-1),mv(k+1)
CM                  av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
CM     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                    if(mv(i-1).gt.0)then
                      aax1=x(i-1)
                      aay1=y(i-1)
                      aaz1=z(i-1)
                    else
                      aax1=ex(i-1)
                      aay1=ey(i-1)
                      aaz1=ez(i-1)
                    endif
                    if(mv(k+1).gt.0)then
                      aax2=x(k+1)
                      aay2=y(k+1)
                      aaz2=z(k+1)
                    else
                      aax2=ex(k+1)
                      aay2=ey(k+1)
                      aaz2=ez(k+1)
                    endif
                    av12=avv(aax1,aay1,aaz1,axi,ayi,azi,
     $                    axk,ayk,azk,aax2,aay2,aaz2)
CM                    print *,'av12=',av12,-bcut_vv
                  if(av12.lt.-bcut_vv)then
CM                    print *,'av12 lt -bcut_vv',i,k,mv(i+1),mv(k-1)
CM                  av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
CM     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                    if(mv(i+1).gt.0)then
                      aax1=x(i+1)
                      aay1=y(i+1)
                      aaz1=z(i+1)
                    else
                      aax1=ex(i+1)
                      aay1=ey(i+1)
                      aaz1=ez(i+1)
                    endif
                    if(mv(k-1).gt.0)then
                      aax2=x(k-1)
                      aay2=y(k-1)
                      aaz2=z(k-1)
                    else
                      aax2=ex(k-1)
                      aay2=ey(k-1)
                      aaz2=ez(k-1)
                    endif
                  av21=avv(axi,ayi,azi,aax1,aay1,aaz1,
     $                    aax2,aay2,aaz2,axk,ayk,azk)
CM                    print *,'av21=',av21,-bcut_vv,EHB5b
                  if(av21.lt.-bcut_vv)then
***   ->
                     fact=abs(bb)*cc !bb->-1,cc->1
                     EHB5b=EHB5b+energyHBb(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk,
     $                    bcut_hh,EHBIJ)
                  endif
                  endif
                  endif
                  endif
                  endif
CM                print *,'if idist gt 4 lt 20=',i,idist,Cr2b,bcut_cc,bb,
CM     +                     bcut_bb,av12,bcut_vv,av21,EHB5b
               elseif(idist.ge.20)then !!!!!!!!!!!!!!!!!!!!!!!!
CM                 print *,'idist ge 20=',sec(k),sec(i)
                  if(sec(k).ne.2.and.sec(i).ne.2)then
CM                    print *,'Cr2 lt Cr2b',Cr2,Cr2b,cc,bcut_cc
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
CM                  print *,'Cr2 lt Cr2b',Cr2b,bcut_cc,
CM     +                     bxi,bxk,byi,byk,bzi,bzk
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk) 
CM                   print *,'bb=',bb,-bcut_bb
                  if(bb.lt.-bcut_bb)then !antiparallel
***   antiparallel-sheet:
CM                     print *,'bb lt -bcut_bb',i,k,mv(i-1),mv(k+1)
CM                     av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
CM     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                    if(mv(i-1).gt.0)then
                      aax1=x(i-1)
                      aay1=y(i-1)
                      aaz1=z(i-1)
                    else
                      aax1=ex(i-1)
                      aay1=ey(i-1)
                      aaz1=ez(i-1)
                    endif
                    if(mv(k+1).gt.0)then
                      aax2=x(k+1)
                      aay2=y(k+1)
                      aaz2=z(k+1)
                    else
                      aax2=ex(k+1)
                      aay2=ey(k+1)
                      aaz2=ez(k+1)
                    endif
                     av12=avv(aax1,aay1,aaz1,axi,ayi,azi,
     $                    axk,ayk,azk,aax2,aay2,aaz2)
CM                     print *,'av12=',av12,-bcut_vv
                     if(av12.lt.-bcut_vv)then
CM                        print *,'av12 lt -bcut_vv',i,k,mv(i+1),mv(k-1)
CM                     av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
CM     $                       aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                    if(mv(i+1).gt.0)then
                      aax1=x(i+1)
                      aay1=y(i+1)
                      aaz1=z(i+1)
                    else
                      aax1=ex(i+1)
                      aay1=ey(i+1)
                      aaz1=ez(i+1)
                    endif
                    if(mv(k-1).gt.0)then
                      aax2=x(k-1)
                      aay2=y(k-1)
                      aaz2=z(k-1)
                    else
                      aax2=ex(k-1)
                      aay2=ey(k-1)
                      aaz2=ez(k-1)
                    endif
                     av21=avv(axi,ayi,azi,aax1,aay1,aaz1,
     $                       aax2,aay2,aaz2,axk,ayk,azk)
CM                     print *,'av21=',av21,-bcut_vv,EHB5b
                     if(av21.lt.-bcut_vv)then
***   ->
                        fact=abs(bb)*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk,
     $                       bcut_hh,EHBIJ)
                     endif
                     endif
CM                     print *,'if bb 1=',i,bb,bcut_bb,av12,bcut_vv,
CM     +                       av21,EHB5b
                  elseif(bb.gt.bcut_bb)then !bb>0, parallel H-bond
***   parallel-sheet:
CM                     print *,'bb gt b_cut_bb',i,k,mv(i-1),mv(k-1),
CM     +                        b_cutbb
CM                     av11=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
CM     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                    if(mv(i-1).gt.0)then
                      aax1=x(i-1)
                      aay1=y(i-1)
                      aaz1=z(i-1)
                    else
                      aax1=ex(i-1)
                      aay1=ey(i-1)
                      aaz1=ez(i-1)
                    endif
                    if(mv(k-1).gt.0)then
                      aax2=x(k-1)
                      aay2=y(k-1)
                      aaz2=z(k-1)
                    else
                      aax2=ex(k-1)
                      aay2=ey(k-1)
                      aaz2=ez(k-1)
                    endif
                     av11=avv(aax1,aay1,aaz1,axi,ayi,azi,
     $                    aax2,aay2,aaz2,axk,ayk,azk)
CM                     print *,'av11=',av11,bcut_vv
                     if(av11.gt.bcut_vv)then
CM                       print *,'av11 gt bcut_vv',i,k,mv(i+1),mv(k+1)
CM                     av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
CM     $                       axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                    if(mv(i+1).gt.0)then
                      aax1=x(i+1)
                      aay1=y(i+1)
                      aaz1=z(i+1)
                    else
                      aax1=ex(i+1)
                      aay1=ey(i+1)
                      aaz1=ez(i+1)
                    endif
                    if(mv(k+1).gt.0)then
                      aax2=x(k+1)
                      aay2=y(k+1)
                      aaz2=z(k+1)
                    else
                      aax2=ex(k+1)
                      aay2=ey(k+1)
                      aaz2=ez(k+1)
                    endif
                     av22=avv(axi,ayi,azi,aax1,aay1,aaz1,
     $                       axk,ayk,azk,aax2,aay2,aaz2)
CM                      print *,'av22=',av22,bcut_vv,EHB5b
                     if(av22.gt.bcut_vv)then
***   ->
                        fact=bb*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk,
     $                       bcut_hh,EHBIJ)
                     endif
                     endif
CM                     print *,'if bb 2=',i,bb,bcut_bb,av11,bcut_vv,
CM     +                       av22,EHB5b
                  endif
                  endif
                  endif
CM                  print *,'if idist gt 20=',i,idist,Cr2b,bcut_cc,bb,
CM     +                     bcut_bb,av12,bcut_vv,av21,EHB5b
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

      print *,'Final EHB',eh1,eh1a,eh1b,eh1c,eh2,eh3,eh4,eh5a,eh5b,
     +         ehb1,ehb1a,ehb1b,ehb1c,ehb2,ehb3,ehb4,ehb5a,ehb5b
      EHB
     $     =eh1*EHB1            !+1/r of Ca-SG
     $     +eh1a*EHB1a          !+1/r for non-parallel of Ca-Ca
     $     +eh1b*EHB1b          !excluded volumn of SG-SG
     $     +eh1c*EHB1c          !pair-wise potential of SG-SG
     $     +eh2*EHB2            !quarsi3 for SG-SG
     $     +eh3*EHB3            !enhance good piece
     $     +eh4*EHB4            !-1/r for parallel contact of Ca-Ca
     $     +eh5a*EHB5a          !H-bond energy (alpha)
     $     +eh5b*EHB5b          !H-bond energy (beta)

CM      print *,'Out EHB',EHB,eh1c,EHB1c,jjjj,kkkk
CM         do ijk = 0,500
CM            print *,'EHB End jbin=',ijk,jbin(ijk)
CM         enddo

c ^^^^^^^^^^ EHB finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate ei5=acops(i,jbin(r15))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ei5(i,idist,Lch,mv,x,y,z,ex,ey,ez,acops,jbin)
      implicit integer(i-z)
      parameter(ndim=1500)

      dimension mv(ndim)
      dimension x(ndim),y(ndim),z(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim)
      dimension JBIN(0:500),acops(ndim,16)


CM      common/chainm/mv(ndim)
CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim)
CM      common/lengths/Lch,Lch1,Lch2
CM      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)

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
c     2, bury potential of SG.
c     3, distmap and contact restrains.
c     4, bias to (prodicted) protein-like structure, panality on crumpling.
c     5, E13,E15,E15.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION ESHORT(iiii,jjjj,istat,seq,mv,x,y,z,ica,gx,gy,gz,ex,ey,
     +                ez,egx,egy,egz,hh,ee,area,afs,mdis,n_dis,kdis,
     +                dist,dev,mdisl,kdisl,mcom,kcom,concut,
     +                McomCA,KcomCA,aweigCA,csr,asr,bsr,sec,ibin,hbx,
     +                hby,hbz,ebx,eby,ebz,cax,cay,caz,ecx,
     +                ecy,ecz,frga,frgb,
     +                ESHORT2,ESHORT3,ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,
     +                ESHORT5b,ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
     +                ESHORT10,ESHORT11,diold,coold,cex,
     +                cey,cez,acrit,distL,aweig,dist_CA_cut,dis,concut2,
     +                Lch,
     +                es2,er1,er3,er4,es3,es3a,es3b,es3c,es4,es5,es6,
     +                er5,er6,er7,conew,colim,codev,codevsum,
     +                didev,didevsum,dinew,jbin)
      IMPLICIT INTEGER(I-Z)
      parameter(ndim=1500)	
      parameter(nvec=312)

      dimension seq(ndim),sec(ndim) 	  
      dimension hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      dimension cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      dimension IBIN(-300:300),asr(ndim,-12:12)
      dimension JBIN(0:500)
      dimension afs(ndim)
      dimension frga(ndim),frgb(ndim)
      dimension csr(ndim,2)
      dimension gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      dimension Mcom(ndim),Kcom(ndim,100)
      dimension Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      dimension concut(0:19,0:19),concut2(0:19,0:19)
      dimension mv(ndim)
      dimension McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim)
      dimension ex(ndim),ey(ndim),ez(ndim) !Ca
      dimension egx(ndim),egy(ndim),egz(ndim) !SG
      dimension ecx(ndim),ecy(ndim),ecz(ndim) !cc
      dimension ebx(ndim),eby(ndim),ebz(ndim) !Hb
      dimension aweig(4000,4000)
      dimension MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      dimension area(ndim)
      dimension bgx(ndim),bgy(ndim),bgz(ndim)
      dimension AA(3,3),EE(3),HH(3,3)
      dimension bsr(ndim,16)
      dimension acops(ndim,16)




CM      common/lengths/Lch,Lch1,Lch2
CM      common/seqe/seq(ndim),sec(ndim) 	  
CM      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
CM      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
CM      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
CM      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
CM      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
CM      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
CM      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
CM      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
CM      common/fr/frga(ndim),frgb(ndim)
CM      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
CM      common/three/angle(nvec,nvec)
CM      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
CM      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
CM      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
CM      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
CM      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
CM      common/msichores/msicho
CM      common/distres/er4,es3c
CM      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
CM      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT4,ESHORT11
CM      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
CM      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT10
CM      common/hopp/eonehw(0:19)
CM      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
CM      common/chainm/mv(ndim)
CM
CM      common/hom1/asr2(ndim,4),asr3(ndim,4),asr4(ndim,14),asr5(ndim,8)
CM      common/hom2/ibb2(0:999),ibb3(0:999),ibb4(-999:999),ibb5(0:999)
CM      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(4000,4000)

CM      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
CM      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
CM      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
CM      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
CM      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
CM      common/zscore/izscore
CM      common/freg/aweig(4000,4000)
CM      common/CAcontact1/dist_CA_cut
CM      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
CM      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
CM      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)
CM      common/E_defo/i_E_defo
CM      common/initialinput/switch,k_cycle,k_phot,N_ann
CM      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
CM      common/expose/mp(20,ndim),area(ndim)
CM
CM      common/sidechain/bgx(ndim),bgy(ndim),bgz(ndim)
      dimension cgx(ndim),cgy(ndim),cgz(ndim)

CM      common/center/cex,cey,cez
CM      common/eigen/AA(3,3),EE(3),HH(3,3)
      
cccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

      print *,'In ESHORT',istat,iiii,jjjj
      print *,'istat eq 10=',istat,conew,colim,codevsum,
     +              coold,didevsum,dinew,diold

CM         do ijk = 0,500
CM            print *,'ESHORT begin jbin=',ijk,jbin(ijk)
CM         enddo

c      ESHORT=0.0
      ESHORT2=0  !bury potential for SG

      ESHORT3=0  !distance restrain from both threading (for C_a)
      ESHORT4=0  !contact restrain from threading (for SG)
      ESHORT4a=0 !deviation of contact restrain

      ESHORT5=0  !bias2,3: v(i)-v(i+4) anti/parallel; c(i)-c(i+2) anit/paralel
      ESHORT5a=0 !crumpling
      ESHORT5b=0 !bias4 to predicted alpha/beta structure.
      ESHORT5c=0 !bias1 to possible alpha/beta structure. 

      ESHORT6=0  !correlation of E13 of Ca
      ESHORT7=0  !correlation of E14, from both common and 2th specific data
      ESHORT8=0  !correlation of E15, from both common and 2th specific data

      ESHORT9=0  !contact restraints of CA
      ESHORT10=0 !Long-range distance restraints of CA

      ESHORT11=0 !RMSD deviation

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
***** Bury energy of SG --------------------->
         gxn=HH(1,1)*(agxi-cex)+HH(1,2)*(agyi-cey)+HH(1,3)*(agzi-cez)
         gyn=HH(2,1)*(agxi-cex)+HH(2,2)*(agyi-cey)+HH(2,3)*(agzi-cez)
         gzn=HH(3,1)*(agxi-cex)+HH(3,2)*(agyi-cey)+HH(3,3)*(agzi-cez)
         fff=gxn**2/EE(1)+gyn**2/EE(2)+gzn**2/EE(3) !=5 for ellipsoid surphase
         if(fff.lt.3)then
            aaa=fff-3
            if(aaa.lt.-1)aaa=-1
c     area(i)<0, propensity to bury, >0, prop. to be exposed
            ESHORT2=ESHORT2-aaa*area(i)
         endif
c^^^^^^^^^^^^^^^^^^ bury energy finished ^^^^^^^^^^^^^^^^^^

CM         do ijk = 0,500
CM            print *,'ESHORT2 middle 1 jbin=',ijk,jbin(ijk)
CM         enddo

c     deepth factor: afs=1, r<r0; [0.5,1], r0<r<2r0; 0.5, r>2r0
         ar=sqrt((axi-cex)**2+(ayi-cey)**2+(azi-cez)**2+0.01) !r
         afs(i)=acrit/ar
         if(afs(i).gt.1)afs(i)=1
         if(afs(i).lt.0.5)afs(i)=0.5

CM         print *,'afs',i,afs(i)
         
****  restrain from 'dist.dat'---------------------->
         N_dis=Mdis(i)          !number of distance restraints on 'i'
         do k=1,N_dis
            j=kdis(i,k)         !i-j restraints
            if(j.lt.i.or.j.gt.jjjj) then !to avoid repeat
               dij2=(axi-aax(j,mv,x,ex))**2+(ayi-aay(j,mv,y,ey))**2+
     +              (azi-aaz(j,mv,z,ez))**2
               dij=sqrt(dij2)
               err=abs(dij-dist(i,k)) !dist: predicted dis
               if(err.gt.dev(i,k)) then !dev: deviation for arca
                  ESHORT3=ESHORT3+1
               endif
            endif
CM            print *,'ESHORT3=',k,ESHORT3
         enddo
c^^^^^^^^^^^^^^^^^ distant restrain finished ^^^^^^^^^^^^^^^^^^^^^^^^^

CM         do ijk = 0,500
CM            print *,'ESHORT3 middle 1 jbin=',ijk,jbin(ijk)
CM         enddo

****  long-range dist-restrain from 'distL.dat'---------------------->
         N_disL=MdisL(i)        !number of long-range dist restraints on 'i'
         do k=1,N_disL
            j=kdisL(i,k)        !i-j restraints
            if(j.lt.i.or.j.gt.jjjj) then !to avoid repeat
               dij2=(axi-aax(j,mv,x,ex))**2+(ayi-aay(j,mv,y,ey))**2+
     +              (azi-aaz(j,mv,z,ez))**2
               dij=sqrt(dij2)
               err=abs(dij-distL(i,k))
               if(err.lt.1)err=1
               ESHORT10=ESHORT10-1/err
            endif
CM            print *,'ESHORT10=',k,ESHORT10
         enddo
c^^^^^^^^^^^^^^^^^distant restrain finished ^^^^^^^^^^^^^^^^^^^^^^^^^
CM         do ijk = 0,500
CM            print *,'ESHORT10 middle 1 jbin=',ijk,jbin(ijk)
CM         enddo

 1    continue
      
CM         do ijk = 0,500
CM            print *,'ESHORT middle 1 jbin=',ijk,jbin(ijk)
CM         enddo

ccc   Contact restrain from 'comb.dat' -------------------------------->
       do 11 i=iiii,jjjj
          N_com=Mcom(i)         !number of contacts on i
          if(N_com.ge.1)then
             iseq=seq(i)
             if(mv(i).gt.0)then
                agxi=x(i)+GX(ica(i-1),ica(i),iseq) !SG
                agyi=y(i)+GY(ica(i-1),ica(i),iseq) !SG
                agzi=z(i)+GZ(ica(i-1),ica(i),iseq) !SG
             else
                agxi=egx(i)     !SG
                agyi=egy(i)
                agzi=egz(i)
             endif
             
             do k=1,N_com
                j=Kcom(i,k)     !k'th contact with i
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
CM                   print *,'aij2=',k,j,aij2,concut2(iseq,jseq)
                   
                   if(aij2.gt.concut2(iseq,jseq))then
                      ESHORT4=ESHORT4+2*aweig(i,j) !no contact
                      if(istat.lt.0)then
                         coold=coold+sqrt(aij2)-concut(iseq,jseq) !penality
                      else
                         conew=conew+sqrt(aij2)-concut(iseq,jseq) !panelity
                      endif
                   endif
                endif
             enddo
CM             print *,'ESHORT4=',ESHORT4
          endif
 11    continue
c^^^^^^^^^^^^^^^^^^ contact restrains finished ^^^^^^^^^^^^^^^^^^^^^

CM         do ijk = 0,500
CM            print *,'ESHORT middle 11 jbin=',ijk,jbin(ijk)
CM         enddo

ccc   Contact restrain from 'combCA.dat' ------------------------------>
c      dist_CA_cut=(6.0/0.87)**2
      do 33 i=iiii,jjjj
         N_com=McomCA(i)        !number of contacts on i
         if(N_com.ge.1)then
            if(mv(i).gt.0)then
               axi=x(i)
               ayi=y(i)
               azi=z(i)
            else
               axi=ex(i)
               ayi=ey(i)
               azi=ez(i)
            endif

            do k=1,N_com
               j=KcomCA(i,k)    !k'th contact with i
               if(j.lt.i.OR.j.gt.jjjj)then
                  aij2=(axi-aax(j,mv,x,ex))**2+(ayi-aay(j,mv,y,ey))**2+
     +                 (azi-aaz(j,mv,z,ez))**2
CM                  print *,'aij2=',i,j,aij2,aweigCA(i,j),dist_CA_cut
                  if(aij2.lt.dist_CA_cut)then
                     ESHORT9=ESHORT9-aweigCA(i,j)
                  endif
               endif
CM               print *,'ESHORT9=',k,j,ESHORT9
            enddo
         endif
 33   continue
c^^^^^^^^^^^^^^^^^^ CAcontact restrains finished ^^^^^^^^^^^^^^^^^^^^^

CM         do ijk = 0,500
CM            print *,'ESHORT middle 33 jbin=',ijk,jbin(ijk)
CM         enddo


*********E13 --------------------------------------->

CM      print *,'Here in ESHORT ha ha ha'
      i1=max(iiii-1,1)
      i2=min(jjjj-1,Lch-2)
CM      print *,'i1 i2',i1,i2
      do i=i1,i2
         ar13=(aax(i,mv,x,ex)-aax(i+2,mv,x,ex))**2+
     +        (aay(i,mv,y,ey)-aay(i+2,mv,y,ey))**2+
     $        (aaz(i,mv,z,ez)-aaz(i+2,mv,z,ez))**2
         if(ar13.lt.48) then    !6.03A
            ESHORT6=ESHORT6+csr(i,1)
         else
            ESHORT6=ESHORT6+csr(i,2)
         endif
CM         print *,'ESHORT6=',i,ESHORT6
      enddo
c^^^^^^^^^^^^^^ E_13 finished ^^^^^^^^^^^^^^^^^^^^^^^^^

CM         do ijk = 0,500
CM            print *,'ESHORT middle i1 i2 jbin=',ijk,jbin(ijk)
CM         enddo

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
CM         print *,'ESHORT7=',r14,hand,ESHORT7
cccccccccE15:
CM         do ijk = 0,500
CM            print *,'ESHORT middle jbin=',ijk,jbin(ijk)
CM         enddo
         ar15=(ax1-ax5)**2+(ay1-ay5)**2+(az1-az5)**2
         r15=nint(ar15)
         print *,'r15=',r15
         if(r15.gt.500) r15=500
CM         do ijk = 0,500
CM            print *,'ESHORT jbin=',ijk,jbin(ijk)
CM         enddo
CM         ijk = JBIN(r15)
CM         print *,'ijk=',ijk
CM         print *,'bsr=',i,r15,jbin(r15)
         ESHORT8=ESHORT8+bsr(i,jbin(r15))
CM         print *,'ESHORT8=',r15,ESHORT8

cccccccccbias1: encourage helix/sheet-like structure
         if(ar15.lt.75)THEN     !7.53A: helix
CM            print *,'ar15 lt 75',ar15,sec(i+1),sec(i+2),sec(i+3)
            if(sec(i+1).ne.4.AND.sec(i+2).ne.4.AND.sec(i+3).ne.4)then
               if(ibin(r14).gt.4.AND.ibin(r14).lt.8)then !right-hand,3A->8A
                  dot13=(ax4-ax3)*(ax2-ax1)
     $                 +(ay4-ay3)*(ay2-ay1)+(az4-az3)*(az2-az1)
CM                  print*,'dot13=',dot13
                  if(dot13.lt.0)then
                     dot24=(ax5-ax4)*(ax3-ax2)
     $                    +(ay5-ay4)*(ay3-ay2)+(az5-az4)*(az3-az2)
CM                     print*,'dot24=',dot24
                     if(dot24.lt.0)then
                        dot14=(ax5-ax4)*(ax2-ax1)
     $                       +(ay5-ay4)*(ay2-ay1)+(az5-az4)*(az2-az1)
CM                        print*,'dot14=',dot14
                        if(dot14.gt.0)then
CM                           print *,'if dot14 gt 0=',i,afs(i),afs(i+1),
CM     +                             afs(i+3),afs(i+4)
                           ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
                           ESHORT5c=ESHORT5c-2-ff
                        endif
                     endif
                  endif
               endif
            endif
         elseif(ar15.gt.160)then !11A, beta
CM            print *,'ar15 gt 160=',ar15,sec(i+1),sec(i+2),sec(i+3)
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
CM               print*,'dot24=',dot24
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
CM                  print*,'dot23=',dot23
                  if(dot23.lt.-0.71)then !angle >135
CM                     print *,'if dot23 lt 0.71=',i,afs(i),afs(i+1),
CM     +                             afs(i+3),afs(i+4)
                     ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
c                     ESHORT5c=ESHORT5c-2-ff
                     ESHORT5c=ESHORT5c-2-ff*2
                  endif
               endif
            endif
         endif
CM         print *,'ESHORT5c=',i,ESHORT5c
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
CM         print *,'b15=',i,b15,seq(i),seq(i+2),seq(i+4)
         if(seq(i).eq.2.and.seq(i+2).eq.2.and.seq(i+4).eq.2)then
            if(b15.gt.0.9)then
               ESHORT5=ESHORT5-ff !alpha
            endif
         else
            if(b15.gt.0.5.or.b15.lt.-0.3)then !beta or turn
c               ESHORT5=ESHORT5-ff
               ESHORT5=ESHORT5-ff*2
            endif
         endif
CM         print *,'ESHORT5=',ESHORT5
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
CM         print *,'c13=',c13
         c13=min(0.71,c13)/0.71 !c13 is the same in [0,45]
         ESHORT5=ESHORT5-ff*c13
CM         print *,'ESHORT5=',ESHORT5
      enddo
c^^^^^^^^^^^^^^^^^^ bias3 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccccccbias4a: to predicted alpha fragment --------------->
      i1=max(iiii-6,1)
      i2=min(jjjj-1,Lch-7)
      do i=i1,i2
         if(frga(i).gt.1)then
            dis=sqrt((aax(i,mv,x,ex)-aax(i+7,mv,x,ex))**2+
     $           (aay(i,mv,y,ey)-aay(i+7,mv,y,ey))**2+
     +               (aaz(i,mv,z,ez)-aaz(i+7,mv,z,ez))**2)
            ESHORT5b=ESHORT5b+abs(dis-frga(i))
         endif
CM         print *,'ESHORT5b=',ESHORT5b
      enddo
ccccccbias4b: to predicted beta fragment --------------->
      i1=max(iiii-5,1)
      i2=min(jjjj-1,Lch-6)
      do i=i1,i2
         if(frgb(i).gt.1)then
            dis=sqrt((aax(i,mv,x,ex)-aax(i+6,mv,x,ex))**2+
     $           (aay(i,mv,y,ey)-aay(i+6,mv,y,ey))**2+
     +               (aaz(i,mv,z,ez)-aaz(i+6,mv,z,ez))**2)
c            ESHORT5b=ESHORT5b+abs(dis-frgb(i))
            ESHORT5b=ESHORT5b+abs(dis-frgb(i))*2
         endif
CM         print *,'ESHORT5b=',ESHORT5b
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
CM         print *,'ESHORT5a=',ESHORT5a
      enddo
c^^^^^^^^^^^^^ penality of bizard structure finished (ESC1) ^^^^^^^^^^^^

c     Further penalize deriviation of restrain, if larger than (colim, dilim):
CM      print *,'istat eq 10=',istat,conew,colim,codevsum,
CM     +              coold,didevsum,dinew,diold

      IF(istat.eq.10) then      !calculate E from the beginning
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
         if(codevsum.gt.colim) ESHORT4a=codevsum-colim
c     total penalty-energy beacuse of restrain-deviation on new conformation:
         if(codev.gt.colim) ESHORT4a=ESHORT4a-(codev-colim)
      endif

      print *,'Final ESHORT=',es2,er1,er3,er4,es3,es3a,es3b,es3c,
     +         es4,es5,es6,er5,er6,er7,ESHORT2,ESHORT3,ESHORT4,
     +         ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c,
     +         ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT10,ESHORT11

      ESHORT=
     $     +es2*ESHORT2
     $     +er1*ESHORT3
     $     +er3*ESHORT4
     $     +er4*ESHORT4a
     $     +es3*ESHORT5
     $     +es3a*ESHORT5a
     $     +es3b*ESHORT5b
     $     +es3c*ESHORT5c
     $     +es4*ESHORT6
     $     +es5*ESHORT7
     $     +es6*ESHORT8
     $     +er5*ESHORT9
     $     +er6*ESHORT10
     $     +er7*ESHORT11

c ^^^^^^^^^^ E_short finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END


ccccccccccccccccweight factor cccccccccccccccccccccccccccccccccccccccccccc
      function weight12(dE,atemp,ichos,aTs)
      parameter(nrep=100)
CM      common/temperature/itemp,atemp
CM      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
CM      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
CM      common/ichos/ichos

CM      print *,'In weight12 dE=',dE,atemp,ichos,aTs

      if(ichos.eq.1)then
         weight12=exp(-dE/atemp)
         ichos=2
      else
CM Removed the square2 function and added line in this function

         weight12=exp(-(dE*dE)/aTs)
         ichos=1
      endif
CM      print *,'Out WT=',ichos,weight12

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

CM      print *,'In avv',avv
c^^^^^^^^^^^^^ v.v finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     H-bond for alpha-helix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energyHBa(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk,
     +                  acut_hh,EHBIJ)
      implicit integer(i-z)
      parameter(ndim=1500)      !maximum length of chain-length
      parameter(nvec=312)       !number of vectors

      dimension EHBIJ(ndim,ndim)

CM      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
CM      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
CM      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh

c     because helix is right-hand, bi=v(i-1)(x)v(i) is always
c     same direction as v(i---k), if (k>i);
c     reverse direction as v(i---k), if (k<i);

      print *,'In energyHBa',acut_hh,EHBIJ(i,k)

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

      print *,'energyHBa=',energyHBa

c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     H-bond for beta-sheet
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energyHBb(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk,
     +                   bcut_hh,EHBIJ)
      implicit integer(i-z)
      parameter(ndim=1500)      !maximum length of chain-length
      parameter(nvec=312)       !number of vectors

      dimension EHBIJ(ndim,ndim)

CM      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
CM      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
CM      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh

      print *,'In energyHBb',bcut_hh,EHBIJ(i,k)

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

      print *,'energyHBb=',energyHBb
c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

