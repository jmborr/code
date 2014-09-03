ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Compile for debugging:
c     pgf77 -Mextex -Bstatic -g -Mbounds -Mchkfpstk -o spkicker spicker.f
c     Compile for production:
c     pgf77 -Mextend -O3 -Bstatic -o spkicker spicker.f
c
c     This program, Structure-PCIKER (SPICKER), is aimed at selecting the 
c     best fold by clustering trajectors produced by CABS model. 
c     Comments and bug reports should be addressed to zhang6@buffalo.edu.
c
c     Input files includes:
c       'rmsinp'---Mandatory, length of protein & piece for RMSD calculation;
c       'seq.dat'--Mandatory, sequence file, for output of PDB models.
c       'tra.in'---Mandatory, list of trajectory names used for clustering.
c       files in 'tra.in'---Mandatory, trajectories of structures.
c       'chain.dat'---Mandatory, containing the template coordinate files
c       'exp.dat'---Mandatory, containing predicted buried surface
c       'CA'-------Optional, native structure, for comparison to native.
c       'TEMP'-----Optional, template file, for structure close to template.
c
c     Output files includes:
c       'rst.dat'-----summary of clustering results;
c       'str.txt'-----list of structure in cluster;
c       'combo*.pdb'--PDB format of cluster centroids;
c       'closc*.pdb'--PDB format of structures closest to centroids;
c       'templ.pdb'---PDB format of structure closest to given template;
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        1         2         3         4         5         6         7 !
c 345678901234567890123456789012345678901234567890123456789012345678901234567
      program cluster
      parameter(ndim=1000)      !Length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(ntr=20000)      !number of trajectories
      parameter(ncl=100)        !number of clusters
      parameter(minloopL=4)     !minimum loop length to consider

      character *256 line       !buffer
      integer nline             !lenght of buffer
      common/whole/Lch,x,y,z
      integer Lch
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures

      common/native/Lch_n,x_n,y_n,z_n
      integer Lch_n
      dimension x_n(ndim),y_n(ndim),z_n(ndim) !native structure

      common/clusteriza/nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,
     $     i_str_cl_ex,rmsd_cl_cut
      integer nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,i_str_cl_ex
      dimension i_cl(ncl) !i_cl(i) is structure center of i'th cluster
      dimension n_str_cl(ncl)!of structures including some used struct in i'th cluster
      dimension n_str_cl_ex(ncl) !n_str_cl_ex(i): # of structures excluding some used struct
                                 ! in i'th cluster
      dimension i_str_cl(ncl,nst) !i_str_cl(i,j) is structure id of j'th structure belonging
                                  !to i'th cluster
      dimension i_str_cl_ex(ncl,nst) !i_str_cl_ex(i,j): structure id of j'th structure
                                     !belonging to i'th cluster
      real rmsd_cl_cut
      dimension rmsd_cl_cut(ncl) !rmsd_cl_cut(i)

c     variables related to exp.dat and seq.dat
      common/expdat/expd,hsc
      real expd(ndim)
      integer hsc(ndim)

c     common/chaindat/ contains variables related to the first template
c     of chain.dat file:
c     isalig(i): residue "i" has coordinates in template
c     nloop: number of template-free regions
c     lbegin(i): position where loop "i" begins
c     lextent(i): length of loop "i"
c     lextents(i): length of extended loop "i"
c     lresid(j,i): residue "j" belongs to extended loop "i"
      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      integer isalig(0:ndim),nloop,lbegin(ndim),lextent(ndim),
     &     lextents(ndim),lresid(ndim,ndim)

c     variables related to loop-derived properties
      common/proploop/lrmsdl,lrmsds,lexp,lhelix,lcoil,lstrand,
     &nrmsdl,nrmsds,lrmsdns
      common/names/protein
      character protein*10

      common/trafiles/filen,n_tra
      character filen(ntr)*72
      integer n_tra

      character c2*6,seq(ndim)*3
      character txt1*70,txt2*70
ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /1000*1.0/
ccc   

      real lrmsdl(ndim),lrmsds(ndim),lexp(ndim),
     &lhelix(ndim),lcoil(ndim),lstrand(ndim),
     &nrmsdl(ndim),nrmsds(ndim),lrmsdns(ndim)

      dimension xtemp(ndim),ytemp(ndim),ztemp(ndim) !structure close to templ
      dimension xt(ndim),yt(ndim),zt(ndim) !temporal coordinate
      dimension amat(nst,nst)   !RMSD matrics
      dimension mark(nst)       !for removing used structures
      dimension n_str_near(nst) !number of neighboring structures
      dimension itra(nst),istr(nst),E(nst)
      dimension E_combo(ncl)
      dimension rmsd_str(nst)
      dimension xs(ndim),ys(ndim),zs(ndim)

      common/globalclusters/nglobal,xc,yc,zc
      integer nglobal
      dimension xc(ncl,ndim),yc(ncl,ndim),zc(ncl,ndim) !combo structures

      common/loop_pdb/loopxc,loopyc,loopzc,
     $     loopxclosc,loopyclosc,loopzclosc,rmsd_combo_closc
      real loopxc,loopyc,loopzc,loopxclosc,loopyclosc,loopzclosc,
     $     rmsd_combo_closc
      dimension loopxc(ncl,ncl,ndim),loopyc(ncl,ncl,ndim), !combo str. for unaligned reg
     $     loopzc(ncl,ncl,ndim) !loopxc(i,j,k) residue k of combo of cluster j for unaligned reg. i
      dimension loopxclosc(ncl,ncl,ndim),loopyclosc(ncl,ncl,ndim), !closc str. for unaligned reg
     $     loopzclosc(ncl,ncl,ndim)
      dimension rmsd_combo_closc(ncl) !rmsd between combo and closc

      dimension xcl(ncl,ndim),ycl(ncl,ndim),zcl(ncl,ndim)
      dimension xc3(ncl,ndim),yc3(ncl,ndim),zc3(ncl,ndim)
      dimension rmsd_close_min(ncl) !RMSD between 'combo.pdb' and 'closm.pdb'
      dimension i1(100),i2(100)
      integer q(1000)
      common/r_/R_in,R_ex,Rc_in,Rc_ex
      dimension R_in(100),R_ex(100),Rc_in(100),Rc_ex(100)

      common/dens/n_str_all,n_str,d_rc_in,drcinmax,cmax,nstr_selected,
     $     str_selected
      real d_rc_in(100),drcinmax
      integer n_str_all,n_str,cmax,nstr_selected,str_selected(nst)

      dimension order(nst)
      dimension ires(ndim)

**********************************************************************
****  input:
      open(1,file='rmsinp',status='old') !read Lch
      open(2,file='seq.dat',status='old') !read SEQ for model.pdb
      open(3,file='CA',status='unknown') !for comparison to native
      open(4,file='tra.in',status='old') !information for trajectories
      open(5,file='TEMP',status='unknown') !for comparison to template
      open(7,file='chain.dat',status='old')!for retrieval of unaligned residues
      open(8,file='exp.dat',status='old')  !for predicted solvent exposed residues
****  output:
      open(20,file='rst.dat',status='unknown')
      open(21,file='str.txt',status='unknown') !structures in cluster
      open(22,file='RMSD.list',status='unknown') !structures in cluster
**********************************************************************

      read(1,*)n1,n2            !calculate RMSD between [n1,n2]
      read(1,*)Lch
      read(1,*)protein

********************************************************************

c     read exposed surface, initializing variable exp
      call readexpdat(Lch)

c     read chain.dat, initializing isalig,lbegin,lextent,nloop
      itemp=1 !read first template by default
      call readchaindat(Lch,itemp)

******** read template file ***************************
      L_temp=0
      do i=1,Lch
         q(i)=0
      enddo
      read(5,*,end=1)L_temp
      k=0
 1    do i=1,L_temp
         read(5,1239)tx,ink,tx,tx,ii,tx,a1,a2,a3
         q(ii)=1
         k=k+1
         r_1(1,k)=a1
         r_1(2,k)=a2
         r_1(3,k)=a3
      enddo
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

***   find out the total number of structures & structure closest to template:
***   The total number of structures will be used for picking clustering struct.
      R_temp_min=10000
      read(4,*)n_tra
      i_str_all=0
      do i_tra=1,n_tra
         i_str_tra=0
         read(4,'72A')filen(i_tra)
         open(10,file=filen(i_tra),status='old')
 10      read(10,*,end=19,err=19)Lch1,energ
         do i=1,Lch1
            read(10,*,end=19,err=19)xt(i),yt(i),zt(i)
         enddo
         i_str_tra=i_str_tra+1
         i_str_all=i_str_all+1
***   find structure closest to template:
         if(L_temp.gt.1)then
            k=0
            do i=1,Lch
               if(q(i).eq.1)then
                  k=k+1
                  r_2(1,k)=xt(i)
                  r_2(2,k)=yt(i)
                  r_2(3,k)=zt(i)
               endif
            enddo
            call u3b(w,r_1,r_2,L_temp,0,rms,u,t,ier)
            armsd=dsqrt(rms/L_temp)
            if(armsd.lt.R_temp_min)then
               R_temp_min=armsd
               do i=1,Lch
                  xtemp(i)=xt(i)
                  ytemp(i)=yt(i)
                  ztemp(i)=zt(i)
               enddo
               i_str_temp=i_str_tra
               i_tra_temp=i_tra
            endif
         endif
***
         goto 10
 19      close(10)
      enddo
      n_str_all=i_str_all
c      write(*,*)'Total number structures=',n_str_all
c^^^^^^^^^^^^^^^^^^^^^^^^ n_str_all done ^^^^^^^^^^^^^^^^^^^^^^^^

cccccccccccc read native structure cccccccccccccccccccc
      i=0
 11   read(3,'a',end=5) line
      if(line(13:16).ne." CA ") goto 11
      if( (line(1:3).eq."TER").or.(line(1:3).eq."END")) goto 5
      read(line,1239)tx,ii,tx,tx,ii,tx,a1,a2,a3
      if(ii.ge.1.and.ii.le.Lch)then
         i=i+1
         ires(i)=ii
         x_n(i)=a1
         y_n(i)=a2
         z_n(i)=a3
         r_1(1,i)=a1
         r_1(2,i)=a2
         r_1(3,i)=a3
      endif
      goto 11
 5    continue
      Lch_n=i                   !length of native
 1239 format(A4,I7,A4,A7,I4,A4,3F8.3)
c^^^^^^^^^^^^^^^ read native done ^^^^^^^^^^^^^^^^^^^^^^

ccc read structures from trajectories and find best structure in 13000 cccccccc
      i_str_tra_a=0 !<i_str_tra>
      rmsd_min_all=100
      rmsd_min_use=100
      n_max=nst                 !maximum structures to handle
      delta=float(n_str_all)/float(n_max) !take a structure in each delta
      if(delta.lt.1)delta=1
      i_str=1			!number of structure used for clustering
      i_str_all=0
      E_min=100000
      E_min_all=100000
      write(22,*)'RMSD   Energy   #str   trajectory'
      do 100 i_tra=1,n_tra
         open(10,file=filen(i_tra),status='old')
         i_str_tra=0            !order number of the structure in this file
 20      read(10,*,end=29,err=29)Lch1,energ
         do i=1,Lch1
            read(10,*,end=29,err=29)xt(i),yt(i),zt(i)
         enddo
         i_str_tra=i_str_tra+1
         i_str_all=i_str_all+1
***   calculate RMSD of structure to native----------->
         if(Lch_n.gt.1)then
            do i=1,Lch_n
               r_2(1,i)=xt(ires(i))
               r_2(2,i)=yt(ires(i))
               r_2(3,i)=zt(ires(i))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n)
            if(rmsd_min_all.gt.armsd)then
               rmsd_min_all=armsd
               E_rma=energ
               itra_rmsd_min_all=i_tra
               istr_rmsd_min_all=i_str_tra
            endif
            if(E_min_all.gt.energ)then
               rmsd_E_min_all=armsd
               E_min_all=energ
               itra_E_min_all=i_tra
               istr_E_min_all=i_str_tra
            endif
         endif
***   select structure for clustering-------->
         if(i_str_all.ge.i_str*delta)then
           if(i_str.gt.n_max) goto 29
            order(i_str_tra)=order(i_str_tra)+1
            i_str_tra_a=i_str_tra_a+i_str_tra
            itra(i_str)=i_tra
            istr(i_str)=i_str_tra
            E(i_str)=energ
            n_str=i_str
            do i=1,Lch
               x(i,i_str)=xt(i)
               y(i,i_str)=yt(i)
               z(i,i_str)=zt(i)
            enddo
            if(Lch_n.gt.1)then
               if(rmsd_min_use.gt.armsd)then
                  rmsd_min_use=armsd
                  E_rmu=energ
                  itra_rmsd_min_use=i_tra
                  istr_rmsd_min_use=i_str_tra
               endif
               if(E_min.gt.energ)then
                  E_min=energ
                  rmsd_E_min=armsd
                  itra_E_min=i_tra
                  istr_E_min=i_str_tra
               endif
               rmsd_str(i_str)=armsd
               write(22,61)armsd,energ,i_str_tra,filen(i_tra)
 61            format(f8.3,1x,f11.2,1x,i10,1x,a18)
            else
               rmsd_str(i_str)=-1
            endif
            i_str=i_str+1
         endif
***
         goto 20
 29      close(10)
 100  continue
      close(22)

c^^^^^^^^^^^^^^^^^read structures finished ^^^^^^^^^^^^^^^^^^^^
c      write(*,*)i_str,n_str,itra(n_str),istr(n_str),n_max

ccc      output distribution of i_str_tra ------->
      i_str_tra_a=float(i_str_tra_a)/float(n_str)
      write(21,*)'delta=',delta
      write(21,*)'<i_str_tra>=',i_str_tra_a
      write(21,*)'distribution------->'
      do i=1,nst
        if(order(i).gt.0.5)then
          write(21,*)i,order(i)
        endif
      enddo


c     subroutine clusterize performs clustering of a particular
c     unaligned region or the whole sequence.
      call select_structures(.false.)
      call clusterize(-1,amat,.false.)

*****************************************************************
c     calculate 'combo.pdb' and E_combo
COMBO************************************************************
      nglobal=nc
      do i=1,nc
         ns=0
         E_combo(i)=0
         rmsd_nat_a=0
         rmsd_cen_a=0
         do ii=1,Lch
            xs(ii)=0
            ys(ii)=0
            zs(ii)=0
         enddo
***
         k=i_cl(i)              !center structure of the cluster
         do ii=1,Lch
            r_2(1,ii)=x(ii,k)   !center structure
            r_2(2,ii)=y(ii,k)
            r_2(3,ii)=z(ii,k)
         enddo
         nn=n_str_cl(i)         !number of structures in cluster (including)
***
         write(21,1200)i,rmsd_str(k),istr(k),filen(itra(k))
 1200    format('#Cluster',i3,f8.3,i8,2x,a15)
         write(21,*)'i_cl   i_str  R_nat   R_cen  E    #str     traj'
         write(21,*)'-----------------------------------------------'
         write(21,*)'Nstr=',nn
         do l=1,nn
            m=i_str_cl(i,l)     !l'th structure in i'th cluster
            write(21,1201)i,l,rmsd_str(m),amat(k,m),E(m),
     $           istr(m),filen(itra(m))
 1201       format(i3,i8,2f7.3,f9.1,i8,2x,a15)
***   E_combo:
            E_combo(i)=E_combo(i)+E(m)
            ns=ns+1
            rmsd_nat_a=rmsd_nat_a+rmsd_str(m)
            rmsd_cen_a=rmsd_cen_a+amat(k,m)
***   rotate nearby structure into the center structure---->
            do n=1,Lch
               r_1(1,n)=x(n,m)
               r_1(2,n)=y(n,m)
               r_1(3,n)=z(n,m)
            enddo
            call u3b(w,r_1,r_2,Lch,1,rms,u,t,ier) !u rotate r_1 to r_2
            do j=1,Lch
               xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
               yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
               zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
               xs(j)=xs(j)+xx
               ys(j)=ys(j)+yy
               zs(j)=zs(j)+zz
            enddo
         enddo
***   averaging:
         do j=1,Lch
            xc(i,j)=xs(j)/float(ns)
            yc(i,j)=ys(j)/float(ns)
            zc(i,j)=zs(j)/float(ns)
         enddo
         E_combo(i)=E_combo(i)/float(ns)
         rmsd_nat_a=rmsd_nat_a/ns
         rmsd_cen_a=rmsd_cen_a/ns
         write(21,*)'-----------------------------------------'
         write(21,1222)rmsd_nat_a,rmsd_cen_a,E_combo(i)
 1222    format('   average=',2f7.3,f9.1)
      enddo
c^^^^^^^^^^^^^^^^^^^^ combo finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
c     calculate <RMSD to centroid>
*****************************************************************
      do i=1,nc
         do k=1,Lch !store center cluster coordinates in temporary array
            r_1(1,k)=xc(i,k)
            r_1(2,k)=yc(i,k)
            r_1(3,k)=zc(i,k)
         enddo
         nn=0
         Rc_in(i)=0         !RMSD to centroid including
         do j=1,n_str_cl(i) !n_str_cl(i): number structures belonging to cluster i
            m=i_str_cl(i,j) !ID of one of the structures
            do k=1,Lch      !store structure coordinates in temporary array
               r_2(1,k)=x(k,m)
               r_2(2,k)=y(k,m)
               r_2(3,k)=z(k,m)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch) !RMSD12
            Rc_in(i)=Rc_in(i)+armsd
            nn=nn+1
         enddo
         Rc_in(i)=Rc_in(i)/nn
         d_rc_in(i)=n_str_cl(i)/Rc_in(i) !related density
      enddo
      do i=1,nc
         do k=1,Lch
            r_1(1,k)=xc(i,k)
            r_1(2,k)=yc(i,k)
            r_1(3,k)=zc(i,k)
         enddo
         nn=0
         Rc_ex(i)=0             !RMSD to centroid excluding
         do j=1,n_str_cl_ex(i)
            m=i_str_cl_ex(i,j)
            do k=1,Lch
               r_2(1,k)=x(k,m)
               r_2(2,k)=y(k,m)
               r_2(3,k)=z(k,m)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch) !RMSD12
            Rc_ex(i)=Rc_ex(i)+armsd
            nn=nn+1
         enddo
         Rc_ex(i)=Rc_ex(i)/nn
      enddo
c     read sequence (seq) and secondary structure prediction (hsc)
      do i=1,Lch
         read(2,*)ii,seq(i),hsc(i)
      enddo
c     determine the cluster with maximum density
      cmax=1              !initialize cluster with maximum Rc_in density
      drcinmax=d_rc_in(1) !initialize density of cluster with maximum density
      do i=2,nc
         if(d_rc_in(i).gt.drcinmax) then
            cmax=i
            drcinmax=d_rc_in(i)
         endif
      enddo
c      write(*,'ai2')'cmax=',cmax

*****************************************************************
c     calculate 'combc.pdb': combining structure with 5 pieces:
COMBO************************************************************
      ncom=5                    !number of pieces for averaging structures.
      k=Lch/ncom
      do i=1,ncom
         i1(i)=k*(i-1)+1
         i2(i)=k*i
      enddo
      i2(ncom)=Lch
      do i=1,nc
         k=i_cl(i)              !center structure of the cluster
         nn=n_str_cl(i)         !number of structures in cluster (including)
         do i_ncom=1,ncom
            m1=i1(i_ncom)       !starting point of the piece
            m2=i2(i_ncom)       !ending point of the piece
            do ii=m1,m2
               xs(ii)=0
               ys(ii)=0
               zs(ii)=0
            enddo
***
            mm=0
            do ii=m1,m2
               mm=mm+1
               r_2(1,mm)=x(ii,k) !center structure
               r_2(2,mm)=y(ii,k)
               r_2(3,mm)=z(ii,k)
            enddo
***
            ns=0
            do l=1,nn
               ns=ns+1
               m=i_str_cl(i,l)  !l'th structure in i'th cluster
***   rotate nearby structure into the center structure---->
               mm=0
               do n=m1,m2
                  mm=mm+1
                  r_1(1,mm)=x(n,m)
                  r_1(2,mm)=y(n,m)
                  r_1(3,mm)=z(n,m)
               enddo
               call u3b(w,r_1,r_2,mm,1,rms,u,t,ier) !u rotate r_1 to r_2
               j=0
            do jj=m1,m2
               j=j+1
               xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
               yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
               zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
               xs(jj)=xs(jj)+xx
               ys(jj)=ys(jj)+yy
               zs(jj)=zs(jj)+zz
            enddo
            enddo
***   averaging:
            do j=m1,m2
               xc3(i,j)=xs(j)/float(ns)
               yc3(i,j)=ys(j)/float(ns)
               zc3(i,j)=zs(j)/float(ns)
            enddo
         enddo
      enddo
c^^^^^^^^^^^^^^^^^^^^ combc finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
c     calculate 'closc.pdb', closest structure to combo
CLOSC************************************************************
      do i=1,nc
         rmsd_close_min(i)=100.
      enddo
      do i_tra=1,n_tra
         open(10,file=filen(i_tra),status='old')
 110     read(10,*,end=119,err=119)Lch1,energ
         do i=1,Lch1
            read(10,*,end=119,err=119)xt(i),yt(i),zt(i)
         enddo
         do i=1,Lch
            r_1(1,i)=xt(i)
            r_1(2,i)=yt(i)
            r_1(3,i)=zt(i)
         enddo
***   check whether this structure close to 'combo.pdb':
         do j=1,nc
            do i=1,Lch
               r_2(1,i)=xc(j,i)
               r_2(2,i)=yc(j,i)
               r_2(3,i)=zc(j,i)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch)
            if(armsd.lt.rmsd_close_min(j))then
               rmsd_close_min(j)=armsd
               do i=1,Lch
                  xcl(j,i)=r_1(1,i)
                  ycl(j,i)=r_1(2,i)
                  zcl(j,i)=r_1(3,i)
               enddo
            endif
         enddo
***   
         goto 110
 119     close(10)
      enddo
      
c^^^^^^^^^^^^^^^^^^^^ closc finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c     call subroutine to calculate some properties of the unaliged
c     regions with average statistics on member of cluster cmax
c     cmax: cluster ID of biggest cluster
c     i_cl(cmax): structure ID of center structure of cluster with
c     cluster ID equal to cmax
c     n_str_cl(cmax): number of structures in cluster with ID cmax
c     i_str_cl(i,l): structure with ID==l in cluster with ID "i"
c     x,y,z(k,m): coordinates of residue "k" in structure with ID==m
c     xcl,ycl,zcl(j,i): coordinates or residue "j" in closc structure of
c     cluster with ID==i
c      write(*,'i4')i_cl(cmax)
      call loop_prop(cmax,i_cl(cmax),n_str_cl(cmax),i_str_cl)

c     output properties of unaliged regions
      call loop_prop_report(Lch_n)

c     compare local and global rmsd values for the template-free
c     fragment (TFF) regions. We choose to compare conformations to
c     cluster centroid conformation.
      call local_global_rms_correlation(cmax,i_cl(cmax),n_str_cl(cmax)
     &     ,i_str_cl)

cccccccccccccc output the nc clusters in PDB format ccccccccccccccc
      do i=1,nc
         if(i.lt.10)then
           c2=char(48+i)//'.pdb'
         else
           c2=char(48+i/10)//char(48+i-10*(i/10))//'.pdb'
         endif
***   'combo.pdb'->average of structures of cluster-------->
         open(10,file='combo'//c2,status='unknown')
         do j=1,Lch
            xx=xc(i,j)
            yy=yc(i,j)
            zz=zc(i,j)
            write(10,1237)'ATOM',j,'CA',seq(j),j,'',xx,yy,zz
          enddo
          write(10,322)
322       format('TER')
         do j=2,Lch
            write(10,1238)'CONECT',j-1,j
         enddo
         close(10)
***   'combc.pdb'->combined structures of cluster-------->
         open(10,file='combc'//c2,status='unknown')
         do j=1,Lch
            xx=xc3(i,j)
            yy=yc3(i,j)
            zz=zc3(i,j)
            write(10,1237)'ATOM',j,'CA',seq(j),j,'',xx,yy,zz
         enddo
          write(10,322)
         do j=2,Lch
            write(10,1238)'CONECT',j-1,j
         enddo
         close(10)
***   'model.pdb'->center structure of cluster-------->
         open(10,file='model'//c2,status='unknown')
         k=i_cl(i)
         do j=1,Lch
            xx=x(j,k)
            yy=y(j,k)
            zz=z(j,k)
            write(10,1237)'ATOM',j,'CA',seq(j),j,'',xx,yy,zz
         enddo
          write(10,322)
         do j=2,Lch
            write(10,1238)'CONECT',j-1,j
         enddo
         close(10)
***   'closc.pdb'->structure closest to combo.pdb-------->
         open(10,file='closc'//c2,status='unknown')
         do j=1,Lch
            xx=xcl(i,j)
            yy=ycl(i,j)
            zz=zcl(i,j)
            write(10,1237)'ATOM',j,'CA',seq(j),j,'',xx,yy,zz
         enddo
          write(10,322)
         do j=2,Lch
            write(10,1238)'CONECT',j-1,j
         enddo
         close(10)
       enddo
 1237 format(A4,I7,A4,A5,I6,A4,3F8.3)
 1238 format(A6,i5,i5)
c^^^^^^^^^^^^^^ output PDB structures finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      write(20,*)'Target: ',protein
      write(20,*)'Modeling Length: ',Lch
      write(20,*)'Native Length: ',Lch_n
      write(20,*)'Number of clusters: ',nc
      write(20,*)'Total number of structures=',n_str_all
      write(20,*)'Number of structure in use=',n_str
      write(20,*)
cccccccccccccc RMSD of cluster to native cccccccccccccccccccccccccc
      if(Lch_n.gt.1)then
         write(20,*)'--------- comparison to native structure---------'
         write(20,*)'  i R_combo   R_combc  R_model R_closc  R_piece '
         write(20,*)'A-------------------------------------------------'
         do i=1,nc
            do j=1,Lch_n
               r_1(1,j)=x_n(j)  !CA
               r_1(2,j)=y_n(j)
               r_1(3,j)=z_n(j)
            enddo
***   combo->
            do j=1,Lch_n
               r_2(1,j)=xc(i,ires(j))
               r_2(2,j)=yc(i,ires(j))
               r_2(3,j)=zc(i,ires(j))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n) !RMSD12
            rmsd_combo=armsd
***   combc->
            do j=1,Lch_n
               r_2(1,j)=xc3(i,ires(j))
               r_2(2,j)=yc3(i,ires(j))
               r_2(3,j)=zc3(i,ires(j))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n) !RMSD12
            rmsd_combc=armsd
***   model->
            k=i_cl(i)
            do j=1,Lch_n
               r_2(1,j)=x(ires(j),k)
               r_2(2,j)=y(ires(j),k)
               r_2(3,j)=z(ires(j),k)
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n) !RMSD12
            rmsd_model=armsd
***   closc->
            do j=1,Lch_n
               r_2(1,j)=xcl(i,ires(j))
               r_2(2,j)=ycl(i,ires(j))
               r_2(3,j)=zcl(i,ires(j))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n) !RMSD12
            rmsd_closc=armsd
***   piece of combo->
            nrm=0
            do j=1,Lch_n
               if(ires(j).ge.n1.and.ires(j).le.n2)then
                  nrm=nrm+1
                  r_1(1,nrm)=x_n(j)
                  r_1(2,nrm)=y_n(j)
                  r_1(3,nrm)=z_n(j)
               endif
            enddo
            nrm=0
            do j=1,Lch_n
               if(ires(j).ge.n1.and.ires(j).le.n2)then
                  nrm=nrm+1
                  r_2(1,nrm)=xc(i,ires(j))
                  r_2(2,nrm)=yc(i,ires(j))
                  r_2(3,nrm)=zc(i,ires(j))
               endif
            enddo
            call u3b(w,r_1,r_2,nrm,0,rms,u,t,ier)
            armsd=dsqrt(rms/nrm) !RMSD12
            rmsd_piece=armsd
***   
            write(20,50)i,rmsd_combo,rmsd_combc,rmsd_model,
     $           rmsd_closc,rmsd_piece
 50         format(i5,10f8.3)
         enddo
         write(20,*)
         txt1='                         '
         txt2='RMSD   Energy   #str    Trajectory'
         write(20,*)txt1(1:23)//txt2(1:35)
         write(20,*)'-------------------------------------------------'
         i=itra_rmsd_min_all
         j=itra_rmsd_min_use
         k=itra_E_min
         l=itra_E_min_all
         write(20,124)rmsd_min_all,E_rma,istr_rmsd_min_all,filen(i)
         write(20,125)rmsd_min_use,E_rmu,istr_rmsd_min_use,filen(j)
         write(20,126)rmsd_E_min,E_min,istr_E_min,filen(k)
         write(20,127)rmsd_E_min_all,E_min_all,istr_E_min_all,filen(l)
 124     format('Minimum RMSD in all=',f8.3,f9.1,i8,2x,a20)
 125     format('Minimum RMSD in use=',f8.3,f9.1,i8,2x,a20)
 126     format('Minimum E    in all=',f8.3,f9.1,i8,2x,a20)
 127     format('Minimum E    in use=',f8.3,f9.1,i8,2x,a20)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

****************************************************************
*     structure closest to template
****************************************************************
      if(L_temp.gt.1)then
***   output:
         open(10,file='templ.pdb',status='unknown')
         write(10,*)'templ',Lch
         do j=1,Lch
            write(10,1237)'ATOM',j,'CA',seq(j),j,''
     $           ,xtemp(j),ytemp(j),ztemp(j)
         enddo
         do j=2,Lch
            write(10,1238)'CONECT',j-1,j
         enddo
         close(10)
         write(20,*)
         write(20,*)'-------- structure closest to template ------'
         write(20,103)i_str_temp,filen(i_tra_temp)
 103     format('The structure in trajectory:',i8,2x,a20)
         write(20,*)'RMSD/aligned-length to templte:',R_temp_min,L_temp
***   comparison to 'CA':
         if(Lch_n.gt.1)then
            do j=1,Lch_n
               r_1(1,j)=x_n(j)
               r_1(2,j)=y_n(j)
               r_1(3,j)=z_n(j)
               r_2(1,j)=xtemp(ires(j))
               r_2(2,j)=ytemp(ires(j))
               r_2(3,j)=ztemp(ires(j))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n) !RMSD12
            rmsd_templ=armsd
            write(20,*)'RMSD to native:',rmsd_templ
         endif
      endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      write(20,*)

cccccccccccccccccc output cluster analysis cccccccccccccccccccc
      write(20,*)'------------ summary of clusers -----------'
      txt1='i Size R_cut density R_cl_co '
      txt2='<E>    E_model  #str Trajectory'
      write(20,*)txt1(1:31)//txt2(1:35)
      write(20,*)'B--------------------------------------------------'
      do i=1,nc
         k=i_cl(i)
         density=float(n_str_cl(i))/rmsd_cl_cut(i)
         write(20,123)i,n_str_cl(i),rmsd_cl_cut(i),density,
     $        rmsd_close_min(i),E_combo(i),E(k),istr(k),filen(itra(k))
      enddo
 123  format(i2,i6,f6.2,f7.0,f6.2,2f9.1,i6,1x,a13)

cccccccccccccccccc output cluster analysis cccccccccccccccccccc
      write(20,*)
      write(20,*)' include used structure  exclude used structre'
      write(20,*)'  ---------------------  ---------------------'
      write(20,*)'i  N_in  <R_in> <Rc_in>   N_ex  <R_ex> <Rc_ex>'
      write(20,*)'C---------------------------------------------'
      do i=1,nc
         write(20,133)i,n_str_cl(i),R_in(i),Rc_in(i),
     $        n_str_cl_ex(i),R_ex(i),Rc_ex(i)
      enddo
 133  format(i2,i6,f7.2,f7.2,i8,f7.2,f7.2)

      write(20,*)
      write(20,*)'RMSD_cut_initial=',RMSD_cut_initial
      write(20,*)'RMSD_cut_min=',RMSD_cut_min
      write(20,*)'RMSD_cut_max=',RMSD_cut_max
      write(20,*)'ratio_min=',ratio1
      write(20,*)'ratio_max=',ratio2
      write(20,*)
      write(20,*)'trajectories used:',n_tra
      do i=1,n_tra
         write(20,*)filen(i)
      enddo
c^^^^^^^^^^^^^^^^^^^^ output done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
c     do spicker on every loop. Do spicker only on structures belonging
c     to most dense cluster 
      call select_structures(.true.) !select only str. of most dense cluster
      do iloop=1,nloop
         if(lextent(iloop).lt.minloopL) goto 12 !loop too short         
         call clusterize(iloop,amat,.true.) !do clustering of loop local rmsd values
         call make_combo(iloop,.true.) !do combo structures for the loop
         call report_clustering(iloop,.true.) !do rst_iloop.dat-like file
         call clusterize(iloop,amat,.false.) !do clustering of loop global rmsd values
         call make_combo(iloop,.false.) !do combo structures for the loop
         call report_clustering(iloop,.false.) !do rst_iloop.dat-like file        
 12   enddo

      call exit(0) !exit program cluster with success
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
 
c     ###########################################

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
      subroutine u4b(w, x, y, start, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m,mm, ier, n, mode,start
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
      mm=start+m
      if (w(m) .lt. 0.0) return 
      wc = wc + w(mm)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(mm) * x(i,mm))
    2 yc(i) = yc(i) + (w(mm) * y(i,mm))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      mm=start+m
      do 4 i = 1, 3
      e0 = e0 + (w(mm) * (((x(i,mm)-xc(i)) ** 2) + ((y(i,mm)-yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(mm) * (y(i,mm) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,mm) - xc(j)))
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
c.....END U4B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end
          
c     ###########################################

c     read exp.dat and store in exp variable
      subroutine readexpdat(Lch)
      parameter(ndim=1000)
      parameter(nst=13000)      !number of used structure, maximum allowed
      common/expdat/expd,hsc
      real expd(ndim)
      integer hsc(ndim)
      integer Lch
      integer i,j,nc,nr,mp(20)
      real maxexp

      
      maxexp=0.85 !maximum predicted surface area
      nc=17 !number of columns in exp.dat devoted to exposed area

      read(8,*) !first line is exp.dat a comment
      do i=1,Lch
         read(8,*)nr,(mp(j),j=1,nc)
         expd(i)=maxexp !initialization to totally exposed
         do j=1,nc
            if(mp(j).eq.0) then
               expd(i)=(maxexp*j)/20
               goto 10
            endif
         enddo
 10   enddo

      end

c     ###########################################

c     read chain.dat file, initializing isalig,lbegin,lextent,nloop,lextents
      subroutine readchaindat(Lch,itemp)
      parameter(ndim=1000)
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(stem=5) !how much to extend either loop
      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      common/end_to_end/end_to_end_dist
      integer Lch
      real end_to_end_dist(ndim)
      integer isalig(0:ndim),lbegin(ndim),lextent(ndim),nloop,
     $     ntemp,itemp,jtemp,ncov,i,j,k,b,e,lextents(ndim),
     $     lresid(ndim,ndim)
      integer nterm_stem,cterm_stem
      real r,x0,y0,z0,trans
      dimension r(3,ndim) !template coordinates
      character buffer*1
c     initialize every residue as not covered by the template
      do i=1,Lch
         isalig(i)=0
      enddo

      read(7,*) ntemp !number of templates

c     find covered residues and their coordinates. We assume a standard
c     format for chain.dat
      jtemp=1
 6    if(jtemp.lt.itemp)then    !jump to next template
         read(7,*) ncov         !number of residues covered by the template
         do i=1,ncov
            read(7,*) buffer 
         enddo
         jtemp=jtemp+1
         goto 6 !iterate until finding itemp template
      endif
      read(7,*) ncov            !number of residues covered by the template
      do i=1,ncov
         read(7,*) j,x0,y0,z0
         isalig(j)=1            !residue at position j is aligned
         r(1,j)=x0              !store template coordinates
         r(2,j)=y0
         r(3,j)=z0
      enddo

c     find number of unaligned regions, termed "loops", as well as where
c     they begin and their length
      nloop=0
      isalig(0)=1 !boundary condition to deal with unaligned N-terminal
      isalig(Lch+1)=1 !boundary condition to deal with unaligned C-terminal
      do i=1,Lch+1
         if(isalig(i).eq.0 .and. isalig(i-1).eq.1)then
            nloop=nloop+1
            lbegin(nloop)=i !position where current loop begins
         elseif(isalig(i).eq.1 .and. isalig(i-1).eq.0)then
            lextent(nloop)=i-lbegin(nloop) !extent of current loop
         endif
      enddo
c     find extended loop for every loop, that is, including the stems
      do i=1,nloop
c     first initialize the extended loop as current loop
         do j=1,lextent(i)
            lresid(j,i)=lbegin(i)+j-1
         enddo
         lextents(i)=lextent(i)
c     next add the stems on both sides, totalling another lextents(i)
c     residues, at most. Stems only included aligned residues
         nterm_stem=0
         cterm_stem=0
         b=lbegin(i)-1 !initialize outer left boundary of the loop
         e=lbegin(i)+lextent(i) !initialize outer right boundary of the loop
 10      if(b.gt.0)then !extend on the left
            if(isalig(b).eq.1) then
               if(nterm_stem.eq.0) nterm_stem=b !first stem of N-terminal
               lextents(i)=lextents(i)+1
               lresid(lextents(i),i)=b
            endif
            b=b-1
         endif
         if(e.le.Lch)then !extend on the right
            if(isalig(e).eq.1) then
               if(cterm_stem.eq.0) cterm_stem=e !first stem of C-terminal
               lextents(i)=lextents(i)+1
               lresid(lextents(i),i)=e               
            endif
            e=e+1
         endif
         if(lextents(i).lt.2*lextent(i) .and. 
     $        (b.gt.0 .or. e.le.Lch) ) goto 10
c     calculate loop end_to_end distance. Since distance between
c     template residues are quite conserved, we use the distance between
c     first stem residues of N-terminal and C-terminal of loop ends
         end_to_end_dist(i)=-1.0 !marks it's not a loop but an unaligned terminus
         if(nterm_stem.ne.0 .and. cterm_stem.ne.0)then
            do j=1,3
               trans=r(j,nterm_stem)-r(j,cterm_stem)
               end_to_end_dist(i)=end_to_end_dist(i)+trans*trans
            enddo
            end_to_end_dist(i)=sqrt( end_to_end_dist(i) )
         endif
      enddo

      end
                  
c     ###########################################

c     Derived properties for each of the template-free fragment
c     (TFF) regions
c     ic: ID of cluster
c     im: ID of the center cluster structure for cluster i
c     nstr number structures belonging to cluster ic
c     i_str_cl(i,j) ID of one of the structures of cluster i
c     x(k,m): x-coordinate of residue k in structure m
c     lrmsdl: local rmsd of the region
c     lrmsds: global rmsd of the region, ie, including stems. We do an
c     rmsd superposition of the TFF plus stems, but only report the
c     rmsd of the TFF, not the rmsd of the loop plus stems. The reason
c     is that stems have very low rmsd because they are part of the
c     template, thus an rmsd value for both TFF and stems will be low
c     and give the impression the whole region is great. In addition,
c     the reported rmsd value is comparable to the local rmsd value,
c     because it only involves TFF residues.
c
c     lrmsdns: new global rmsd of the region. Given a TFF of size "L",
c     we extend the region with stems of size L/2 on each terminus. Then
c     we do pairwise rmds between cluster member and cluster center but
c     only taking into the rmsd calculation the two L/2 stems, nor the
c     TFF region itself, neither the rest of the protein. After aligning
c     the stems of two snapshots, we look at the conformation of the TFF
c     on each snapshot. In particular, we first superimpose the center
c     of mass of the two conformations by translation and then calculate
c     1/2*sum_{i=1}^{L}( {vec(r_i^{I}) - vec(r_i^{II})|^2) , which is
c     the rmsd function but without doing any rotation to minimize its
c     value. Note this new global rmsd value is calculated only for the
c     TFF fragment, thus it scales with "L" same as the local rmsd.
c
c     lrmsdns_ha, lrmsdl_ha: same as lrmsdns and lrmsdl, but only when
c     rmsd between stems is below a restrictive cut-off
c     =1+(6.5-1)*(L-4)/(100-4), where L TFF length. The point of this
c     measurement is to observe how local and global rmsd correlate when
c     the stems are relatively fixed, to see if a global change in rmsd
c     must neccessarily destroy the local conformation of the loop, thus
c     increasing the local rmsd

c     lextentg: length of the unaligned regions plus stems
c     lexp: average predicted exposure of the region
c     lhelix:fraction of residues predicted in the helical state
c     lstrand:fraction of residues predicted in the extended state
c     lcoil=1-lhelix-lstrand
      subroutine loop_prop(ic,im,nstr,i_str_cl)
      parameter(ndim=1000)      !maximum sequence length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(ncl=100)        !number of clusters
      parameter(minloopL=4)     !minimum loop length to consider

      common/native/Lch_n,x_n,y_n,z_n
      integer Lch_n
      dimension x_n(ndim),y_n(ndim),z_n(ndim) !native structure

      common/expdat/expd,hsc
      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      common/proploop/lrmsdl,lrmsds,lexp,lhelix,lcoil,lstrand,
     &nrmsdl,nrmsds,lrmsdns,nrmsdns
      common/whole/Lch,x,y,z
      integer i,ii,j,k,l,b,e,isalig(0:ndim),lbegin(ndim),lextent(ndim),
     &nloop,Lch,ic,im,i_str_cl(ncl,nst),ier,
     &hsc(ndim),lextents(ndim),lresid(ndim,ndim)

      real lexp(ndim),armsd,lrmsdl(ndim),lrmsds(ndim),expd(ndim),
     &lrmsdns(ndim),lhelix(ndim),lstrand(ndim),lcoil(ndim),
     &nrmsdl(ndim),nrmsds(ndim),rmsdns,residual,residual2,
     &nrmsdns(ndim)

      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000),
     &r_4(3,ndim)
      double precision u(3,3),t(3),rms,drms
      data w /1000*1.0/
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures

c      write(*,*) 'nloop=',nloop
      do j=1,nloop         
         Lstem=lextents(j)-lextent(j) 
         if(lextent(j).lt.minloopL) goto 10 !too short a loop. Don't process
         lrmsdl(j)=0.0  !initialize average local rmsd
         lrmsds(j)=0.0  !initialize average stem rmsd
         lrmsdns(j)=0.0 !initialize average new stem rmsd
         lexp(j)=0.0   !initialize average predicted solvent exposure
         lhelix(j)=0.0 !initialize helix content
         lstrand(j)=0.0 !initialize extended content
         do i=lbegin(j),lbegin(j)+lextent(j)-1
            lexp(j)=lexp(j)+expd(i) !expd(i):predicted exposed surface of residue "i"
            if(hsc(i).eq.2) then !hsc: predicted secondary structure of residue "i"
               lhelix(j)=lhelix(j)+1.
            else if(hsc(i).eq.4) then
               lstrand(j)=lstrand(j)+1.
            endif
         enddo
         lexp(j)=lexp(j)/lextent(j)    !normalize by loop length
         lhelix(j)=lhelix(j)/lextent(j)
         lstrand(j)=lstrand(j)/lextent(j)
         lcoil(j)=1-lhelix(j)-lstrand(j)

         do i=1,nstr            !cycle over all structures belonging to cluster         
            k=i_str_cl(ic,i)    !ID of one of the structures belonging to cluster ic
            do b=1,lextents(j)  !cycle over all residues belonging to loop "j"
               l=lresid(b,j)   !residue number for one of the residues
               r_1(1,b)=x(l,im)!coordinates from the cluster center
               r_1(2,b)=y(l,im)
               r_1(3,b)=z(l,im)
               r_2(1,b)=x(l,k) !coordinates from the particular structure
               r_2(2,b)=y(l,k)
               r_2(3,b)=z(l,k)
            enddo
c     rmsd for loop between center cluster "im" and structure "k". We
c     take advantage that the first lextent(j) residues in r_1 and r_2
c     do correspond to the loop withouth the stems
            call u3b(w,r_1,r_2,lextent(j),0,rms,u,t,ier)
            lrmsdl(j)=lrmsdl(j)+dsqrt(rms/lextent(j))
c     rmsd for loop including stems between center cluster "im" and structure "k"
            call u3b(w,r_1,r_2,lextents(j),1,rms,u,t,ier)
            lrmsds(j)=lrmsds(j)+residual2(r_1,r_2,lextent(j),u,t)

         enddo         
         lrmsdl(j) =lrmsdl(j)/nstr !normalize by the number of structures belonging to cluster
         lrmsds(j) =lrmsds(j)/nstr

         k=1
         do k2=lextent(j)+1,lextents(j)     !cycle over stem residues belonging to loop "j"
            l=lresid(k2,j)       !residue number for one of the residues
            r_1(1,k)=x(l,im)    !coordinates from the cluster center
            r_1(2,k)=y(l,im)
            r_1(3,k)=z(l,im)
            k=k+1
         enddo
         do k=1,lextent(j) !cycle over local residues belonging to loop "j"
            l=lresid(k,j)      !residue number for one of the residues
            r_3(1,k)=x(l,im)    !coordinates from the cluster center
            r_3(2,k)=y(l,im)
            r_3(3,k)=z(l,im)           
         enddo
         do ii=1,nstr !cycle over all structures belonging to cluster ic
            i=i_str_cl(ic,ii)   !ID of one of the structures belonging to cluster ic
c            write(*,*)'i=',i
            k=1
            do k2=lextent(j)+1,lextents(j) !cycle over stem residues belonging to loop "j"
               l=lresid(k2,j)    !residue number for one of the residues
               r_2(1,k)=x(l,i)  !coordinates from the particular structure
               r_2(2,k)=y(l,i)
               r_2(3,k)=z(l,i)               
               k=k+1
            enddo
            call u3b(w,r_1,r_2,Lstem,1,rms,u,t,ier) !align stems-->rot. matrix u
            do k=1,lextent(j)   !cycle over local residues belonging to loop "j"
               l=lresid(k,j)    !residue number for one of the residues
               r_4(1,k)=x(l,i) !coordinates from one of the structures
               r_4(2,k)=y(l,i)
               r_4(3,k)=z(l,i)           
            enddo
            rmsdns=residual(r_3,r_4,lextent(j),u,t)
c            write(*,'f5.2') rmsdns
            lrmsdns(j)=lrmsdns(j)+rmsdns !lrmsdns provided previous rot. u
         enddo
         lrmsdns(j)=lrmsdns(j)/nstr

c         write(*,*)'j=',j,' Lstem=',Lstem,' lrmsdl=',lrmsdl(j),
c     &        ' lrmsds=',lrmsds(j),' lrmsdns=',lrmsdns(j)

         if(Lch_n.gt.1)then     !comparison to native
            do b=1,lextents(j)
               l=lresid(b,j)    !residue number for one of the residues
               r_1(1,b)=x_n(l)  !coordinates from the native state
               r_1(2,b)=y_n(l)
               r_1(3,b)=z_n(l)
               r_2(1,b)=x(l,im) !coordinates of closc for cluster ic
               r_2(2,b)=y(l,im)
               r_2(3,b)=z(l,im)
            enddo
            call u3b(w,r_1,r_2,lextent(j),0,rms,u,t,ier)
            nrmsdl(j)=dsqrt(rms/lextent(j))
            call u3b(w,r_1,r_2,lextents(j),1,rms,u,t,ier)
            nrmsds(j)=residual2(r_1,r_2,lextent(j),u,t)            
            k=1
            do k2=lextent(j)+1,lextents(j) !cycle over stem residues belonging to loop "j"
               l=lresid(k2,j)   !residue number for one of the residues
               r_1(1,k)=x_n(l)  !native coordinates
               r_1(2,k)=y_n(l)
               r_1(3,k)=z_n(l)
               r_2(1,k)=x(l,im) !coordinates from the cluster center
               r_2(2,k)=y(l,im)
               r_2(3,k)=z(l,im)
c               write(*,'3(f8.3,x)')r_2(1,k),r_2(2,k),r_2(3,k)
               k=k+1
            enddo
            call u3b(w,r_1,r_2,Lstem,1,rms,u,t,ier)
c            write(*,*)rms
            do k=1,lextents(j) !cycle over stem residues belonging to loop "j"
               l=lresid(k,j)   !residue number for one of the residues
               r_1(1,k)=x_n(l)  !native coordinates
               r_1(2,k)=y_n(l)
               r_1(3,k)=z_n(l)
               r_2(1,k)=x(l,im) !coordinates from the cluster center
               r_2(2,k)=y(l,im)
               r_2(3,k)=z(l,im)
            enddo
            nrmsdns(j)=residual(r_1,r_2,lextent(j),u,t)
c            write(*,'i4,x,i2,x,f5.2')im,j,nrmsdns(j)
         endif

 10   enddo

      end
                  
c     ###########################################

      subroutine loop_prop_report(Lch_n)
      parameter(ndim=1000)      !maximum sequence length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(minloopL=4)     !minimum loop length to consider
      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      common/proploop/lrmsdl,lrmsds,lexp,lhelix,lcoil,lstrand,
     &nrmsdl,nrmsds,lrmsdns,nrmsdns

      integer i,j,k ,isalig(0:ndim),nloop,lbegin(ndim),
     &lextent(ndim),lextents(ndim),lresid(ndim,ndim),Lch_n,
     &Lstem

      real lrmsdl(ndim),lrmsds(ndim),lexp(ndim),lhelix(ndim),
     &lcoil(ndim),lstrand(ndim),nrmsdl(ndim),nrmsds(ndim),
     &lrmsdns(ndim),nrmsdns(ndim)

      common/whole/Lch,x,y,z
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures
      integer Lch

      open(1,file='loop.dat',status='unknown') !This unit is local to the subroutine

      write(1,'A')'#properties for unaligned regions'
      write(1,'A')'#lid: loop id'
      write(1,'A')'#l_st: loop plus stems length'
      write(1,'A')'#<rmsd>: average local rmsd'
      write(1,'A')'#<r_gl1>:average global rmsd. Align both loop and'
      write(1,'A')'#        stems, then report only rmsd of loop region'
      write(1,'A')'#<r_gl2>:average global rmsd. Align only stems and'//
     $     'obtain rotation'
      write(1,'A')'#        matrix which we use to calculate rmsd of '//
     $     'loop region'
      write(1,'A')'#lid begin end l <exp>   h     e     c  <rmsd> l_st'//
     &     ' <r_gl1> <r_gl2>'

 10   format(i2,x,3(x,i3),4(x,f5.3),x,f5.2,x,i3,2x,f5.2,3x,f5.2)
 11   format(i2,x,3(x,i3))
      do i=1,nloop
         if(lextent(i).lt.minloopL)then !too short a loop. Don't process
            write(1,11)i,lbegin(i),lbegin(i)+lextent(i)-1,lextent(i)
         else
            write(1,10)i,lbegin(i),lbegin(i)+lextent(i)-1,lextent(i),
     &           lexp(i),lhelix(i),lstrand(i),lcoil(i),lrmsdl(i),
     &           lextents(i),lrmsds(i),lrmsdns(i)            
         endif
 15   enddo

c     report comparison to native
 20   format(i2,x,3(x,i3),3(x,f5.2))
 21   format(i2,x,3(x,i3))
      if(Lch_n.gt.1)then
         write(1,'A')'#comparison to native'
         write(1,'A')'#lid: loop id'
         write(1,'A')'#rmsd: local rmsd to native'
         write(1,'A')'#r_gl1: global rmsd to native. Same method as'//
     $        'in <r_gl 1>'
         write(1,'A')'#r_gl2: global rmsd to native. Same method as'//
     $        'in <r_gl2>'
         write(1,'A')'#lid begin end l  rmsd r_gl1 r_gl2'         
         do i=1,nloop
            if(lextent(i).lt.minloopL) then!too short a loop. Don't process
               write(1,21)i,lbegin(i),lbegin(i)+lextent(i)-1,lextent(i)
            else
               write(1,20)i,lbegin(i),lbegin(i)+lextent(i)-1,lextent(i),
     &              nrmsdl(i),nrmsds(i),nrmsdns(i)
            endif
 25      enddo
      endif

      end

c     ###########################################

c     compare local and global rmsd values for the template-free
c     fragment (TFF) regions. The aim is to observe if TASSER has to
c     melt the local conformation of the TFF in order to globally move
c     the TFF. If we see high global-rmsd values concomitant to low
c     local-rmsd values, then we now that TASSER is able to move the TFF
c     as a whole. We choose to compare conformations to cluster centroid
c     conformation.

      subroutine local_global_rms_correlation(ic,im,nstr,i_str_cl)

      parameter(ndim=1000)      !maximum sequence length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(ncl=100)        !number of clusters
      parameter(minloopL=4)     !minimum loop length to consider
      common/whole/Lch,x,y,z
      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid

      integer isalig(0:ndim),nloop,lbegin(ndim),lextent(ndim),
     $     lextents(ndim),lresid(ndim,ndim)

      integer i,j,k,k2,l,ic,im,i_str_cl(ncl,nst),ier,Lstem
      integer Lch
      real stem_rmsd,lrmsdl,rmsdns

      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim),
     &r_4(3,ndim)
      double precision u(3,3),t(3),rms,drms
      data w /1000*1.0/
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures

      character buf*80

 10   format(3(x,f5.2))

      do j=1,nloop              !cycle over all unaligned regions
         Lstem=lextents(j)-lextent(j)
         if(lextent(j).lt.minloopL) goto 20 !too short a loop. Don't process
         if(j.lt.10)then !open file
            write(unit=buf,'ai1')'0',j
         else
            write(unit=buf,'i2')j
         endif
         open(1,file='loop_lg_'//buf(1:2)//'.dat',status='unknown')       
         write(1,'A')'#ref. structure is center of first cluster' !write header
         write(1,'A')'#global-rmsd does not include stems'
         write(1,'Ai3')'#Lstem=',Lstem
         write(1,'A')'#stem-rmsd local-rmsd global-rmsd' 
         k=1
         do k2=lextent(j)+1,lextents(j)     !cycle over stem residues belonging to loop "j"
            l=lresid(k2,j)       !residue number for one of the residues
            r_1(1,k)=x(l,im)    !coordinates from the cluster center
            r_1(2,k)=y(l,im)
            r_1(3,k)=z(l,im)
            k=k+1
         enddo
         do k=1,lextent(j) !cycle over local residues belonging to loop "j"
            l=lresid(k,j)      !residue number for one of the residues
            r_3(1,k)=x(l,im)    !coordinates from one of the structures
            r_3(2,k)=y(l,im)
            r_3(3,k)=z(l,im)           
         enddo
         do ii=1,nstr !cycle over all structures belonging to cluster ic
            i=i_str_cl(ic,ii) !ID of one of the structures belonging to cluster ic
            k=1
            do k2=lextent(j)+1,lextents(j) !cycle over stem residues belonging to loop "j"
               l=lresid(k2,j)    !residue number for one of the residues
               r_2(1,k)=x(l,i)   !coordinates from the particular structure
               r_2(2,k)=y(l,i)
               r_2(3,k)=z(l,i)               
               k=k+1
            enddo
            call u3b(w,r_1,r_2,Lstem,1,rms,u,t,ier) !align stems-->rot. matrix u
            stem_rmsd=dsqrt(rms/Lstem)
            do k=1,lextent(j)   !cycle over local residues belonging to loop "j"
               l=lresid(k,j)    !residue number for one of the residues
               r_4(1,k)=x(l,i) !coordinates from one of the structures
               r_4(2,k)=y(l,i)
               r_4(3,k)=z(l,i)           
            enddo
            rmsdns=residual(r_3,r_4,lextent(j),u,t) !global rmsd of the TFF
            call u3b(w,r_3,r_4,lextent(j),0,rms,u,t,ier)
            lrmsdl=dsqrt(rms/lextent(j)) !local rmsd of the TFF
            write(1,10)stem_rmsd,lrmsdl,rmsdns
         enddo
         close(1)
 20   enddo

      end


c     ###########################################

c     calculate residual \sqrt(1/2 Sum_i (U\vec(r_i^I)-vec(r_i^II))^2),
c     where U is rotation matrix that we pass along, intead of being
c     calculated to minimize the residuas as it is done in the typical
c     rmsd calculation
c     overwrite "t" with new translation

      function residual(x,y,n,u,t)
      parameter(ndim=1000)      !maximum sequence length     
      double precision x(3,1),y(3,1),u(3,3),xx(3,ndim),xcm(3),ycm(3),
     $     t(3)
      real residual
      integer i,n
c     rotate x according to u
      do j=1,n
         do i=1,3
            xx(i,j)=u(i,1)*x(1,j)+u(i,2)*x(2,j)+u(i,3)*x(3,j)
         enddo
      enddo
c      superimpose center of mass
      do i=1,3
         xcm(i)=0.0
         ycm(i)=0.0
         do j=1,n
            xcm(i)=xcm(i)+xx(i,j)
            ycm(i)=ycm(i)+y(i,j)
        enddo
         t(i)=( ycm(i)-xcm(i) )/n
      enddo

      do j=1,n
         do i=1,3
            xx(i,j)=t(i)+xx(i,j)
         enddo
      enddo
c      calculate residual
      residual=0.0
      do j=1,n
         do i=1,3
            residual=residual+(xx(i,j)-y(i,j))*(xx(i,j)-y(i,j))
         enddo
      enddo
      residual=sqrt(residual/n)
      end

c     ###########################################

      function residual2(x,y,n,u,t)
      parameter(ndim=1000)      !maximum sequence length     
      double precision x(3,1),y(3,1),u(3,3),xx(3,ndim),xcm(3),ycm(3)
     &,t(3)
      real residual2
      integer i,n
c     pitch x according to u,t
      do j=1,n
         do i=1,3
            xx(i,j)=t(i)+u(i,1)*x(1,j)+u(i,2)*x(2,j)+u(i,3)*x(3,j)
         enddo
      enddo
c      calculate residual
      residual2=0.0
      do j=1,n
         do i=1,3
            residual2=residual2+(xx(i,j)-y(i,j))*(xx(i,j)-y(i,j))
         enddo
      enddo
      residual2=sqrt(residual2/n)
      end


c     ###########################################
c     calculate RMSD matrics for a particular loop
c     input parameters are:
c     amat: matrix to store rmsd values
c     iloop: loop index identifying the loop we are to
c            cluster. iloop.le.0 means whole structure, not a loop
c     islocal: logical variable. True for local rmsd, False for global rmsd
      subroutine calc_rmsd_matrix(amat,iloop,islocal)
      parameter(ndim=1000)      !Length
      parameter(nst=13000)      !number of used structure, maximum allowed

      logical islocal
      integer iloop,b,e,ex,looplength,Lstem

      common/whole/Lch,x,y,z
      integer Lch
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures
      real amat
      dimension amat(nst,nst)   !RMSD matrics

      common/dens/n_str_all,n_str,d_rc_in,drcinmax,cmax,nstr_selected,
     $     str_selected
      real d_rc_in(100),drcinmax
      integer n_str_all,n_str,cmax,nstr_selected,str_selected(nst)

      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      integer isalig(0:ndim),lbegin(ndim),lextent(ndim),nloop,hsc(ndim),
     &lextents(ndim),lresid(ndim,ndim)

ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000),
     $     r_4(3,ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /1000*1.0/

c      write(*,'ai2') 'iloop=',iloop
c     init b,e,ex
      b=1 !begining sequence index
      e=Lch !ending sequence index
      ex=e-b+1 !extent
      if(iloop.gt.0)then
         looplength=lextent(iloop)
         Lstem=lextents(iloop)-looplength
         if(local.eq..true.)then
            b=lbegin(iloop)
            e=lbegin(iloop)+looplength-1
            ex=looplength
         endif
      endif
c      write(*,'4(xai4)') 'b=',b,'e=',e,'ex=',ex,
c     $     'nstr_selected=',nstr_selected

      do i=1,nstr_selected-1
         ii=str_selected(i)     !structure ID
         amat(ii,ii)=0.0
         do j=i+1,nstr_selected
            jj=str_selected(j)  !structure ID
            if(iloop.lt.1 .or. islocal.eq..true.)then !whole sequence or local loop rmsd
               ik=1
               do k=b,e
                  r_1(1,ik)=x(k,ii)
                  r_1(2,ik)=y(k,ii)
                  r_1(3,ik)=z(k,ii)
                  r_2(1,ik)=x(k,jj)
                  r_2(2,ik)=y(k,jj)
                  r_2(3,ik)=z(k,jj)
                  ik=ik+1
               enddo
               call u3b(w,r_1,r_2,ex,0,rms,u,t,ier)
               armsd=dsqrt(rms/ex) !if(iloop.gt.0) write(*,*)armsd
            else !global loop rmsd
               ik=1
               do k=looplength+1,lextents(iloop) !cycle over stem res. belonging to loop "iloop"
                  l=lresid(k,iloop) !residue number for one of the residues
                  r_1(1,ik)=x(l,ii) !coordinates from the cluster center
                  r_1(2,ik)=y(l,ii)
                  r_1(3,ik)=z(l,ii)
                  ik=ik+1
               enddo
               do k=1,looplength !cycle over local residues belonging to loop "iloop"
                  l=lresid(k,iloop) !residue number for one of the residues
                  r_3(1,k)=x(l,ii) !coordinates from the cluster center
                  r_3(2,k)=y(l,ii)
                  r_3(3,k)=z(l,ii)           
               enddo
               ik=1
               do k=looplength+1,lextents(iloop) !cycle over stem res. belonging to loop "iloop"
                  l=lresid(k,iloop) !residue number for one of the residues
                  r_2(1,k)=x(l,jj) !coordinates from the particular structure
                  r_2(2,k)=y(l,jj)
                  r_2(3,k)=z(l,jj)               
                  ik=ik+1
               enddo
               call u3b(w,r_1,r_2,Lstem,1,rms,u,t,ier) !align stems-->rot. matrix u
               do k=1,looplength !cycle over local residues belonging to loop "iloop"
                  l=lresid(k,iloop) !residue number for one of the residues
                  r_4(1,k)=x(l,jj) !coordinates from one of the structures
                  r_4(2,k)=y(l,jj)
                  r_4(3,k)=z(l,jj)           
               enddo
               armsd=residual(r_3,r_4,looplength,u,t) !calculate "rmsd" value
c               write(*,'f5.2')armsd
            endif
            amat(ii,jj)=armsd
            amat(jj,ii)=armsd
         enddo
      enddo
      amat(nstr_selected,nstr_selected)=0.0
      end

c     ###########################################
c     subroutine clusterize performs clustering of a particular
c     unaligned region or the whole sequence.
c     input parameters are:
c     iloop: loop index identifying the loop we are to
c            cluster. iloop.le.0 means whole structure, not a loop
c     amat: matrix to store rmsd values
      subroutine clusterize(iloop,amat,islocal)
      parameter(ndim=1000)      !Length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(ntr=20000)       !number of trajectories
      parameter(ncl=100)            !number of clusters
      parameter(minloopL=4)     !minimum loop length to consider
      logical islocal
      common/clusteriza/nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,
     $     i_str_cl_ex,rmsd_cl_cut
      integer nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,i_str_cl_ex
      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      common/end_to_end/end_to_end_dist
      integer isalig(0:ndim),lbegin(ndim),lextent(ndim),nloop,hsc(ndim),
     &lextents(ndim),lresid(ndim,ndim)
      integer mark,n_str_near
      real amat,end_to_end_dist,R_in,R_ex
      dimension end_to_end_dist(ndim)
      dimension mark(nst)       !for removing used structures
      dimension amat(nst,nst)    !RMSD matrics
      dimension n_str_near(nst) !number of neighboring structures
      dimension i_cl(ncl) !i_cl(i) is structure center of i'th cluster
      dimension n_str_cl(ncl)!of structures including some used struct in i'th cluster
      dimension n_str_cl_ex(ncl) !n_str_cl_ex(i): # of structures excluding some used struct
                                 ! in i'th cluster
      dimension i_str_cl(ncl,nst) !i_str_cl(i,j) is structure id of j'th structure belonging
                                  !to i'th cluster
      dimension i_str_cl_ex(ncl,nst) !i_str_cl_ex(i,j): structure id of j'th structure
                                     !belonging to i'th cluster
      real rmsd_cl_cut
      dimension rmsd_cl_cut(ncl) !rmsd_cl_cut(i)

      common/r_/R_in,R_ex,Rc_in,Rc_ex
      dimension R_in(100),R_ex(100),Rc_in(100),Rc_ex(100)

      common/dens/n_str_all,n_str,d_rc_in,drcinmax,cmax,nstr_selected,
     $     str_selected
      real d_rc_in(100),drcinmax
      integer n_str_all,n_str,cmax,nstr_selected,str_selected(nst)

      call calc_rmsd_matrix(amat,iloop,islocal) !calculate RMSD matrics

c     find out the biggest cluster and decide the RMSD_cut_min
      if(iloop.gt.0)then
         if(islocal.eq..true.)then !loop local rmsd
            rmsd_cut_ini=random_rmsd(lextent(iloop),
     $           end_to_end_dist(iloop),.true.) !expected random rmsd for loop local
         else !loop global rmsd
            rmsd_cut_ini=random_rmsd(lextent(iloop),
     $           end_to_end_dist(iloop),.false.) !expected random rmsd for loop global
         endif
         rmsd_cut_max=2.0*rmsd_cut_ini
         rmsd_cut_min=0.5*rmsd_cut_ini
         rmsd_cut=rmsd_cut_max  !start with generous cut-off
      else
         rmsd_cut_ini=8.0
         rmsd_cut_max=12.0
         rmsd_cut_min=3.5
         rmsd_cut=rmsd_cut_ini
      endif
      ratio1=0.7
      ratio2=0.15
      nc_max=10

 140  do j=1,nstr_selected
         jj=str_selected(j)
         mark(jj)=1
         n_str_near(jj)=0        !number of structure in cluster
         do k=1,nstr_selected
            kk=str_selected(k)
            if(amat(jj,kk).lt.rmsd_cut)then
               n_str_near(jj)=n_str_near(jj)+1
            endif
         enddo
c         if(iloop.gt.0)write(*,'ai4ai4')'n_str_near(',jj,')=',
c     $        n_str_near(jj)
      enddo

c     find out the biggest cluster
      n_str_cl_max=0
      do j=1,nstr_selected
         jj=str_selected(j)
         if(n_str_near(jj).gt.n_str_cl_max)then
            n_str_cl_max=n_str_near(jj)
         endif
      enddo
c     check the size of the cluster
      ratio=float(n_str_cl_max)/float(n_str_selected)
      if(ratio.gt.ratio1.and.rmsd_cut.gt.rmsd_cut_min)then
         rmsd_cut=rmsd_cut-0.1 !too big a cluster, thus restrict rmsd cutoff
         goto 140
      endif
      if(ratio.lt.ratio2.and.rmsd_cut.lt.rmsd_cut_max)then
         rmsd_cut=rmsd_cut+0.2 !too small a cluster, thus increase rmsd cutoff
         goto 140
      endif
      if(rmsd_cut.lt.rmsd_cut_min)rmsd_cut=rmsd_cut_min
      if(rmsd_cut.gt.rmsd_cut_max)rmsd_cut=rmsd_cut_max
c     find out nc_max clusters
      nc=0                      !number of clusters
      do i=1,nc_max
c     calculate nc(j), number of neightboring structures ----------------->
         do j=1,nstr_selected
            jj=str_selected(j)
            n_str_near(jj)=0     !number of structure in cluster
            do k=1,nstr_selected
               kk=str_selected(k)
               if(mark(jj).eq.1.and.mark(kk).eq.1)then
                  if(amat(jj,kk).lt.rmsd_cut)then
                     n_str_near(jj)=n_str_near(jj)+1
                  endif
               endif
            enddo
         enddo
c     find out the biggest cluster
         n_str_cl_max=0         !maximum number of structures in cluster
         do j=1,nstr_selected
            jj=str_selected(j)
            if(n_str_near(jj).gt.n_str_cl_max)then
               n_str_cl_max=n_str_near(jj)
               i_cl(i)=jj       !structure center of i'th cluster
            endif
         enddo
c         write(*,'ai4')'n_str_cl_max=',n_str_cl_max

c     check the size of the cluster 
         if(n_str_cl_max.lt.1) goto 41 !no remaining clusters.
c     remove the structures in the cluster
         n_str_cl(i)=0          !# of structure including some used struct
         n_str_cl_ex(i)=0       !# of structure excluding used structure
         R_in(i)=0            !average pairwise RMSD including
         R_ex(i)=0            !average pairwise RMSD excluding
         do j=1,nstr_selected
            jj=str_selected(j)
            if(amat(jj,i_cl(i)).lt.rmsd_cut)then
               if(mark(jj).ne.0)then
                  n_str_cl_ex(i)=n_str_cl_ex(i)+1
                  R_ex(i)=R_ex(i)+amat(jj,i_cl(i))
                  i_str_cl_ex(i,n_str_cl_ex(i))=jj
               endif
               mark(jj)=0
               n_str_cl(i)=n_str_cl(i)+1
               i_str_cl(i,n_str_cl(i))=jj
               R_in(i)=R_in(i)+amat(jj,i_cl(i))
            endif
         enddo
         R_in(i)=R_in(i)/n_str_cl(i)
         R_ex(i)=R_ex(i)/n_str_cl_ex(i)
         rmsd_cl_cut(i)=rmsd_cut
         nc=i
c         write(*,'ai2,3(af5.2),ai4')
c     $        'i=',i,' rmsd_cl_cut(i)=',rmsd_cl_cut(i),
c     $        ' R_in(i)=',R_in(i),' R_ex(i)=',R_ex(i),
c     $        ' n_str_cl(i)=',n_str_cl(i)
      enddo
 41   continue    

      end !end of subroutine clusterize()


c     ###########################################
c     function random_rmsd provides with expected random rmsd for a loop
c     of given length and end-to-end distance, or for a terminus of
c     given length
c     input parameters are:
c     L: loop length
c     d: end-to-end distance. If negative, means the loop is a terminus
      function random_rmsd(L,d,islocal)
      logical islocal
      if(islocal.eq..true.)then !local rmsd
         if(L.lt.15)then
            random_rmsd=2+(4-2)/(15-5)*(L-5)
         else if(L.lt.30)then
            random_rmsd=4+(6.5-4)/(30-15)*(L-15)
         else
            random_rmsd=6.5
         endif
      else !global rmsd
         if(L.lt.15)then
            random_rmsd=3+(5-3)/(15-5)*(L-5)
         else if(L.lt.30)then
            random_rmsd=5+(8-5)/(30-15)*(L-15)
         else
            random_rmsd=8
         endif         
      endif
      end !end of function random_rmsd


c     ###########################################
c     subroutine make_combo produces combo 
c     and closc files for each unaligned region
c     either for local or global rmsd
      subroutine make_combo(iloop,islocal)
      parameter(ndim=1000)      !Length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(ntr=20000)       !number of trajectories
      parameter(ncl=100)            !number of clusters
      parameter(minloopL=4)     !minimum loop length to consider

      logical islocal
c     local variables to the subroutine      
      integer Lch1,i_tra,looplength,Lstem
      real energ
      real xt(ndim),yt(ndim),zt(ndim) !temporal coordinates

      common/trafiles/filen,n_tra
      character filen(ntr)*72
      integer n_tra

      common/whole/Lch,x,y,z
      integer Lch
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures      

      common/clusteriza/nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,
     $     i_str_cl_ex,rmsd_cl_cut
      integer nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,i_str_cl_ex
      dimension i_cl(ncl) !i_cl(i) is structure center of i'th cluster
      dimension n_str_cl(ncl)!of structures including some used struct in i'th cluster
      dimension n_str_cl_ex(ncl) !n_str_cl_ex(i): # of structures excluding some used struct
                                 ! in i'th cluster
      dimension i_str_cl(ncl,nst) !i_str_cl(i,j) is structure id of j'th structure belonging
                                  !to i'th cluster
      dimension i_str_cl_ex(ncl,nst) !i_str_cl_ex(i,j): structure id of j'th structure
                                     !belonging to i'th cluster
      real rmsd_cl_cut
      dimension rmsd_cl_cut(ncl) !rmsd_cl_cut(i)

      common/loop_pdb/loopxc,loopyc,loopzc,
     $     loopxclosc,loopyclosc,loopzclosc,rmsd_combo_closc
      real loopxc,loopyc,loopzc,loopxclosc,loopyclosc,loopzclosc,
     $     rmsd_combo_closc
      dimension loopxc(ncl,ncl,ndim),loopyc(ncl,ncl,ndim), !combo str. for unaligned reg
     $     loopzc(ncl,ncl,ndim)!loopxc(i,j,k) residue k of combo of cluster j for unaligned reg. i
      dimension loopxclosc(ncl,ncl,ndim),loopyclosc(ncl,ncl,ndim), !closc str. for unaligned reg
     $     loopzclosc(ncl,ncl,ndim)
      dimension rmsd_combo_closc(ncl) !rmsd between combo and closc

      dimension xs(ndim),ys(ndim),zs(ndim)
      common/r_/R_in,R_ex,Rc_in,Rc_ex
      dimension R_in(100),R_ex(100),Rc_in(100),Rc_ex(100)

      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      integer isalig(0:ndim),lbegin(ndim),lextent(ndim),nloop,hsc(ndim),
     &lextents(ndim),lresid(ndim,ndim)

c     nstr_selected: number of structures to do spicker on loops
c     str_selected(i): structure id for structure "i"
      common/dens/n_str_all,n_str,d_rc_in,drcinmax,cmax,nstr_selected,
     $     str_selected
      real d_rc_in(100),drcinmax
      integer n_str_all,n_str,cmax,nstr_selected,str_selected(nst)
ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000),
     $     r_4(3,ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /1000*1.0/

      do i=1,nc !cycle through the clusters
c     initialize centroid coordinates
         do ii=1,Lch
            xs(ii)=0
            ys(ii)=0
            zs(ii)=0
         enddo
c     retrieve coordinates for the center structure of the cluster
         k=i_cl(i) !structure id for center structure
         Rc_in(i)=0         !RMSD to centroid including

         looplength=lextent(iloop)
         do j=1,looplength      !cycle over local residues belonging to loop "iloop"
            l=lresid(j,iloop)   !residue number for one of the residues
            r_2(1,j)=x(l,k)   
            r_2(2,j)=y(l,k)
            r_2(3,j)=z(l,k)
         enddo
         if(islocal.eq..false.)then !if loop global rmsd
            Lstem=lextents(iloop)-looplength
            ik=1
            do j=looplength+1,lextents(iloop) !cycle over stem res. belonging to loop "iloop"
               l=lresid(j,iloop) !residue number for one of the residues
               r_4(1,ik)=x(l,k) !coordinates from the cluster center
               r_4(2,ik)=y(l,k)
               r_4(3,ik)=z(l,k)
               loopxc(iloop,i,j)=x(l,k) !also store in loopxc for future use
               loopyc(iloop,i,j)=y(l,k)
               loopzc(iloop,i,j)=z(l,k)
               ik=ik+1            
            enddo
         endif

         do l=1,n_str_cl(i) !cycle over structures of cluster "i"
            m=i_str_cl(i,l)     !l'th structure in i'th cluster
***   rotate nearby structure into the center structure---->
            do n=1,looplength   !cycle over residues in loop
               nn=lresid(n,iloop) 
               r_1(1,n)=x(nn,m)
               r_1(2,n)=y(nn,m)
               r_1(3,n)=z(nn,m)
            enddo
            if(islocal.eq..true.)then
               call u3b(w,r_1,r_2,looplength,1,rms,u,t,ier) !u rotate r_1 to r_2
               armsd=dsqrt(rms/looplength) !    write(*,*) armsd
            else !global loopp rmsd
               ik=1
               do n=looplength+1,lextents(iloop) !cycle over stem res. belonging to loop "iloop"
                  nn=lresid(n,iloop) !residue number for one of the residues
                  r_3(1,ik)=x(nn,m) !coordinates from the cluster center
                  r_3(2,ik)=y(nn,m)
                  r_3(3,ik)=z(nn,m)
                  ik=ik+1            
               enddo
               call u3b(w,r_3,r_4,Lstem,1,rms,u,t,ier) !align stems-->rot. matrix u
               armsd=residual(r_1,r_2,looplength,u,t) !calculate "rmsd" value and new "t"
c               write(*,'f5.2') armsd
            endif            
            Rc_in(i)=Rc_in(i)+armsd
            do j=1,looplength
               xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
               yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
               zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
               xs(j)=xs(j)+xx
               ys(j)=ys(j)+yy
               zs(j)=zs(j)+zz
            enddo
         enddo

         do j=1,looplength !normalize
            loopxc(iloop,i,j)=xs(j)/float( n_str_cl(i) )
            loopyc(iloop,i,j)=ys(j)/float( n_str_cl(i) )
            loopzc(iloop,i,j)=zs(j)/float( n_str_cl(i) )
         enddo
         Rc_in(i)=Rc_in(i)/float( n_str_cl(i) )
         d_rc_in(i)=n_str_cl(i)!initialize to big value (Rc_in(i)=1Anstrom)
         if(n_str_cl(i).gt.1) d_rc_in(i)=float(n_str_cl(i))/Rc_in(i) !related density
      enddo

c     calculate closc structure
      do i=1,nc
         rmsd_combo_closc(i)=1000.0
      enddo
      do i_tra=1,n_tra !cycle look through all trajectory files
         open(10,file=filen(i_tra),status='old')
 110     read(10,*,end=119,err=119)Lch1,energ!cycle through all snapshots in the trajectory file
         do i=1,Lch1
            read(10,*,end=119,err=119)xt(i),yt(i),zt(i)
         enddo

         do j=1,looplength      !select those coordinates of the unaliged region
            l=lresid(j,iloop)
            r_1(1,j)=xt(l)
            r_1(2,j)=yt(l)
            r_1(3,j)=zt(l)
         enddo
         if(islocal.eq..false.)then
            ik=1
            do j=looplength+1,lextents(iloop) !cycle over stem res. belonging to loop "iloop"
               l=lresid(j,iloop) !residue number for one of the residues
               r_3(1,ik)=xt(l)
               r_3(2,ik)=yt(l)
               r_3(3,ik)=zt(l)
               ik=ik+1
            enddo
         endif
         do i=1,nc !compare to each of the combo structures
            do j=1,looplength !select those coordinates of the unaliged region
               r_2(1,j)=loopxc(iloop,i,j)
               r_2(2,j)=loopyc(iloop,i,j)
               r_2(3,j)=loopzc(iloop,i,j)
            enddo
            if(islocal.eq..true.)then
               call u3b(w,r_1,r_2,looplength,0,rms,u,t,ier)
               armsd=dsqrt(rms/looplength)
            else
               ik=1
               do j=looplength+1,lextents(iloop) !cycle over stem res. belonging to loop "iloop"
                  r_4(1,ik)=loopxc(iloop,i,j) !these are actual coordinates of the cluster center
                  r_4(2,ik)=loopyc(iloop,i,j)
                  r_4(3,ik)=loopzc(iloop,i,j)
                  ik=ik+1
               enddo
               call u3b(w,r_3,r_4,Lstem,1,rms,u,t,ier) !align stems-->rot. matrix u
               armsd=residual(r_1,r_2,looplength,u,t) !calculate "rmsd" value and new "t"
            endif
            if(armsd.lt.rmsd_combo_closc(i))then
               rmsd_combo_closc(i)=armsd
               do j=1,looplength !update corresponding closc structure for loop residues
                  loopxclosc(iloop,i,j)=r_1(1,j)
                  loopyclosc(iloop,i,j)=r_1(2,j)
                  loopzclosc(iloop,i,j)=r_1(3,j)
               enddo
               if(islocal.eq..false.)then
                  ik=looplength+1
                  do j=1,Lstem !update corresponding closc structure for stem residues
                     loopxclosc(iloop,i,ik)=r_3(1,j)
                     loopyclosc(iloop,i,ik)=r_3(2,j)
                     loopzclosc(iloop,i,ik)=r_3(3,j)
                     ik=ik+1
                  enddo                  
               endif
            endif
         enddo
         goto 110 !read next snapshot
 119     close(10) !close current trajectory file
      enddo
         
c     determine the cluster with maximum density
      cmax=1              !initialize cluster with maximum Rc_in density
      drcinmax=d_rc_in(1) !initialize density of cluster with maximum density
      do i=2,nc
         if(d_rc_in(i).gt.drcinmax) then
            cmax=i
            drcinmax=d_rc_in(i)
         endif
      enddo
c      write(*,'ai2') 'cmax=',cmax

c     determine averages for excluding
      do i=1,nc
         do j=1,looplength !coordinates of centroid
            r_1(1,j)=loopxc(iloop,i,j)
            r_1(2,j)=loopyc(iloop,i,j)
            r_1(3,j)=loopzc(iloop,i,j)
         enddo
         if(islocal.eq..false.)then
            do j=1+looplength,lextents(iloop) !coordinates of centroid for stem residues
               r_3(1,j)=loopxc(iloop,i,j)
               r_3(2,j)=loopyc(iloop,i,j)
               r_3(3,j)=loopzc(iloop,i,j)
            enddo
         endif
         Rc_ex(i)=0             !RMSD to centroid excluding
         do l=1,n_str_cl_ex(i)  !cycle over structures of cluster "i"
            m=i_str_cl_ex(i,l)  !l'th structure in i'th cluster
            do n=1,looplength   !cycle over residues in loop
               nn=lresid(n,iloop) 
               r_2(1,n)=x(nn,m)
               r_2(2,n)=y(nn,m)
               r_2(3,n)=z(nn,m)
            enddo
            if(islocal.eq..true.)then
               call u3b(w,r_1,r_2,looplength,0,rms,u,t,ier) !u rotate r_1 to r_2
               armsd=dsqrt(rms/looplength) !RMSD12
            else
               ik=1
               do n=1+looplength,lextents(iloop) !cycle over stem residues in loop
                  nn=lresid(n,iloop) 
                  r_4(1,ik)=x(nn,m)
                  r_4(2,ik)=y(nn,m)
                  r_4(3,ik)=z(nn,m)
                  ik=ik+1
               enddo
               call u3b(w,r_3,r_4,Lstem,1,rms,u,t,ier)
               armsd=residual(r_1,r_2,looplength,u,t)
            endif
            Rc_ex(i)=Rc_ex(i)+armsd
         enddo
         Rc_ex(i)=Rc_ex(i)/float( n_str_cl(i) )
c         write(*,'2(af5.2)ai4af5.0')
c     $        'Rc_in(i)=',Rc_in(i),' Rc_ex(i)=',Rc_ex(i),
c     $        ' n_str_cl(i)=',n_str_cl(i),' d_rc_in(i)=',d_rc_in(i)         
      enddo

c      do i=1,nc
c         do j=1,lextent(iloop) !normalize
c            write(*,'3(x,3f8.3)') loopxc(iloop,i,j),loopyc(iloop,i,j),
c     $           loopzc(iloop,i,j)
c         enddo
c         write(*,*) ''
c      enddo

      end                       !end of subroutine make_combo


c     ###########################################
c     subroutine report_clustering produces an rst.dat file
c     for each unaligned region
      subroutine report_clustering(iloop,islocal)
      parameter(ndim=1000)      !Length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(ntr=20000)       !number of trajectories
      parameter(ncl=100)            !number of clusters
      parameter(minloopL=4)     !minimum loop length to consider

      logical islocal
c     local variables to the subroutine
      integer bl,cl,looplength,Lstem             !current length of buffer
      character *10000 b,c        !buffer
      character *32 rstfile
      real armsd,rmsd_combo,rmsd_combc,rmsd_model,rmsd_closc,
     $     rmsd_piece,density
      double precision u(3,3),t(3),rms,drms
      double precision r_1(3,ndim),r_2(3,ndim),w(1000),r_3(3,ndim),
     $     r_4(3,ndim)
      data w /1000*1.0/

c     global variables
      common/loop_pdb/loopxc,loopyc,loopzc,
     $     loopxclosc,loopyclosc,loopzclosc,rmsd_combo_closc
      real loopxc,loopyc,loopzc,loopxclosc,loopyclosc,loopzclosc,
     $     rmsd_combo_closc
       dimension loopxc(ncl,ncl,ndim),loopyc(ncl,ncl,ndim), !combo str. for unaligned reg
     $     loopzc(ncl,ncl,ndim)!loopxc(i,j,k) residue k of combo of cluster j for unaligned reg. i
      dimension loopxclosc(ncl,ncl,ndim),loopyclosc(ncl,ncl,ndim), !closc str. for unaligned reg
     $     loopzclosc(ncl,ncl,ndim)
      dimension rmsd_combo_closc(ncl) !rmsd between combo and closc

      common/globalclusters/nglobal,xc,yc,zc
      integer nglobal
      dimension xc(ncl,ndim),yc(ncl,ndim),zc(ncl,ndim) !combo structures

      common/whole/Lch,x,y,z
      integer Lch
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures

      common/native/Lch_n,x_n,y_n,z_n
      integer Lch_n
      dimension x_n(ndim),y_n(ndim),z_n(ndim) !native structure

      common/clusteriza/nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,
     $     i_str_cl_ex,rmsd_cl_cut
      integer nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,i_str_cl_ex

      dimension i_cl(ncl) !i_cl(i) is structure center of i'th cluster
      dimension n_str_cl(ncl)!of structures including some used struct in i'th cluster
      dimension n_str_cl_ex(ncl) !n_str_cl_ex(i): # of structures excluding some used struct
                                 ! in i'th cluster
      dimension i_str_cl(ncl,nst) !i_str_cl(i,j) is structure id of j'th structure belonging
                                  !to i'th cluster
      dimension i_str_cl_ex(ncl,nst) !i_str_cl_ex(i,j): structure id of j'th structure
                                     !belonging to i'th cluster
      real rmsd_cl_cut
      dimension rmsd_cl_cut(ncl) !rmsd_cl_cut(i)

      common/r_/R_in,R_ex,Rc_in,Rc_ex
      dimension R_in(100),R_ex(100),Rc_in(100),Rc_ex(100)

      common/dens/n_str_all,n_str,d_rc_in,drcinmax,cmax,nstr_selected,
     $     str_selected
      real d_rc_in(100),drcinmax
      integer n_str_all,n_str,cmax,nstr_selected,str_selected(nst)

      common/names/protein
      character protein*10

      common/chaindat/isalig,nloop,lbegin,lextent,lextents,lresid
      integer isalig(0:ndim),nloop,lbegin(ndim),lextent(ndim),
     &     lextents(ndim),lresid(ndim,ndim)

c     write header
      write(b,'a')'Target: '//protein//'\n'
      bl=lnblnk(b) !returns last non-blank character in a string
      write(b,'ai3a') b(1:bl)//'Modeling Length: ',0,'\n'
      bl=lnblnk(b)
      write(b,'ai4a') b(1:bl)//'Native Length: ',lextent(iloop),'\n'   
      bl=lnblnk(b)
      write(b,'ai2a') b(1:bl)//'Number of clusters: ',nc,'\n'
      bl=lnblnk(b)
      write(b,'ai6a') b(1:bl)//'Total number of structures=',
     $     n_str_all,'\n'
      bl=lnblnk(b)
      write(b,'ai5a') b(1:bl)//'Number of structure in use=',
     $     nstr_selected,'\n\n'

      looplength=lextent(iloop)
c     section comparison to native structure
      if(Lch_n.gt.1)then
         bl=lnblnk(b)
         write(b,'a') b(1:bl)//
     $        '--------- comparison to native structure---------\n'//
     $        '  i R_combo   R_combc  R_model R_closc  R_piece \n'//
     $        'A-------------------------------------------------\n'

         Lstem=lextents(iloop)-looplength
         do j=1,looplength      !native coordinates
            l=lresid(j,iloop)
            r_1(1,j)=x_n(l)     
            r_1(2,j)=y_n(l)
            r_1(3,j)=z_n(l)
         enddo
         if(islocal.eq..false.)then
            ik=1
            do j=1+looplength,lextents(iloop)   !native coordinates
               l=lresid(j,iloop)
               r_3(1,ik)=x_n(l)     
               r_3(2,ik)=y_n(l)
               r_3(3,ik)=z_n(l)
c               write(*,'3(f8.3,x)')r_3(1,ik),r_3(2,ik),r_3(3,ik)
               ik=ik+1
            enddo            
         endif

         do i=1,nc !cycle over the clusters
            do j=1,looplength !combo structure
               r_2(1,j)=loopxc(iloop,i,j)
               r_2(2,j)=loopyc(iloop,i,j)
               r_2(3,j)=loopzc(iloop,i,j)
            enddo
            if(islocal.eq..true.)then
               call u3b(w,r_1,r_2,looplength,0,rms,u,t,ier)
               rmsd_combo=dsqrt(rms/looplength)
            else
               ik=1
               do j=1+looplength,lextents(iloop) 
                  r_4(1,ik)=loopxc(iloop,i,j)
                  r_4(2,ik)=loopyc(iloop,i,j)
                  r_4(3,ik)=loopzc(iloop,i,j)
                  ik=ik+1
               enddo                   
               call u3b(w,r_3,r_4,Lstem,1,rms,u,t,ier)
               rmsd_combo=residual(r_1,r_2,looplength,u,t)
            endif
            do j=1,looplength !closc structure
               r_2(1,j)=loopxclosc(iloop,i,j)!residue "l" of combo "i" for unaligned region "iloop"
               r_2(2,j)=loopyclosc(iloop,i,j)
               r_2(3,j)=loopzclosc(iloop,i,j)
            enddo
            if(islocal.eq..true.)then
               call u3b(w,r_1,r_2,looplength,0,rms,u,t,ier)
               rmsd_closc=dsqrt(rms/looplength)
            else
               ik=1
               do j=1+looplength,lextents(iloop) !native coordinates
                  r_4(1,ik)=loopxclosc(iloop,i,j)
                  r_4(2,ik)=loopyclosc(iloop,i,j)
                  r_4(3,ik)=loopzclosc(iloop,i,j)
                  ik=ik+1
               enddo    
               call u3b(w,r_3,r_4,Lstem,1,rms,u,t,ier)
               rmsd_closc=residual(r_1,r_2,looplength,u,t)
            endif
            rmsd_combc=-1.0!for the moment we do not calculate rmsd to other average structures
            rmsd_model=-1.0
            rmsd_piece=-1.0

            bl=lnblnk(b)        !output to buffer
c            write(*,*) i,rmsd_combo,rmsd_combc,
c     $           rmsd_model,rmsd_closc,rmsd_piece
            write(b,'a,i5,5f8.3,a')
     $           b(1:bl)//'',i,rmsd_combo,rmsd_combc,
     $           rmsd_model,rmsd_closc,rmsd_piece,'\n'            
         enddo
         bl=lnblnk(b)
         write(b,'a')b(1:bl)//'\n'
      endif


c     section summary of clusters'
      bl=lnblnk(b)
      write(b,'a') b(1:bl)//
     $     ' ------------ summary of clusters -----------\n'//
     $     ' i Size R_cut density\n'//
     $     ' B--------------------------------------------------\n'
      do i=1,nc
         k=i_cl(i)
         density=float(n_str_cl(i))/rmsd_cl_cut(i)
         bl=lnblnk(b)
         write(b,132)b(1:bl)//'',i,n_str_cl(i),rmsd_cl_cut(i),
     $        density,'\n'
      enddo
 132  format(a,i2,i6,f6.2,f9.2,a)

c     section on include and exclude averages
      bl=lnblnk(b)
      write(b,'a') b(1:bl)//
     $     '  include used structure  exclude used structre\n'//
     $     '   ---------------------  ---------------------\n'//
     $     ' i  N_in  <R_in> <Rc_in>   N_ex  <R_ex> <Rc_ex>\n'//
     $     ' C---------------------------------------------\n'


 133  format(a,i4,i6,f7.2,f7.2,i6,f7.2,f7.2,a)
      do i=1,nc
         bl=lnblnk(b)
         write(b,133)b(1:bl)//'',i,n_str_cl(i),R_in(i),Rc_in(i),
     $        n_str_cl_ex(i),R_ex(i),Rc_ex(i),'\n'
      enddo

c     section on similarity between combo structures for the loop region
c     and loop region from the global combo structures
      bl=lnblnk(b)
      write(b,'a')b(1:bl)//
     $     '\nsimilarity between combo loop structures and'//
     $     ' loop region from the global combo structures\n'//
     $     'horizontal-axis: combo loop, vertical-axis: global combo\n'
      do i=1,nglobal !cycle over global combo structures
         write(c,'i7')i
         do j=1,looplength      !combo structure
            l=lresid(j,iloop)
            r_1(1,j)=xc(i,l)
            r_1(2,j)=yc(i,l)
            r_1(3,j)=zc(i,l)
         enddo
         if(islocal.eq..false.)then
            ik=1
            do j=1+looplength,lextents(iloop)   !combo structure
               l=lresid(j,iloop)
               r_3(1,ik)=xc(i,l)     
               r_3(2,ik)=yc(i,l)
               r_3(3,ik)=zc(i,l)
               ik=ik+1
            enddo            
         endif
         do j=1,nc !cycle over combo loop regions
            do k=1,looplength !combo structure
               r_2(1,k)=loopxc(iloop,j,k)
               r_2(2,k)=loopyc(iloop,j,k)
               r_2(3,k)=loopzc(iloop,j,k)
            enddo
            if(islocal.eq..true.)then
               call u3b(w,r_1,r_2,looplength,0,rms,u,t,ier)
               armsd=dsqrt(rms/looplength)
            else
               ik=1
               do k=1+looplength,lextents(iloop) !native coordinates
                  r_4(1,ik)=loopxc(iloop,i,k)
                  r_4(2,ik)=loopyc(iloop,i,k)
                  r_4(3,ik)=loopzc(iloop,i,k)
                  ik=ik+1
               enddo    
               call u3b(w,r_3,r_4,Lstem,1,rms,u,t,ier)
               armsd=residual(r_1,r_2,looplength,u,t)
            endif
            cl=lnblnk(c)
            write(c,'a,f6.2')c(1:cl)//' ',armsd
         enddo
         bl=lnblnk(b)
         cl=lnblnk(c)
         write(b,'a')b(1:bl)//c(1:cl)//'\n'
      enddo

c     section on similarity between combo structures for the loop region
      bl=lnblnk(b)
      write(b,'a')b(1:bl)//
     $     '\nsimilarity between combo loop structures\n'
      write(c,'a')'-------'
      do i=1,nc
         cl=lnblnk(c)
         write(c,'a,i6')c(1:cl)//' ',i
      enddo
      bl=lnblnk(b)
      cl=lnblnk(c)
      write(b,'a')b(1:bl)//c(1:cl)//'\n'
      
      do i=1,nc-1 !cycle over clusters
c     insert as many blank spaces as neccessary
         cl=0
         write(c,'i7')i
         do j=1,i
            cl=cl+7
            write(c,'a')c(1:cl)//'       '
         enddo
c     calculate the rmsd values between combos
         do j=1,looplength  !combo structure
            r_1(1,j)=loopxc(iloop,i,j)
            r_1(2,j)=loopyc(iloop,i,j)
            r_1(3,j)=loopzc(iloop,i,j)
         enddo         
         if(islocal.eq..false.)then
            ik=1
            do j=1+looplength,lextents(iloop)   !combo structure
               r_3(1,ik)=loopxc(iloop,i,j)     
               r_3(2,ik)=loopyc(iloop,i,j)
               r_3(3,ik)=loopzc(iloop,i,j)
               ik=ik+1
            enddo            
         endif
         do j=i+1,nc !cycle over clusters
            do k=1,looplength !combo structure
               r_2(1,k)=loopxc(iloop,j,k)
               r_2(2,k)=loopyc(iloop,j,k)
               r_2(3,k)=loopzc(iloop,j,k)
            enddo
            if(islocal.eq..true.)then
               call u3b(w,r_1,r_2,looplength,0,rms,u,t,ier)
               armsd=dsqrt(rms/looplength)
            else
               ik=1
               do k=1+looplength,lextents(iloop) !native coordinates
                  r_4(1,ik)=loopxc(iloop,j,k)
                  r_4(2,ik)=loopyc(iloop,j,k)
                  r_4(3,ik)=loopzc(iloop,j,k)
                  ik=ik+1
               enddo    
               call u3b(w,r_3,r_4,Lstem,1,rms,u,t,ier)
               armsd=residual(r_1,r_2,looplength,u,t)
            endif
            cl=cl+7
            write(c,'a,f5.2')c(1:cl)//'  ',armsd
         enddo
       bl=lnblnk(b)
       cl=lnblnk(c)
       write(b,'a')b(1:bl)//c(1:cl)//'\n'       
      enddo

c     rst_XX.dat file name
      if(islocal.eq..true.)then
         if(iloop.lt.10)then
            write(rstfile,'ai1a')'rst_local_0',iloop,'.dat'
         else
            write(rstfile,'ai2a')'rst_local_',iloop,'.dat'
         endif
      else
         if(iloop.lt.10)then
            write(rstfile,'ai1a')'rst_global_0',iloop,'.dat'
         else
            write(rstfile,'ai2a')'rst_global_',iloop,'.dat'
         endif
      endif

c     write buffer to file
      open(unit=20,file=rstfile,status='unknown')
      bl=lnblnk(b)
      write(20,'a') b(1:bl)//'\n'
      close(20)

      end


c     ###########################################
c     subroutine select_structures initializes str_selected, which
c     identifies the structures belonging to most dense cluster
      subroutine select_structures(isloop)
      parameter(ndim=1000)      !Length
      parameter(nst=13000)      !number of used structure, maximum allowed
      parameter(ntr=20000)       !number of trajectories
      parameter(ncl=100)            !number of clusters
      parameter(minloopL=4)     !minimum loop length to consider
      
      logical isloop !.true. if we will do spicker on an unaligned region

      common/dens/n_str_all,n_str,d_rc_in,drcinmax,cmax,nstr_selected,
     $     str_selected
      real d_rc_in(100),drcinmax
      integer n_str_all,n_str,cmax,nstr_selected,str_selected(nst)
      
      common/clusteriza/nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,
     $     i_str_cl_ex,rmsd_cl_cut
      integer nc,i_cl,n_str_cl,i_str_cl,n_str_cl_ex,i_str_cl_ex
      dimension i_cl(ncl) !i_cl(i) is structure center of i'th cluster
      dimension n_str_cl(ncl)!of structures including some used struct in i'th cluster
      dimension n_str_cl_ex(ncl) !n_str_cl_ex(i): # of structures excluding some used struct
                                 ! in i'th cluster
      dimension i_str_cl(ncl,nst) !i_str_cl(i,j) is structure id of j'th structure belonging
                                  !to i'th cluster
      dimension i_str_cl_ex(ncl,nst) !i_str_cl_ex(i,j): structure id of j'th structure
                                     !belonging to i'th cluster
      real rmsd_cl_cut
      dimension rmsd_cl_cut(ncl) !rmsd_cl_cut(i)

      if(isloop.eq. .true.)then !select structures belonging to most dense cluster
         nstr_selected=n_str_cl(cmax)
         do i=1,nstr_selected
            str_selected(i)=i_str_cl(cmax,i)
c           write(*,*)str_selected(i)
         enddo
c        write(*,'ai2ai5') 'cmax=',cmax,' nstr_selected=',nstr_selected
      else !select all structures
         nstr_selected=n_str
         do i=1,nstr_selected
            str_selected(i)=i
         enddo
      endif

      end                       !end of subroutine select_structures
c     ###########################################
      
 
