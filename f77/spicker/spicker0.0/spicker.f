ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This program, Structure-PCIKER (SPICKER), is aimed at selecting the 
c     best fold by clustering trajectors produced by CABS model. 
c     Comments and bug reports should be addressed to zhang6@buffalo.edu.
c
c     Input files includes:
c       'rmsinp'---Mandatory, length of protein & piece for RMSD calculation;
c       'seq.dat'--Mandatory, sequence file, for output of PDB models.
c       'tra.in'---Mandatory, list of trajectory names used for clustering.
c       files in 'tra.in'---Mandatory, trajectories of structures.
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
c      parameter(nst=100)      !number of used structure, maximum allowed
      parameter(ntr=20000)       !number of trajectories
      parameter(ncl=100)            !number of clusters
      character filen(ntr)*72,c2*6,seq(ndim)*3,protein*10
      character txt1*70,txt2*70
ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(1000)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /1000*1.0/
ccc   
      dimension xtemp(ndim),ytemp(ndim),ztemp(ndim) !structure close to templ
      dimension xt(ndim),yt(ndim),zt(ndim) !temporal coordinate
      dimension x_n(ndim),y_n(ndim),z_n(ndim) !native structure
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures
      dimension amat(nst,nst)    !RMSD matrics
      dimension mark(nst)       !for removing used structures
      dimension n_str_near(nst) !numbef of neighboring structures
      dimension itra(nst),istr(nst),E(nst)
      dimension n_str_cl(ncl),n_str_cl_ex(ncl)
      dimension i_cl(ncl),rmsd_cl_cut(ncl)
      dimension E_combo(ncl)
      dimension i_str_cl(ncl,nst),i_str_cl_ex(ncl,nst)
      dimension rmsd_str(nst)
      dimension xs(ndim),ys(ndim),zs(ndim)
      dimension xc(ncl,ndim),yc(ncl,ndim),zc(ncl,ndim)
      dimension xcl(ncl,ndim),ycl(ncl,ndim),zcl(ncl,ndim)
      dimension xc3(ncl,ndim),yc3(ncl,ndim),zc3(ncl,ndim)
      dimension rmsd_close_min(ncl) !RMSD between 'combo.pdb' and 'closm.pdb'
      dimension i1(100),i2(100)
      integer q(1000)
      dimension R_in(100),R_ex(100),Rc_in(100),Rc_ex(100)
      dimension order(nst)
      dimension ires(ndim)
      
**********************************************************************
****  input:
      open(1,file='rmsinp',status='old') !read Lch
      open(2,file='seq.dat',status='old') !read SEQ for model.pdb
      open(3,file='CA',status='unknown') !for comparison to native
      open(4,file='tra.in',status='old') !information for trajectories
      open(5,file='TEMP',status='unknown') !for comparison to template

****  output:
      open(20,file='rst.dat',status='unknown')
      open(21,file='str.txt',status='unknown') !structures in cluster
      open(22,file='RMSD.list',status='unknown') !structures in cluster
**********************************************************************

      read(1,*)n1,n2            !calculate RMSD between [n1,n2]
      read(1,*)Lch
      read(1,*)protein

******** following parameter has been optimized ********************
      RMSD_cut_initial=8
      RMSD_cut_min=3.5
      RMSD_cut_max=12
      ratio1=0.7
      ratio2=0.15
      nc_max=10
********************************************************************

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
***   The total number of structure will be used for picking clustering struct.
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
 11   read(3,1239,end=5)tx,ii,tx,tx,ii,tx,a1,a2,a3
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

cccccccccccccc calculate RMSD matrics ccccccccccccccccccccccccccc
c      write(*,*)'number of used structures=',n_str
      do i=1,n_str
         mark(i)=1
         do j=i,n_str
            do k=1,Lch
               r_1(1,k)=x(k,i)
               r_1(2,k)=y(k,i)
               r_1(3,k)=z(k,i)
               r_2(1,k)=x(k,j)
               r_2(2,k)=y(k,j)
               r_2(3,k)=z(k,j)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch) !RMSD12
            amat(i,j)=armsd
            amat(j,i)=armsd
         enddo
      enddo
c^^^^^^^^^^^^ RMSD matrics finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**************************************************************
*     find out the biggest cluster and decide the RMSD_cut_min
**************************************************************
      RMSD_cut=RMSD_cut_initial
 140  do j=1,n_str
         n_str_near(j)=0        !number of structure in cluster
         do k=1,n_str
            if(amat(j,k).lt.RMSD_cut)then
               n_str_near(j)=n_str_near(j)+1
            endif
         enddo
      enddo
***   find out the biggest cluster ----------------->
      n_str_cl_max=0            !maximum number of structures in cluster
      do j=1,n_str
         if(n_str_near(j).gt.n_str_cl_max)then
            n_str_cl_max=n_str_near(j)
         endif
      enddo
***   check the size of the cluster ------------->
      ratio=float(n_str_cl_max)/float(n_str)
      if(ratio.gt.ratio1.and.RMSD_cut.gt.RMSD_cut_min)then
         RMSD_cut=RMSD_cut-0.1
         goto 140
      endif
      if(ratio.lt.ratio2.and.RMSD_cut.lt.RMSD_cut_max)then
         RMSD_cut=RMSD_cut+0.2
         goto 140
      endif
      if(RMSD_cut.lt.RMSD_cut_min)RMSD_cut=RMSD_cut_min
      if(RMSD_cut.gt.RMSD_cut_max)RMSD_cut=RMSD_cut_max
*****RMSD_cut of first cluster, i.e. RMSD_cut_min, is decided ***

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc find out nc_max clusters cccccccccccccccccccccccc      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nc=0                      !number of clusters
      do i=1,nc_max
***   calculate nc(j), number of neightboring structures ----------------->
         do j=1,n_str
            n_str_near(j)=0     !number of structure in cluster
            do k=1,n_str
               if(mark(j).eq.1.and.mark(k).eq.1)then
                  if(amat(j,k).lt.RMSD_cut)then
                     n_str_near(j)=n_str_near(j)+1
                  endif
               endif
            enddo
         enddo
***   find out the biggest cluster ----------------->
         n_str_cl_max=0         !maximum number of structures in cluster
         do j=1,n_str
            if(n_str_near(j).gt.n_str_cl_max)then
               n_str_cl_max=n_str_near(j)
               i_cl(i)=j        !structure center of i'th cluster
            endif
         enddo
***   check the size of the cluster ------------->
         if(n_str_cl_max.lt.1) goto 41 !no remaining clusters.
***   remove the structures in the cluster-------->
         n_str_cl(i)=0          !# of structure including some used struct
         n_str_cl_ex(i)=0       !# of structure excluding used structure
         R_in(i)=0            !average pairwise RMSD including
         R_ex(i)=0            !average pairwise RMSD excluding
         do j=1,n_str
c     if(mark(j).eq.1)then
            if(amat(j,i_cl(i)).lt.RMSD_cut)then
               if(mark(j).ne.0)then
                  n_str_cl_ex(i)=n_str_cl_ex(i)+1
                  R_ex(i)=R_ex(i)+amat(j,i_cl(i))
                  i_str_cl_ex(i,n_str_cl_ex(i))=j
               endif
               mark(j)=0
               n_str_cl(i)=n_str_cl(i)+1
               i_str_cl(i,n_str_cl(i))=j
               R_in(i)=R_in(i)+amat(j,i_cl(i))
            endif
c     endif
         enddo
         R_in(i)=R_in(i)/n_str_cl(i)
         R_ex(i)=R_ex(i)/n_str_cl_ex(i)
         rmsd_cl_cut(i)=rmsd_cut
         nc=i
c         write(*,*)i,n_str_cl_ex(i),n_str_cl_max
      enddo
c^^^^^^^^^^^^^^^^ 5 cluster find out ^^^^^^^^^^^^^^^^^^^^^^
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 41   continue

*****************************************************************
c     calculate 'combo.pdb' and E_combo
COMBO************************************************************
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
         do k=1,Lch
            r_1(1,k)=xc(i,k)
            r_1(2,k)=yc(i,k)
            r_1(3,k)=zc(i,k)
         enddo
         nn=0
         Rc_in(i)=0             !RMSD to centroid including
         do j=1,n_str_cl(i)
            m=i_str_cl(i,j)
            do k=1,Lch
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
      enddo
      do i=1,nc
         do k=1,Lch
            r_1(1,k)=xc(i,k)
            r_1(2,k)=yc(i,k)
            r_1(3,k)=zc(i,k)
         enddo
         nn=0
         Rc_ex(i)=0             !RMSD to centroid including
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

cccccccccccccc output the 5 clusters in PDB format ccccccccccccccc
      do i=1,Lch
         read(2,*)ii,seq(i)
      enddo
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
         density=n_str_cl(i)/rmsd_cl_cut(i)
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

      stop
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
