c        calculate par.dat (pair pot) comb.dat (contact restr ) for SG
c           dist.dat distL.dat (dist restr)  combCA.dat (contact restr ) for CA
	implicit real*8 (a-h,o-z)
        parameter (maxres=600,maxmd=200,iwin=9,maxp=20000,
     &  maxd=maxres*maxmd,maxdp=20000000,maxq=600)
        character whole*80,ctmp6*6,frgfil*255,exclfil*255,
     &  resn(20)*3,mypdb*10,rname*3,line*80,mytemp(maxp)*6
     &  ,exclist(maxp)*6,datbase*255,atmn*3,resnyz(20)*3
     &  ,tempin*255,temppdb(maxp)*6,whole1(maxres)*80
        integer ifrg(maxres,2),itemp(maxres,maxmd,2),myseq(maxres),
     &  tempid(maxp),tempres(maxp),tempnum(maxres),templst(maxres,maxmd)
     &  ,tempused(maxp),temppos(maxp),tempmap(maxd,2),excind(maxp),
     &  excnum,pdbselected(maxp),isort(maxp),ibase(maxdp),frgpos(maxdp),
     &  invmap(maxres),resmap(20),ndist(maxq,maxq),rnum,
     &  ndistA(maxq,maxq),tempseq(maxres),exclid(maxp)
        real*8  tempscore(maxres,maxmd),tempxyz(maxres,maxmd,iwin,3),
     &  cacnt(maxres,maxres),zscore(maxp),cacnt0(maxres,maxres)
     &  ,sgcut(20,20),cacut(20,20),sgcut_dev(20,20),cacut_dev(20,20)
     &  ,frgsc(maxdp),scmx(maxres,maxmd),scmy(maxres,maxmd),
     &   scmz(maxres,maxmd),cax(maxres,maxmd),cay(maxres,maxmd),
     &   caz(maxres,maxmd),cx,cy,cz,freq_comb(maxq,maxq),
     &   freq0_par(maxres,maxres),dist(maxq,maxq,maxmd),
     &   distA(maxq,maxq,maxmd),freq_combCA(maxq,maxq),ee(maxq,maxq)
     &   ,disave(maxp),disdelta(maxp),cax1(0:maxres),cay1(0:maxres),
     &    caz1(0:maxres),scmx1(maxres),scmy1(maxres),scmz1(maxres)
     &   ,pij(maxres,maxres)
         integer irestyp(maxres),jtemp(maxres), 
     &   mee(maxq,maxq),iid(maxp),jjd(maxp)
        data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/
        data resnyz/'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &  'MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE',
     &  'TYR','TRP'/

	 open(unit=10,file='contact.comm',status='old')
          read(10,*) 
          do i=1,20
             read(10,*) rname,(sgcut(i,j),j=1,20)                     
          enddo
          read(10,*) 
          do i=1,20
             read(10,*) rname,(sgcut_dev(i,j),j=1,20)                     
          enddo
          do i=1,20
              do j=1,20
              sgcut(i,j)=sgcut(i,j)+sgcut_dev(i,j)*2.5
              enddo
          enddo
          read(10,*) 
          do i=1,20
             read(10,*) rname,(cacut(i,j),j=1,20)                     
          enddo
          read(10,*) 
          do i=1,20
             read(10,*) rname,(cacut_dev(i,j),j=1,20)                     
          enddo
          do i=1,20
              do j=1,20
              cacut(i,j)=cacut(i,j)+cacut_dev(i,j)*2.5
              enddo
          enddo
          close(10)



               do i=1,maxq
                   do j=1,maxq
                    ndist(i,j)=0
                    ndistA(i,j)=0
                    freq_comb(i,j)=0.0
                    freq0_par(i,j)=0.0
                   pij(i,j)=0.0
                   enddo
               enddo
               pcut=0.20


        read(*,*) ntemp,irest


              ib0=0 
c           open(unit=40,file='seq.dat',status='old') 
c	read in template  coordinates:
c              do i=1,irest 
c              read(40,*) j,rname
c              jtemp(i)=0
c                 do k=1,20
c                 if(rname.eq.resn(k))   irestyp(i)=k
c                 enddo
c              enddo
c             close(40)
        do im=1,ntemp
          

c           read(*,'(a255)' ) tempin
           read(*,* ) tempin
           open(unit=10,file=tempin,status='old')
           ibs=index(tempin,' ')-1
        write(*,*) 'reading ',im,' ', tempin(1:ibs)
           rewind(10)
           do i=1,irest
           frgpos(i+ib0)=-1
           cax1(i)=0
           cay1(i)=0
           caz1(i)=0
           enddo
        id=1
        ires0=-100
        rnum=0
15      read(10,'(a80)',end=10) whole1(id)
        if(whole1(id)(1:4).ne.'ATOM') goto 15
        if(whole1(id)(14:15).ne.'CA') goto 15  ! CA only

        READ (whole1(id)(23:26),'(i4)') ires1
        if(ires1.ne.ires0) then
             if(ires0.lt.0) ist=ires1
             ires0=ires1             
             rname=whole1(id)(18:20)
             do ll=1,20
              if(resnyz(ll).eq.rname) tempseq(ires1)=ll
              if(resn(ll).eq.rname) irestyp(ires1)=ll
             enddo
            frgpos(ires1+ib0)=ires1          
            jtemp(ires1)=jtemp(ires1)+1
        endif             

             READ (whole1(id)(31:54),'(3F8.3)') Xx,Yy,Zz
c            write(*,*) ires1,xx,yy,zz
              cax1(ires1)=xx
              cay1(ires1)=yy
              caz1(ires1)=zz
                id=id+1
             goto 15
 10	close(10)
            iend=ires0

        write(*,*) ist,iend
c          do i=ist,iend
c         write(*,*) i, irestyp(i)
c          enddo
        call getsg(irest,ist,iend,irestyp,cax1,cay1,caz1,scmx1,
     &  scmy1,scmz1,frgpos)          
	 
c         do i=ist,iend
c         write(*,*) cax1(i),cay1(i),caz1(i),' CA ',resn(irestyp(i))
c         write(*,*) scmx1(i),scmy1(i),scmz1(i),' SG'
c         enddo

c       SG contact and pair number's
                 do l=1,irest-1
                    l1=frgpos(l+ib0)
                    do j=l+1,irest
                    j1=frgpos(j+ib0)
                    if(l1.gt.0.and.j1.gt.0) then
                    dis=sqrt((scmx1(l)-scmx1(j))**2
     &                      +(scmy1(l)-scmy1(j))**2
     &                      +(scmz1(l)-scmz1(j))**2)
                    freq0_par(l,j)=freq0_par(l,j)+1
                    if(dis.lt.sgcut(tempseq(l1),tempseq(j1))) then
                    freq_comb(l,j)=freq_comb(l,j)+1
                    endif

                    endif
                    enddo
                 enddo



c       CA contact and distances
                 do l=1,irest-1
                    l1=frgpos(l+ib0)
                    do j=l+1,irest
                    j1=frgpos(j+ib0)
                    if(l1.gt.0.and.j1.gt.0) then
                    dis=sqrt((cax1(l)-cax1(j))**2
     &                      +(cay1(l)-cay1(j))**2
     &                      +(caz1(l)-caz1(j))**2)
                    ndista(l,j)=ndista(l,j)+1
                    dista(l,j, ndista(l,j))=dis
c       short range distance:
                    if((j-l).le.6) then
                    do iu=l,j
                    if(frgpos(iu+ib0).lt.0) goto 301 
                    enddo

                    do iu=l+1,j
                    dis2=sqrt((cax1(iu-1)-cax1(iu))**2
     &                      +(cay1(iu-1)-cay1(iu))**2
     &                      +(caz1(iu-1)-caz1(iu))**2)
                    if(dis2.gt.4.1) goto 301 
                    enddo
                    

                    ndist(l,j)=ndist(l,j)+1
                    dist(l,j, ndist(l,j))=dis

                    endif
c       end if short range
 301		    continue

                    endif
                    enddo
                 enddo

            enddo  ! im
c           do i=1,irest
c             jtemp(i)=ntemp
c           enddo 
c	    ntemp=ntemp/nchunk

c       output .dat files:

c       output comb.dat for SG
           nres=0
           do i=1,irest-5
              do j=i+5,irest
c              if(jtemp(i).gt.0.1.and.jtemp(j).gt.0.1) then
c               pij(i,j)=freq_comb(i,j)/sqrt(float(jtemp(i)*jtemp(j)))
c              endif
         if(freq0_par(i,j).gt.0) 
     &    pij(i,j)=freq_comb(i,j)/freq0_par(i,j)
c           if(i.eq.1) write(*,*) freq_comb(i,j),' nnnn'
              if(pij(i,j).gt.pcut) nres=nres+1
              enddo
           enddo
           open(unit=20,file='comb.dat',status='replace')
           write(20,*) nres
           do i=1,irest-5
              do j=i+5,irest
              if(pij(i,j).gt.pcut) 
     &     write(20,130) i,j,freq_comb(i,j)
c     &   /sqrt(float(jtemp(i)*jtemp(j)))
     &   /freq0_par(i,j)
              enddo
 	   enddo   

 130	   format(5x,i5,3x,i5,3x,f5.3)
           close(20)

c       output par.dat for SG:
c          xtemp=1./ntemp*0.8
          ee_a=0
          n_ee=0
        do i=1,irest-1
           do j=i+1,irest
           mee(i,j)=-1
           if(freq0_par(i,j).gt.0) then


           frequence=freq_comb(i,j)/(freq0_par(i,j)
     &     /(0.8*sqrt(float(jtemp(i)*jtemp(j)))))
c          write(*,*) freq0_par(i,j),' freq0 ',frequence
           if(frequence.gt.0.0001) then
           ee(i,j)=-log(frequence)
c            write(*,*) -log(frequence),frequence
           else
           ee(i,j)=0.0
           endif    
           n_ee=n_ee+1
           mee(i,j)=1
           ee_a=ee_a+ee(i,j)
           endif

           enddo
        enddo
c       shift ee
         ee_a=ee_a/n_ee
	 do i=1,irest-1
              do j=i+1,irest
              if(mee(i,j).gt.0) then
                ee(i,j)=ee(i,j)-ee_a
              else
                ee(i,j)=0.0
              endif
              ee(j,i)=ee(i,j)
              enddo
         enddo
           open(unit=20,file='par.dat',status='replace')
           do i=1,irest
           write(20,'(5x,i5,1x,15a)') i,"==============="
           write(20,'(1x, 10f8.3)')(ee(i,j),j=1,irest)
           enddo
           close(20)
               
              

c       output combCA.dat for CA
           nres=0
           do i=1,irest-5
              do j=i+5,irest
                 nn=ndistA(i,j)
                 freq_combCA(i,j)=0
                do k=1,nn
                 if(distA(i,j,k).lt.6.0) then 
                     freq_combCA(i,j)=freq_combCA(i,j)+1
                 endif
                enddo
c               if(jtemp(i).gt.0.1.and.jtemp(j).gt.0.1) then
c                pij(i,j)=freq_combCA(i,j)/sqrt(float(jtemp(i)*jtemp(j)))
c               endif
               if(freq0_par(i,j).gt.0)
     &       pij(i,j)=freq_combCA(i,j)/freq0_par(i,j)
                if(pij(i,j).gt.pcut) nres=nres+1
              enddo
           enddo

           open(unit=20,file='combCA.dat',status='replace')
           write(20,*) nres
           do i=1,irest-5
              do j=i+5,irest
              if(pij(i,j).gt.pcut) 
     &     write(20,130) i,j,freq_combCA(i,j)
c     &     /sqrt(float(jtemp(i)*jtemp(j)))
     &     /freq0_par(i,j)
              enddo
 	   enddo   

           close(20)

c       output dist.dat  for CA
         nres=0
         xdist0=2
         if(ntemp.lt.2) xdist0=1
         do i=1,irest-1
            jm=i+6
            if(jm.gt.irest) jm=irest       
            do j=i+2,jm
              kk0=ndist(i,j)
c            write(*,*) kk0,' kk0'
              if(kk0.ge.xdist0) then
              dist_a=0
              dist2_a=0
                 do k=1,kk0
                     dist_a=dist_a+dist(i,j,k)   
                     dist2_a=dist2_a+dist(i,j,k)**2   
                 enddo
                  dist_a=dist_a/kk0
                  dist2_a=dist2_a/kk0
                  delta2=dist2_a-dist_a**2
                  if(delta2.lt.0.00001) then
                     delta=0.0
                  else          
                     delta=sqrt(delta2)
                  endif
                  nres=nres+1
                  iid(nres)=i                  
                  jjd(nres)=j                  
                  disave(nres)=dist_a
                  disdelta(nres)=delta

              endif
            enddo
         enddo

           open(unit=20,file='dist.dat',status='replace')
           write(20,'(5x,i5)') nres
           do i=1,nres
           write(20,'(1x,i5,1x,i5,1x,i5,1x,f8.3,1x,f8.3)') 
     &     iid(i),jjd(i),ndist(iid(i),jjd(i)),disave(i),disdelta(i)
           enddo
           close(20)
           

c       output distL.dat  for CA
           nres=0
           do i=1,irest-10
              do j=i+10,irest,10
                 nn=ndistA(i,j)
                do k=1,nn
                     nres=nres+1
                enddo
              enddo
           enddo

           open(unit=20,file='distL.dat',status='replace')
           write(20,*) nres
           do i=1,irest-10
              do j=i+10,irest,10
                 nn=ndistA(i,j)
               do k=1,nn
           write(20,140) i,j,distA(i,j,k)
               enddo
              enddo
 	   enddo   
 140	   format(5x,i5,3x,i5,3x,f7.3)
           close(20)




	stop
	end
          
c	###########################################
	subroutine setcacnt(mypdb,cacnt,zs_i)
        parameter(maxres=600)
	real*8 cacnt(maxres,maxres),zs_i
        character mypdb*6
           

            ipos=index(mypdb,' ')-1
            if(ipos.le.0) ipos=6
            open(unit=20,file=mypdb(1:ipos)//'_sp3temp.cnt',
     &      status='old')

10	    read(20,*,end=50) i1,i2,ca
            cacnt(i1,i2)=zs_i
            cacnt(i2,i1)=zs_i
            goto 10
50	close(20)
           return
           end
c     ########################################
          subroutine readexc(excfil,exclist,excnum,excind)
          character excfil*255,exclist(*)*6
          integer excnum,excind(*)

         open(unit=70,file=excfil,status='old')
         excnum=1
40        read(70,'(a6)',end=50) exclist(excnum)
          idd=index(exclist(excnum),' ')-1
          if(idd.le.0) idd=6
         excind(excnum)=idd
         excnum=excnum+1
         goto 40
50      close(70)
              excnum=excnum-1

          return
        end

***********************************************************************
        subroutine sortzscore(npro,zscore,zscore_sort)
        real*8 zscore(*)
        integer zscore_sort(*)
c       use unix command sort

                open(unit=50,file='_sort.in',status='replace')
                   do i=1,npro
                   write(50,*) i,zscore(i)
                   enddo
                close(50)
                call system('sort -nr +1 _sort.in > _sort.out')
                open(unit=50,file='_sort.out',status='old')
                do i=1,npro
                 read(50,*) zscore_sort(i)
                enddo
                close(50)



         return
        end


c       #########################
          subroutine topnzscore(npro,zscore,isort,nn)
        parameter (maxres=600,maxp=20000)
           real*8 zscore(*),xx,yy
           integer isort(*),idel(maxp)
           do i=1,npro
            idel(i)=i
           enddo

	   do i=1,nn
                 isort(i)=i
                 xx=zscore(i)
                 do j=1,npro
                 if(idel(j).gt.0) then
                    yy=zscore(j)
                 if(yy.gt.xx) then 
                    j0=j
                    xx=yy
                 endif   
                 endif
                 enddo  
                 isort(i)=j0
                 idel(j0)=-1

           enddo

          return
	  end


c###########################################
        subroutine getsg(rnum,ist,iend,irestyp,cax,cay,caz,scmx,
     &  scmy,scmz,frgpos)
         parameter(maxa=100000,ntyp=20,maxres=600,nvec=312)
         integer rnum, seq(maxres),irestyp(*),frgpos(*)
         real*8 cax(0:maxres),cay(0:maxres),caz(0:maxres),
     &   scmx(*),scmy(*),scmz(*)
         character base*250,aa(-1:20)*3,resn(20)*3
           data init/0/
      data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
c                -1    0     1     2     3     4     5     6
     &     'PRO','MET','ASP','ASN','LEU',
c            7     8     9    10    11
     &     'LYS','GLU','GLN','ARG',
c           12    13    14    15
     &     'HIS','PHE','TYR','TRP','CYX'/
c           16    17    18    19    20

      common/beta3/seq 

      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)

             if(init.eq.0) then
             init=1
c     read in SG data:
      open(unit=4,
c     & file=base(1:ibs)//'/data/sidechain.comm',
     & file='sidechain.comm',
     & status='old') !for Sc position.

      do k=0,19
         read(4,*)
         read(4,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
      enddo
      CLOSE(4)                  !sidecent.comm



        endif  ! init


             do i=ist,iend
               do l=0,19
                if(resn(irestyp(i)).eq.aa(l)) seq(i)=l
               enddo 
             enddo    


c     side-chain coordinate from C_a to side-chain ---------------->

       
         cax(0)=cax(1)+(cax(2)-cax(3))
         cay(0)=cay(1)+(cay(2)-cay(3))
         caz(0)=caz(1)+(caz(2)-caz(3))


         cax(rnum+1)=cax(rnum)+(cax(rnum-1)-cax(rnum-2))         
         cax(rnum+1)=cay(rnum)+(cay(rnum-1)-cay(rnum-2))         
         cax(rnum+1)=caz(rnum)+(caz(rnum-1)-caz(rnum-2))         

         do i=ist,iend   
c         write(*,*) cax(i),cay(i),caz(i),' CA ',i
           j=i
           jm=i-1
           jp=i+1
           amx=cax(j)-cax(jm)
           amy=cay(j)-cay(jm)
           amz=caz(j)-caz(jm)

           aaa=sqrt(amx**2+amy**2+amz**2)
           amx=amx/aaa      
           amy=amy/aaa
           amz=amz/aaa
           apx=cax(jp)-cax(j)
           apy=cay(jp)-cay(j)
           apz=caz(jp)-caz(j)
           aaa=sqrt(apx**2+apy**2+apz**2)
           apx=apx/aaa      
           apy=apy/aaa
           apz=apz/aaa

           ang=acos(-(amx*apx+amy*apy+amz*apz))*180/3.1415926 !


              aaax=amx+apx
              aaay=amy+apy
              aaaz=amz+apz
              aaa=sqrt(aaax**2+aaay**2+aaaz**2)
              ax=aaax/aaa
              ay=aaay/aaa
              az=aaaz/aaa


               ccx=amx-apx
               ccy=amy-apy
               ccz=amz-apz
               aaa=sqrt(ccx**2+ccy**2+ccz**2)
               cx=ccx/aaa
               cy=ccy/aaa
               cz=ccz/aaa


               bbx=amy*apz-amz*apy
               bby=amz*apx-amx*apz
               bbz=amx*apy-amy*apx
               aaa=sqrt(bbx**2+bby**2+bbz**2)
               bx=bbx/aaa
               by=bby/aaa
               bz=bbz/aaa

            k=seq(i)
            if(ang.lt.105) then ! alpha-helix or turn like
              dx=(axalf(k)*ax+ayalf(k)*bx+azalf(k)*cx)
              dy=(axalf(k)*ay+ayalf(k)*by+azalf(k)*cy)
              dz=(axalf(k)*az+ayalf(k)*bz+azalf(k)*cz)
            else                ! beta-sheet
              dx=(axbet(k)*ax+aybet(k)*bx+azbet(k)*cx)
              dy=(axbet(k)*ay+aybet(k)*by+azbet(k)*cy)
              dz=(axbet(k)*az+aybet(k)*bz+azbet(k)*cz)
            endif


           scmx(i)=cax(j)+dx 
           scmy(i)=cay(j)+dy
           scmz(i)=caz(j)+dz 

         enddo ! i->rnum


           return
           end
