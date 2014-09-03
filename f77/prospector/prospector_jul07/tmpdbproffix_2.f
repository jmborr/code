c 	uses TMalign.v2c.f
c	*************************************************************************
*     			TM-align version 2.0
*
c*	******************************************************************
*     Previous Bug fixing by Yang.
*     TM-align version 0.1. A small bug of two-point superposition was 
*     fixed at June 1, 2005, compared with previous version.
*
*     This program is to identify the best alignment of two protein 
*     structures to give the best TM-score. By default, TM-score is 
*     normalized by the second protein. The program can be freely 
*     copied or modified or redistributed.
*
*     Reference:
*     Yang Zhang, Jeffrey Skolnick, Nucl. Acid Res. 2005 33: 2303-9
*     (For comments, please email to: zhang6@buffalo.edu)
*cccc******************************************************************
* 
*     The new fixed version of TM-align V.2. has been modified extensively to 
*     identify alignment with highest TM-score between given two structures. 
*     The changes are:
*
*     1. The gapless threading is modified with changes in step size and
*        the way of scanning.
*     2. In gapless threading, we have used two scores, GL and dpout, to 
*	 identify best possible alignment. We infact take first 40 or 80
* 	 alignments given by gapless threading and process them to identify
*	 the best possible alignment using RMSD scoring and TMscore matrices.
*     3. We need atleast 3 residues for superposition in the RMSD calculation.
*     4. With the above mentioned modifications, the program runs relatively
*        slow in comparison to initial TMalign. However, it ensures alignment 
*        with better TMscore. The next step would involve increasing the speed of
*	 program if possible.
*     5. The lower threshold value d0_search is changes to 3.0 instead of 4.5!!
*
*     6. The fast version of TMalign is first run and if we fail to get TMscore
*        > 0.85, we switch to slow version.
*
*  

*     July, 2007. Shashi B. Pandit. Report bugs to spandit3@mail.gatech.edu
*******************************************************************************
c	uses tm cutoff of 0.65 for profile alignments

*     Reference:
*     Yang Zhang, Jeffrey Skolnick, Nucl. Acid Res. 2005 33: 2303-9
*     (For comments, please email to: zhang6@buffalo.edu)
*************************************************************************
      
      program compares
      PARAMETER(NMAX=1000)
      PARAMETER(NCASES=8700)

      logical*1 homology(ncases)

      COMMON/BACKBONE/XA(3,NMAX,0:1)
      common/seq/nresseq,nseq_ori

      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002

      character*4 genome
      character*5 namep,NAME(0:NCASES),name2
      character*6 prodin


      dimension jseq(0:nmax),icode(0:nmax)
      dimension m1(nmax),m2(nmax)
      dimension mm1(nmax),mm2(nmax)
      dimension mm1a(nmax,ncases),mm2a(nmax,ncases)
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      dimension nres(ncases),atm(ncases),ntm(ncases),trms(ncases),item(ncases)
      dimension x(3,nmax,ncases)
c	===============================
      common/init/invmap_i(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d8/d8

      common/dpcommon/dpout,ivalue !value in TMscore search will define what 
                                   !search we are doing, a quick or exhaustive
                                   ! 0 = quick search and 1 = exhaustive search


ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   
	 character*255 templatedat
	 logical*1 tncbi
	 dimension ix(0:22)
	 character*1 ancbi(1:22),aa1(0:19)
	 dimension tncbi(0:22),ncbi(22)
	 dimension xm(0:59,0:nmax)
	 dimension xmav(0:59,nmax),anhist(nmax)

        DATA AA1/ 'G','A','S','C','V','T','I','P',
     &               'M','D','N','L',
     &               'K','E','Q','R',
     &               'H','F','Y','W'/
        DATA ancbi/'A','B','C','D','E','F','G','H','I','K','L','M','N','P',
     &  'Q','R','S','T','V','W','X','Y'/

c      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
c    &     'PRO','MET','ASP','ASN','LEU',
c    &     'LYS','GLU','GLN','ARG',
c    &     'HIS','PHE','TYR','TRP','CYX'/
c     character*1 slc(-1:20)
c      data slc/'X','G','A','S','C','V','T','I',
c    &     'P','M','D','N','L','K','E','Q','R',
c    &     'H','F','Y','W','C'/

		 do jres=1,22
		 tncbi(jres)=.false.
		 do ires=0,19
		 if(aa1(ires).eq.ancbi(jres))then
		 tncbi(jres)=.true.
			  ncbi(jres)=ires
			  go to 22
		 end if
		 end do
22      continue
		 end do
		Open(unit=11,file='namesize')
		read(11,*)nsize
		close(11)
        read(5,13)prodin
13      format(a6)
        close(11)

c	&&&&&&& templatedir:    &&&&&&&&&&
c		open(unit=1,file='templatedir')
c		read(unit=1,fmt='A255')templatedat
c		nadr=0
c		do i=255,1,-1
c		if(templatedat(i:i) .eq. ' ')nadr=i
c		end do
c		nadr=nadr-1
c	&&&&&&& templatedir:    &&&&&&&&&&

        OPEN(unit=33, file='LIST.jul2007good35')
	read(33,*)ndata
	do 1100 ik=1,ndata
	read(33,33)name(IK)
33	format(a5)
1100 	CONTINUE
	DO iK=1,ndata
	OPEN(unit=30,file='/local/images/templates-2007081400/orien6/'//name(ik)//'.SEQ')
	rewind(30)
	read(30,2)name2,nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
2	format(a5,1x,i5)
	close(30)

	END DO
c	=================================
        OPEN(UNIT=9,FILE='stat_tmpdbproffix_2'//prodin)
	do ik=1,ndata
	open(unit=45,file='/local/images/templates-2007081400/cafiles6/'//name(ik)//'.real')
        read(45,*)nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
	do i=1,nres(ik)
	read(45,*)id,(x(j,i,ik),j=1,3)
	end do
	END DO

        OPEN(UNIT=19,FILE='av_'//prodin)
        OPEN(unit=33,file='LIST.tm'//prodin)
        read(33,*)ndata2,genome
	do 1001 ik=1,ndata2
	read(33,33)namep
        OPEN(unit=45,file='/local/images/templates-2007081400/cafiles6/'//NAMEp//'.real')
	rewind(45)
	read(45,*)MRES
	open(unit=11,file='../test/'//namep)
	write(11,*)namep,mres
	close(11)
	if(MRES.gt.nmax)mres=nmax
        L2=MRES
        do i=1,L2
        read(45,*)mm1(i),xa(1,i,0),xa(2,i,0),xa(3,i,0)
        end do
        close(45)
	OPEN(unit=30,file='/local/images/templates-2007081400//orien6/'//NAMEp//'.SEQ')
	rewind(30)
	read(30,2)name2,mres
	if(mres.gt.nmax)mres=nmax
	read(30,*)(jseq(i),i=1,mres)
	close(30)
	do i=1,mres
	if(jseq(i).gt.19)jseq(i)=0
	end do

	! read in the PDB library/
	nseq1=mres
c				==============================
c   initialization of structure based average array
		do ires=0,59
			do i=1,nseq1
			xmav(ires,i)=0.
			end do
		enddo
		do i=1,nseq1
			anhist(i)=0.
			end do
c				==============================

c       do jk=1,ndata
c       homology(jk)=.true.
c       end do


	do 1000 jk=1,ndata
	atm(jk)=-3000.
	ntm(jk)=0
	trms(jk)=0.
	item(jk)=jk
c	if(homology(jk))then
	do j=1,nres(jk)
	mm2(j)=j
	do l=1,3
	xa(l,j,1)=x(l,j,jk)
	end do
	end do
	nseq2=nres(jk)

*!!!  Scale of TM-score in search is based on the smaller protein --------->
      d0_min=0.5
c      if(m_d0_min.eq.1)then
c         d0_min=d0_min_input
c      endif
      anseq_min=min(nseq1,nseq2)
      anseq=anseq_min           !length for defining TMscore in search
      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
c      if(d00.lt.4.5)d00=4.5 updated
      if(d00.lt.3.0)d00=3.0

!!! change to d00=3.5
      d002=d00**2
      nseq=max(nseq1,nseq2)
      do i=1,nseq
         n1(i)=i
         n2(i)=i
      enddo
      
***** do alignment **************************
      CALL super_align          !to find invmap(j)
      
************************************************************
***   resuperpose to find residues of dis<d8 ------------------------>
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
            m1(n_al)=i          !for recording residue order
            m2(n_al)=j
         endif
      enddo
      d0_input=d0
      call TMscore8(d0_input,n_al,xtm1,ytm1,ztm1,n1,n_al,
     &     xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !TM-score with dis<d8 only
	
*!!!  Output TM-score is based on the second protein------------------>
c      anseq=nseq2               !length for defining final TMscore
      anseq=nseq1               !length for defining final TMscore
c      if(m_fix.eq.1)anseq=L_fix !input length
c     if(m_ave.eq.1)anseq=(nseq1+nseq2)/2.0 !<L>
      if(anseq.lt.anseq_min)anseq=anseq_min
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      
***   remove dis>d8 in normal TM-score calculation for final report----->
      j=0
      n_eq=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)
            m1(j)=m1(i)
            m2(j)=m2(i)
c            if(ss1(m1(i)).eq.ss2(m2(i)))then
c              n_eq=n_eq+1
c           endif
         endif
      enddo
c      seq_id=float(n_eq)/(n_al+0.00000001)
      n8_al=j
      d0_input=d0
      call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &     xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm) !normal TMscore
      rmsd8_al=Rcomm
c      TM8=TM8*n8_al/anseq       !TM-score after cutoff
! normalized by the length of the query sequence nseq1
        TM8=TM8*float(n8_al)/nseq1       !TM-score after cutoff
      	atm(jk)=tm8
	ntm(jk)=LCOMM
	if(LCOMM .ne. n8_al)then
	write(9,*)' *** LCOMM N8_al are unequal', lcomm,n8_al
	ntm(jk)=n8_al
	end if
	trms(jk)=RCOMM	

	do i=1,n8_al
! 	output format: template, target
	mm2a(i,jk)=mm2(m2(i))
	Mm1a(i,jk)=mm1(m1(i))
	end do
c	END IF
1000	continue
        do jk=1,ndata
        do kk=jk+1,ndata
        if(atm(jk).lt.atm(kk))then
        etem=atm(jk)
        atm(jk)=atm(kk)
        atm(kk)=etem

        etem=trms(jk)
        trms(jk)=trms(kk)
        trms(kk)=etem

	itt=item(jk)
        item(jk)=item(kk)
        item(kk)=itt

        itt=ntm(jk)
        ntm(jk)=ntm(kk)
        ntm(kk)=itt
        end if
        end do
        end do
	av=0.
	av2=0.
	nq=0
        do jk=1,ndata
	jjk=item(jk)
c	if(homology(jjk))then		
	nq=nq+1
	av=av+atm(jk)
	av2=av2+atm(jk)**2
c	end if
	end do
	av=av/nq
	av2=(av2/nq)
	std=sqrt(av2-av*av)
	write(9,99)'number of cases= ',nq, ' <tm> = ',av,' std =',std
99	format(a17,1x,i4,a8,1x,1f7.3,1x,a6,1f7.3)
	write(19,19)nq,av,std
19	format(i4,1x,2(1f7.3,1x))
	nth2=0
!	here nth is the number of templates that have a tm score >0.
	do jk=1,ndata
	jjk=item(jk)
c	if(homology(jjk))then
	if(atm(jk).gt.0.65)then
	OPEN(unit=30,file='/local/images/templates-2007081400/orien6/'//name(jjk)//'.SEQ')
	rewind(30)
	read(30,*)
	read(30,*)(icode(i),i=1,nres(jjk))
	close(30)

        open(unit=59,file='/local/images/templates-2007081400/templatesprofilejul07/'//name(jjk)//'.profile')
        read(59,*)nseq11
        if(nseq1.gt.nmax)nseq11=nmax
        do i=1,nseq11
        read(59,59)ii,(xm(ires,i),ires=0,59)
c59      format(i4,1x,20(1f7.3,1x),/,2(5x,20(1f7.3,1x),/))
        end do
        close(59)


	do i=1,nres(jjk)
	if(icode(i).gt.19)icode(i)=0.
	end do
	nth2=nth2+1
	pid=0.
	do ii=1,ntm(jk)
	j=mm2a(ii,jjk)
	i=mm1a(ii,jjk)
	anhist(i)=anhist(i)+1
	do ires=0,59
	xmav(ires,i)=xmav(ires,i)+xm(ires,j)
	end do
	if(jseq(i).eq.icode(j))pid=pid+1
	end do
	pid=pid/nseq1
	write(9,93)name(jjk),pid,atm(jk),ntm(jk),trms(jk)	
	end if
c	end if
	end do

93	format(a5,1x,2(1f9.3),1x,i4,1x,2(1f7.3,1x))
! this is fixed
	open(unit=49,file='../struhistjul07/'//namep//'.hist2')
	NTH=0
	do jk=1,ndata
	if(atm(jk).gt.0.40)NTH=NTH+1
	end do
	write(49,*)NTH,av,std
	do jk=1,ndata
	if(atm(jk).gt.0.40)THEN
	jjk=item(jk)
	open(unit=30,file='/local/images/templates-2007081400/orien6/'//name(jjk)//'.SEQ')
	rewind(30)
	read(30,*)
	read(30,*)(icode(i),i=1,nres(jjk))
	pid=0.
	do ii=1,ntm(jk)
	j=mm2a(ii,jjk)
	i=mm1a(ii,jjk)
	if(jseq(i).eq.icode(j))pid=pid+1
	end do
	pid=pid/nseq1
	write(49,93)name(jjk),pid,atm(jk),ntm(jk),trms(jk)	
	end if
	end do
	close(49)
	write(9,*)namep,nth
	do i=1,nseq1
		do ires=0,59
		xmav(ires,i)=xmav(ires,i)/anhist(i)
		end do
		end do
c	******* CALCULATE THE AVERAGE PROFILE *******
	open(unit=59,file='../struprofilejul07/'//namep//'.profile2')
	write(59,*)nseq1,nth2,pid
	do i=1,nseq1
	write(59,59)i,(xmav(ires,i),ires=0,59)
59	format(i4,1x,20(1f7.3,1x),/,2(5x,20(1f7.3,1x),/))
	end do
	close(59)
1001	continue
9999 	END
c       *********************************************************!
***********************************************************************
***********************************************************************
*     Structure superposition
***********************************************************************
***********************************************************************
***********************************************************************
      SUBROUTINE super_align
      PARAMETER(nmax=1000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      dimension gapp(100)
      common/dpcommon/dpout,ivalue
      common/sec/isec(nmax),jsec(nmax)

      TMmax=0
      n_gapp=2
      gapp(1)=-0.6
      gapp(2)=0

c      n_gapp=11
c      do i=1,n_gapp
c         gapp(i)=-(n_gapp-i)
c      enddo

*11111111111111111111111111111111111111111111111111111111
*     get initial alignment from gapless threading
**********************************************************
      call get_initial          !gapless threading
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, matrix score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

*****************************************************************
*       initerative alignment, for different gap_open:
*****************************************************************
      DO 111 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 222 id=1,30         !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 33
            endif
            TM_old=TM
 222     continue
 33      continue
 111  continue

*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*     get initial alignment from secondary structure alignment
**********************************************************
      call get_initial2         !DP for secondary structure
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 1111 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 2222 id=1,30	!maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
c     write(*,21)gap_open,rmsd_al,n_al,TM
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 333
            endif
            TM_old=TM
 2222    continue
 333     continue
 1111 continue
      
*333333333333333333333333333333333333333333333333333333333333
*     get initial alignment from invmap0+SS
*************************************************************
      call get_initial3         !invmap0+SS
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 1110 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 2220 id=1,30	!maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
c     write(*,21)gap_open,rmsd_al,n_al,TM
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 330
            endif
            TM_old=TM
 2220    continue
 330     continue
 1110 continue


       if(TMmax.lt.0.85) then

        call get_initial2
        call get_initial0(aTM)
           
        if(aTM.gt.TMmax) then
         do j=1,nseq2
          invmap0(j)=invmap_i(j)
         enddo
        endif

       endif


c^^^^^^^^^^^^^^^ best alignment invmap0(j) found ^^^^^^^^^^^^^^^^^^
      RETURN
      END

**************************************************************
*     get initial alignment invmap0(i) from gapless threading
**************************************************************
      subroutine get_initial
      PARAMETER(nmax=1000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)

      aL=min(nseq1,nseq2)
      idel=aL/2.5               !minimum size of considered fragment
      if(idel.le.5)idel=5
      n1=-nseq2+idel
      n2=nseq1-idel
      GL_max=0
      do ishift=n1,n2
         L=0
         do j=1,nseq2
            i=j+ishift
            if(i.ge.1.and.i.le.nseq1)then
               L=L+1
               invmap(j)=i
            else
               invmap(j)=-1
            endif
         enddo
         if(L.ge.idel)then
            call get_GL(GL)
            if(GL.gt.GL_max)then
               GL_max=GL
               do i=1,nseq2
                  invmap_i(i)=invmap(i)
               enddo
            endif
         endif
      enddo

      return
      end

**************************************************************
*     get initial alignment invmap0(i) from secondary structure
**************************************************************
      subroutine get_initial2
      PARAMETER(nmax=1000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)

********** assign secondary structures ***************
c     1->coil, 2->helix, 3->turn, 4->strand
      do i=1,nseq1
         isec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq1)then
            dis13=diszy(0,j1,j3)
            dis14=diszy(0,j1,j4)
            dis15=diszy(0,j1,j5)
            dis24=diszy(0,j2,j4)
            dis25=diszy(0,j2,j5)
            dis35=diszy(0,j3,j5)
            isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
      do i=1,nseq2
         jsec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq2)then
            dis13=diszy(1,j1,j3)
            dis14=diszy(1,j1,j4)
            dis15=diszy(1,j1,j5)
            dis24=diszy(1,j2,j4)
            dis25=diszy(1,j2,j5)
            dis35=diszy(1,j3,j5)
            jsec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
      call smooth               !smooth the assignment

********** score matrix **************************
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=1
            else
               score(i,j)=0
            endif
         enddo
      enddo

********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(NSEQ1,NSEQ2)      !produce alignment invmap(j)
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

**************************************************************
*     get initial alignment invmap0(i) from secondary structure 
*     and previous alignments
**************************************************************
      subroutine get_initial3
      PARAMETER(nmax=1000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)

********** score matrix **************************
      do i=1,nseq2
         invmap(i)=invmap0(i)
      enddo
      call get_score1           !get score(i,j) using RMSD martix
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=0.5+score(i,j)
            else
               score(i,j)=score(i,j)
            endif
         enddo
      enddo

********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(NSEQ1,NSEQ2)      !produce alignment invmap(j)
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

**************************************************************
*     smooth the secondary structure assignment
**************************************************************
      subroutine smooth
      PARAMETER(nmax=1000)
      common/sec/isec(nmax),jsec(nmax)
      common/length/nseq1,nseq2

***   smooth single -------------->
***   --x-- => -----
      do i=3,nseq1
         if(isec(i).eq.2.or.isec(i).eq.4)then
            j=isec(i)
            if(isec(i-2).ne.j)then
               if(isec(i-1).ne.j)then
                  if(isec(i+1).ne.j)then
                     if(isec(i+1).ne.j)then
                        isec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=3,nseq2
         if(jsec(i).eq.2.or.jsec(i).eq.4)then
            j=jsec(i)
            if(jsec(i-2).ne.j)then
               if(jsec(i-1).ne.j)then
                  if(jsec(i+1).ne.j)then
                     if(jsec(i+1).ne.j)then
                        jsec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo

***   smooth double -------------->
***   --xx-- => ------
      do i=1,nseq1-5
         if(isec(i).ne.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
         if(isec(i+3).eq.2)then
         if(isec(i+4).ne.2)then
         if(isec(i+5).ne.2)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(isec(i).ne.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
         if(isec(i+3).eq.4)then
         if(isec(i+4).ne.4)then
         if(isec(i+5).ne.4)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo
      do i=1,nseq2-5
         if(jsec(i).ne.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
         if(jsec(i+3).eq.2)then
         if(jsec(i+4).ne.2)then
         if(jsec(i+5).ne.2)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(jsec(i).ne.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
         if(jsec(i+3).eq.4)then
         if(jsec(i+4).ne.4)then
         if(jsec(i+5).ne.4)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo

***   connect -------------->
***   x-x => xxx
      do i=1,nseq1-2
         if(isec(i).eq.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
            isec(i+1)=2
         endif
         endif
         endif

         if(isec(i).eq.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
            isec(i+1)=4
         endif
         endif
         endif
      enddo
      do i=1,nseq2-2
         if(jsec(i).eq.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
            jsec(i+1)=2
         endif
         endif
         endif

         if(jsec(i).eq.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
            jsec(i+1)=4
         endif
         endif
         endif
      enddo

      return
      end

**************************************************************
*      New function to get best possible alignment
* Added by: Shashi
**************************************************************
      subroutine get_initial0(aTM)
      PARAMETER(nmax=1000)
      parameter(top=100)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      common/dpcommon/dpout,ivalue
      common/kela/iali(nmax,top),tmal(top),iswap(top),dpo(top)
     &            ,jali(nmax,top) ! for alignment
      common/sec/isec(nmax),jsec(nmax)
      dimension invmapA(nmax),invmapB(nmax) ! storing alignment
      dimension score1(nmax,nmax),score2(nmax,nmax)
      common/lcdim/iiali(nmax,top),jjali(nmax,top),gap(10)
      dimension map(nmax),mapa(nmax)
      logical check

      aL=min(nseq1,nseq2)
      aL1=max(nseq1,nseq2)
      fr=aL/aL1

      if(aL.lt.170) then
        idel=aL/2.5
        minoverlap=max(12,idel)
        if(fr.lt.0.16)then
          istep=3
        elseif(fr.ge.0.16.and.fr.lt.0.25) then
         istep=(minoverlap/3)-1
        elseif(fr.ge.0.25.and.fr.lt.0.30) then
         istep=(minoverlap/3)-3
        elseif(fr.ge.0.30.and.fr.lt.0.40) then
         istep=(minoverlap/3)-2
        elseif(fr.ge.0.40.and.fr.lt.0.50)then
         istep=minoverlap/3
        else
         if(aL1.lt.100)then
          istep=1
         else
          istep=3
         endif
        endif

        n1=minoverlap-nseq1
        n2=nseq2-minoverlap
       if(istep.eq.0)istep=1

       if(nseq2.lt.nseq1) then
         kk0=n2
         iadd=0
         do while (mod(kk0,istep).ne.0)
          iadd=iadd+1
          kk0=n2+iadd
         enddo
        inum=(nseq1-minoverlap-iadd)
        imult=inum/istep

        inadd=inum-(imult*istep)
        n1=n1+inadd
       endif

       ntempM=60
       ntempD=60

      else

       if(fr.lt.0.15)then
       idel=aL/2.5
       istep=4
        minoverlap=max(12,idel)
       else
        idel=(aL/2.5)
        minoverlap=max(12,idel)
        istep=3
       endif

        n1=minoverlap-nseq1
        n2=nseq2-minoverlap

       if(nseq2.lt.nseq1) then
         kk0=n2
         iadd=0
         do while (mod(kk0,istep).ne.0)
          iadd=iadd+1
          kk0=n2+iadd
         enddo
        inum=(nseq1-minoverlap-iadd)
        imult=inum/istep

        inadd=inum-(imult*istep)
        n1=n1+inadd
       endif

       ntempM=20
       ntempD=20

      endif
      ntot=0    ! total no. of scans
      itemp=0    ! no. of alignment to be considered
      itemp1=0   ! no. of alignment for dpout

      gap_open=0.0

      do ishift=n1,n2,istep
          LL=0
         do j=1,nseq2
            i=j-ishift
            if(i.ge.1.and.i.le.nseq1)then
               invmap(j)=i
               mapa(j)=i
               LL=LL+1
            else
               invmap(j)=-1
               mapa(j)=-1
            endif
         enddo

        call get_GL(GL)

            ntot=ntot+1

             if(ntot.lt.(ntempM+1)) then
              tmal(ntot)=GL
               do j=1,nseq2
                 iali(j,ntot)=invmap(j)
               enddo

                  if(ntot.eq.ntempM) then
                   itemp=ntot
                   call sort(itemp,1)
                     do i=1,itemp
                      do j=1,nseq2
                       iiali(j,i)=iali(j,iswap(i))
                      enddo
                     enddo
                  endif
             else
               do i=itemp,1,-1
                if(GL.gt.tmal(i))then
                  do ii=1,i-1
                    do j=1,nseq2
                     iiali(j,ii)=iiali(j,ii+1)
                    enddo
                     tmal(ii)=tmal(ii+1)
                  enddo
                    do j=1,nseq2
                     iiali(j,i)=invmap(j)
                    enddo
                     tmal(i)=GL
                  goto 11
                endif
               enddo
  11         continue
             endif

           call get_score1
           call dp(nseq1,nseq2)
            do j=1,nseq2
             invmap(j)=mapa(j)
            enddo

             if(ntot.lt.(ntempD+1)) then
                dpo(ntot)=dpout
               do j=1,nseq2
                 jali(j,ntot)=invmap(j)
               enddo
                  if(ntot.eq.ntempD) then
                   itemp1=ntot
                   call sort(itemp1,2)
                     do i=1,itemp1
                      do j=1,nseq2
                       jjali(j,i)=jali(j,iswap(i))
                      enddo
                     enddo
                  endif
             else
               do i=itemp1,1,-1
                if(dpout.gt.dpo(i))then
                  do ii=1,i-1
                    do j=1,nseq2
                     jjali(j,ii)=jjali(j,ii+1)
                    enddo
                     dpo(ii)=dpo(ii+1)
                  enddo
                    do j=1,nseq2
                     jjali(j,i)=invmap(j)
                    enddo
                     dpo(i)=dpout
                  goto 12
                endif
               enddo
  12         continue
             endif
      enddo

        if(ntot.lt.ntempD) then
          itemp=ntot
          call sort(itemp,1)
             do i=1,itemp
              do j=1,nseq2
               iiali(j,i)=iali(j,iswap(i))
              enddo
             enddo
           itemp1=ntot
           call sort(itemp1,2)
             do i=1,itemp1
              do j=1,nseq2
               jjali(j,i)=jali(j,iswap(i))
              enddo
             enddo
        endif

cccc  Finding the best alignment from top 40/80 alignments using dpout function.
cccc  We use a brute force method and scan all alignment with combination of 
cccc  RMSD and TMscore matrices with different gap opening penalties!!
cccc  

      ivalue=1
      b1TMmax=0.0
      ist=itemp1-5 !scanning only best 5 alignments

      do 212 ial=ist,itemp1
       do kkk=1,2
        if(kkk.eq.1) then
         gap_open=-0.6
        else
         gap_open=0.0
        endif
        do j=1,nseq2
         invmap(j)=jjali(j,ial)
        enddo

         call get_score
         call dp(nseq1,nseq2)

         if(check()) goto 211
         
         call iter(4,b1TMmax,2)
          do j=1,nseq2
           map(j)=invmap_i(j)
          enddo
       enddo
  211 continue
  212 continue

      do 214 ial=1,itemp1
        do j=1,nseq2
          invmap(j)=jjali(j,ial)
        enddo
        call get_score1
        call iter(6,b1TMmax,1)
         do j=1,nseq2
          map(j)=invmap_i(j)
         enddo

         do kkk=1,2
          do j=1,nseq2
            invmap(j)=jjali(j,ial)
          enddo
          if(kkk.eq.1) then
           gap_open=-0.6
          else
           gap_open=0.0
          endif
           call get_score
           call dp(nseq1,nseq2)
           if(check())goto 215
           call get_score1
           call iter(6,b1TMmax,1)
            do j=1,nseq2
             map(j)=invmap_i(j)
            enddo
         enddo
  215 continue
  214 continue


cccc Scanning of best alignment from dpout score is done. Next we scan alignments from
cccc score GL. We use a combination of RMSD/TMscore matrices to find best alignment.
cccc There is finer optimization of scoring matrices to incorporate the secondary structure
cccc in process of finding highest TMscore alignment.

       nid=20
       b2TMmax=0.0
       ivalue=1
       n_gapp=2
       gap(1)=-0.6
       gap(2)=0.0

       do ial=1,itemp
        do j=1,nseq2
         invmap(j)=iiali(j,ial)
        enddo
        call get_score

           do 112 i_gapp=1,n_gapp
             gap_open=gap(i_gapp)

                 if(fr.gt.0.25)then
                    if(TM.lt.0.20.or.TM.gt.0.50) then
                        do i=1,nseq1
                         do j=1,nseq2
                          if(isec(i).eq.jsec(j))then
                             score(i,j)=0.3+score(i,j)
                          else
                            score(i,j)=score(i,j)
                          endif
                        enddo
                       enddo
                    else
                       do i=1,nseq1
                         do j=1,nseq2
                          if(isec(i).eq.jsec(j))then
                            score(i,j)=0.4+score(i,j)
                          else
                            score(i,j)=score(i,j)
                          endif
                         enddo
                       enddo
                   endif
                 endif

                    if(fr.le.0.25)then
                     if(TM.lt.0.20.or.TM.gt.0.50)then
                       do i=1,nseq1
                         do j=1,nseq2
                          if(isec(i).eq.jsec(j))then
                            score(i,j)=0.2+score(i,j)
                          else
                           score(i,j)=score(i,j)
                         endif
                         enddo
                       enddo
                      else
                       do i=1,nseq1
                         do j=1,nseq2
                         if(isec(i).eq.jsec(j))then
                            score(i,j)=0.4+score(i,j)
                          else
                            score(i,j)=score(i,j)
                          endif
                         enddo
                       enddo
                     endif
                   endif

               do 114 id=1,nid
                  call dp(nseq1,nseq2)
                  if(check()) goto 115
                  call get_score
                  if (TM.gt.b2TMmax) then
                     b2TMmax=TM
                      do j=1,nseq2
                       invmapA(j)=invmap(j)
                      enddo
                  endif

                    if(id.gt.1)then
                      diff=abs(TM-TM_old)
                      if(diff.lt.0.000001)goto 115
                    endif
                 TM_old=TM
  114          continue
  115      continue
  112     continue
       enddo     

cccc Extension of first alignment from dpout using iter routine.

        do j=1,nseq2
          invmap(j)=map(j)
          invmap_i(j)=invmap(j)
        enddo

        call iter(20,b1TMmax,0)

         do j=1,nseq2
          map(j)=invmap_i(j)
         enddo

cccccccccccccc The Alignment with best GL is taken for extension and iteration

       ivalue=1     ! Quick search

cccc THE invmap has BEST GL alignment stored now. This point onward invmapB will have
cccc the best alignment resulting from progression of best GL alignment.
cccc Mixing of SS with RMSD rotation matrix by any two alignments.
cccc the invmap(j) can be for best GL OR from best alignment 
cccc from previous iteration. The first alignment is from best GL
cccc and second alignment is from best of all alignment

**** CHANGED Jul 25, invmapB is not longer used. we have only two alignments best of GL
**** and best of DPOUT

       iRUN=2       
       if(b1TMmax.eq.b2TMmax) then
        iRUN=1
       endif

       do 217 istep=1,iRUN ! 1-best alignment, 2- best DP alignment

         if(istep.eq.1) then
          do j=1,nseq2
            invmap(j)=invmapA(j)
            invmap_i(j)=invmap(j)
          enddo
         elseif(istep.eq.2) then
          do j=1,nseq2
            invmap(j)=map(j)
            invmap_i(j)=invmap(j)
          enddo
         endif

cc SS based alignment
           do i=1,nseq1
            do j=1,nseq2
              if(isec(i).eq.jsec(j))then
                score(i,j)=1
              else
                score(i,j)=0
              endif
            enddo
           enddo

           gap_open=-1.0
           call dp(nseq1,nseq2)
           if(check()) goto 216
           call get_score

             if(istep.eq.1)then
                if(TM.gt.b2TMmax) then
                 b2TMmax=TM
                 do jj=1,nseq2
                  invmapA(jj)=invmap(jj)
                 enddo
                endif
               call iter(20,b2TMmax,1)
                do jj=1,nseq2
                  invmapA(jj)=invmap_i(jj)
                  invmap(jj)=invmap_i(jj)
                enddo
             elseif(istep.eq.2)then
                if(TM.gt.b1TMmax) then
                 b1TMmax=TM
                 do jj=1,nseq2
                  map(jj)=invmap(jj)
                 enddo
                endif
               call iter(20,b1TMmax,1)
                do jj=1,nseq2
                  map(jj)=invmap_i(jj)
                  invmap(jj)=invmap_i(jj)
                enddo
             endif

cccc Increasing score matrix using SS alignment 

             call get_score1
               do i=1,nseq1
                do j=1,nseq2
                 if(isec(i).eq.jsec(j))then
                   score(i,j)=0.5+score(i,j)
                 else
                   score(i,j)=score(i,j)
                 endif
                enddo
               enddo

                gap_open=-1.0
                call dp(nseq1,nseq2)
                if(check()) goto 216
                call get_score

                  if(istep.eq.1)then
                     if(TM.gt.b2TMmax) then
                      b2TMmax=TM
                        do jj=1,nseq2
                         invmapA(jj)=invmap(jj)
                        enddo
                     endif
                       call iter(20,b2TMmax,1)
                        do jj=1,nseq2
                         invmapA(jj)=invmap_i(jj)
                         invmap(jj)=invmap_i(jj)
                         enddo
                  elseif(istep.eq.2)then
                     if(TM.gt.b1TMmax) then
                      b1TMmax=TM
                        do jj=1,nseq2
                         map(jj)=invmap(jj)
                        enddo
                     endif
                       call iter(20,b1TMmax,1)
                         do jj=1,nseq2
                           map(jj)=invmap_i(jj)
                           invmap(jj)=invmap_i(jj)
                         enddo
                  endif

  216  continue
  217  continue

       imax=0 ! if imax =1 -> invmapA, imax=2->map
       bmax=0
       if(b2TMmax.gt.bmax)then
        imax=1
        bmax=b2TMmax
       endif
       if(b1TMmax.gt.bmax)then
        imax=2
        bmax=b1TMmax
       endif

       if(imax.eq.1) then
        do j=1,nseq2
         invmap_i(j)=invmapA(j)
        enddo
       elseif(imax.eq.2) then
        do j=1,nseq2
         invmap_i(j)=map(j)
        enddo
       endif
       aTM=bmax

      return
      end

cccc Check routine to see if the invmap(i) has sufficient number of 
cccc residues for the next step.

      logical function check
      parameter(nmax=1000)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/length/nseq1,nseq2
      
      check=.false.
      LL=0
      do jj=1,nseq2
       if(invmap(jj).ge.1)LL=LL+1
      enddo

      if(LL.lt.4)check=.true.

      return
      end

cccc  Routine to do iteration of the alignment given the  
cccc  no. of iteration (nid), to find highest TMscore (aTMmax)
cccc  and which score to use (itt). itt=0, score matrix need to 
cccc  be filled, itt=1, score matrix is filled already and 
cccc  itt=2, fill score matrix in special way!! (optimized condition)

      subroutine iter(nid,aTMmax,itt)
      parameter (nmax=1000)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/init/invmap_i(nmax)
      common/TM/TM,TMmax
      common/dpcommon/dpout,ivalue
      dimension gap(10)
      logical check

      ivalue=1

      if(check()) return 

      if(itt.eq.0) then
       call get_score
        if (TM.gt.aTMmax) then
         aTMmax=TM
          do j=1,nseq2
           invmap_i(j)=invmap(j)
          enddo
        endif
      endif
       
      n_gapp=2
      gap(1)=-0.6
      gap(2)=0.0

      do 1111 i_gapp=1,n_gapp
        gap_open=gap(i_gapp)
      
       if(itt.eq.2) then
        call get_score1
         call dp(nseq1,nseq2)
          if(check()) goto 333
          call get_score
       endif
         do 2222 id=1,nid
           call dp(nseq1,nseq2)
           if(check()) goto 333
           call get_score
             if(TM.gt.aTMmax)then
               aTMmax=TM
               do j=1,nseq2
                 invmap_i(j)=invmap(j)
               enddo
             endif
                    if(id.gt.1)then
                      diff=abs(TM-TM_old)
                      if(diff.lt.0.000001)goto 333
                    endif
           TM_old=TM
 2222     continue
  333   continue
 1111  continue
 
 1112 continue
      return
      end

cccc Sorting the alignment based on dpout or GL scores
cccc temp is the number of alignment and kind=1, for GL and
cccc kind=2 for dpout scores.

      subroutine sort(nos,kind)
      parameter(top=100)
      parameter(nmax=1000)
      common/length/nseq1,nseq2
      common/kela/iali(nmax,top),tmal(top),iswap(top),dpo(top)
     &            ,jali(nmax,top) ! for alignment

       do i=1,nos
        iswap(i)=i
       enddo

      if(kind.eq.1)then
       do i=1,nos
         do j=1,nos-i
          if(tmal(j+1).lt.tmal(j))then
           tmp=tmal(j)
           tmal(j)=tmal(j+1)
           tmal(j+1)=tmp
           ii=iswap(j)
           iswap(j)=iswap(j+1)
           iswap(j+1)=ii
          endif
         enddo
       enddo
      elseif(kind.eq.2)then
       do i=1,nos
         do j=1,nos-i
          if(dpo(j+1).lt.dpo(j))then
           tmp=dpo(j)
           dpo(j)=dpo(j+1)
           dpo(j+1)=tmp
           ii=iswap(j)
           iswap(j)=iswap(j+1)
           iswap(j+1)=ii
          endif
         enddo
       enddo

      endif

      return
      end
cccc

*************************************************************
*     assign secondary structure:
*************************************************************
      function diszy(i,i1,i2)
      PARAMETER(nmax=1000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     &     +(xa(2,i1,i)-xa(2,i2,i))**2
     &     +(xa(3,i1,i)-xa(3,i2,i))**2)
      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end

****************************************************************
*     quickly calculate TM-score with given invmap(i) in 3 iterations
****************************************************************
      subroutine get_GL(GL)
      PARAMETER(nmax=1000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d00/d00,d002

      dimension xo1(nmax),yo1(nmax),zo1(nmax)
      dimension xo2(nmax),yo2(nmax),zo2(nmax)
      dimension dis2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
            xo1(n_al)=xa(1,i,0)
            yo1(n_al)=xa(2,i,0)
            zo1(n_al)=xa(3,i,0)
            xo2(n_al)=xa(1,j,1)
            yo2(n_al)=xa(2,j,1)
            zo2(n_al)=xa(3,j,1)
         endif
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      GL=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         GL=GL+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002
 21   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3)then
         d002t=d002t+.5
         goto 21
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G2=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G2=G2+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002+1
 22   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3)then
         d002t=d002t+.5
         goto 22
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G3=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G3=G3+1/(1+dis2(i)/(d0**2))
      enddo
      if(G2.gt.GL)GL=G2
      if(G3.gt.GL)GL=G3

c^^^^^^^^^^^^^^^^ GL done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

****************************************************************
*     with invmap(i) calculate TM-score and martix score(i,j) for rotation 
****************************************************************
      subroutine get_score
      PARAMETER(nmax=1000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/dpcommon/dpout,ivalue

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for TM-score:
            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
         endif
      enddo
***   calculate TM-score for the given alignment----------->
      d0_input=d0
      call TMscore8_search(d0_input,n_al,xtm1,ytm1,ztm1,n1,
     &     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !simplified search engine
      TM=TM*n_al/anseq          !TM-score
***   calculate score matrix score(i,j)------------------>
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d0**2)
         enddo
      enddo
      
c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

****************************************************************
*     with invmap(i) calculate score(i,j) using RMSD rotation 
****************************************************************
      subroutine get_score1
      PARAMETER(nmax=1000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      common/d0min/d0_min
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo

***   calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d02)
         enddo
      enddo

c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
*** dis<8, simplified search engine
      subroutine TMscore8_search(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=1000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)
      dimension iL0(nmax)
      common/dpcommon/dpout,ivalue

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.3.0)d0_search=3.0
    
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************

      if(ivalue.eq.0) n_init=2

      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        k=0
        do i=1,iL_max,40        !this is the simplification!
          k=k+1
          iL0(k)=i
        enddo
        if(iL0(k).lt.iL_max)then
          k=k+1
          iL0(k)=iL_max
        endif
        n_shift=k
        do 300 i_shift=1,n_shift
           iL=iL0(i_shift)
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun8    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
***   dis<8, but same search engine
      subroutine TMscore8(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=1000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.3.0)d0_search=3.0
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max    !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun8    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun8
      PARAMETER(nmax=1000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score,score_max
      common/d8/d8

      d_cp=d
  116 continue

      n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d_cp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         if(dis.le.d8)then
            score_sum=score_sum+1/(1+(dis/d0)**2)
         endif
      enddo

       if(n_cut.gt.0.and.n_cut.lt.3) then
          d_cp=d_cp+0.5
         goto 116
       endif
      score=score_sum/float(nseqB) !TM-score

      return
      end

*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
***  normal TM-score:
      subroutine TMscore(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=1000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
c      d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.3.0)d0_search=3.0
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun
      PARAMETER(nmax=1000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score

      d_cp=d
         
  10  continue

      n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d_cp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         score_sum=score_sum+1/(1+(dis/d0)**2)
      enddo

      score=score_sum/float(nseqB) !TM-score

      if(n_cut.gt.0.and.n_cut.lt.3) then
        d_cp=d_cp+0.5
        goto 10
      endif

      return
      end

********************************************************************
*     Dynamic programming for alignment.
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
********************************************************************
      SUBROUTINE DP(NSEQ1,NSEQ2)
      PARAMETER(nmax=1000)
      LOGICAL*1 DIR
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      dimension DIR(0:nmax,0:nmax),VAL(0:nmax,0:nmax)
      REAL H,V
      common/dpcommon/dpout,ivalue
      
***   initialize the matrix:
      val(0,0)=0
      do i=1,nseq1
        dir(i,0)=.false.
        val(i,0)=0
      enddo
      do j=1,nseq2
        dir(0,j)=.false.
        val(0,j)=0
        invmap(j)=-1
      enddo

***   decide matrix and path:
      DO j=1,NSEQ2
        DO i=1,NSEQ1
          D=VAL(i-1,j-1)+SCORE(i,j)
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          
          IF((D.GE.H).AND.(D.GE.V)) THEN
            DIR(I,J)=.true.
            VAL(i,j)=D
          ELSE
            DIR(I,J)=.false.
            if(V.GE.H)then
              val(i,j)=v
            else
              val(i,j)=h
            end if
          ENDIF
        ENDDO
      ENDDO
      
***   extract the alignment:
      i=NSEQ1
      j=NSEQ2
      DO WHILE((i.GT.0).AND.(j.GT.0))
        IF(DIR(i,j))THEN
          invmap(j)=i
          i=i-1
          j=j-1
        ELSE
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          IF(V.GE.H) THEN
            j=j-1
          ELSE
            i=i-1
          ENDIF
        ENDIF
      ENDDO

      j=nseq2
      do while(invmap(j).lt.1)
         j=j-1
      enddo
      dpout=val(invmap(j),j)
      
c^^^^^^^^^^^^^^^Dynamical programming done ^^^^^^^^^^^^^^^^^^^
      return
      END

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
      PARAMETER(nmax=1000)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(nmax), x(3, nmax), y(3, nmax), u(3, 3), t(3), rms, sigma
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

