	OPEN(UNIT=9,FILE='stat_verygood')
	CALL r14STAT()
	open(unit=1,file='finish_verygood')
	rewind(1)
	write(1,*)'1'
	close(1)
	STOP
	END


	
	SUBROUTINE r14STAT()
	PARAMETER(NMAX=3000)
	PARAMETER(NCASES=4200)
	PARAMETER(MAXRES=1800)
c	PARAMETER(REJECT=25.)	
	CHARACTER*5 NAME(0:NCASES),NAMESELECT
	CHARACTER*5 NAMETARG
	CHARACTER*3 NAME3,TER
	character*20 a33
	character*1 a3,qdir
	CHARACTER*6 a6
	CHARACTER*60 head
	CHARACTER*5 NAMEH(500),namel
	CHARACTER*255 NAMEP,NAMEV(0:ncases),namep5(0:ncases)
	CHARACTER*1 AA(0:19),ad
	CHARACTER*21 TEXT
	CHARACTER*5 name5
	character*4 genome
	CHARACTER*6 prodin
	LOGICAL*1 DIR,kselect,hommap(NCASES),TM1(0:NMAX,5)
        logical*1 tgood

	DIMENSION TGOOD(0:nmax),TGOOD2(0:nmax)
	double precision u,t,umin(3,3),tm(3)
	common/matrix/u(3,3),t(3)
c	DIMENSION AMUT(-1:19,0:NMAX)	
	DIMENSION XSC(3,NMAX),XSCT(3,NMAX)
	DIMENSION XSC2(3,NMAX),XSCT2(3,NMAX)	
	COMMON/C/NDATA,NDATA2
	DIMENSION XA(3,NMAX,0:NCASES)
        COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *   LIGNMENT(0:2*MAXRES), GAP_OPEN, GAP_EXTN
	COMMON/DIR2/DIR(0:NMAX,0:NMAX)    
	COMMON/SEQALL/IMUT(-1:19,-1:19)	
	COMMON/ALIGNMENT/INVMAP(0:NMAX,0:NCASES,4)
	COMMON/MAPPINGS2/JK,LST
	COMMON/nam/name3
	COMMON/MAPPINGS/JCODE(-1:MAXRES),IMAP(0:MAXRES)     	
	dimension isz(3000)	
	DIMENSION ICODE(NMAX,0:NCASES)
	DIMENSION NRES(0:NCASES)
	DIMENSION IMAPOLD(-5:NMAX,NCASES)
	DIMENSION AMUT(-1:19,0:NMAX)
	dimension ist(0:2),ist2(0:2)
	common/rotate/mres,xs(3,NMAX)
        dimension energ(0:ncases),nm(0:ncases),item(0:ncases)
	dimension zminrms(0:ncases),sa(3),sb(3)
	dimension zminrmso(0:ncases)
	dimension ires1(0:nmax),jres1(0:nmax)
        dimension nc(0:nmax), icon(15,0:nmax)
        dimension jc(0:nmax), jcon(15,0:nmax)
	dimension  map(0:nmax)
	dimension conpred3(-4:nmax,-4:nmax)
	dimension fav(0:4),nstat(0:4),ratioc(0:4)
        dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
        dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax),x3(3)
	dimension i1(nmax),i2(nmax)
c	=================================
c	directories
	character*255 home
c	=================================

        common/para/d,d0	
         parameter  (maxseq =3000)	
        common/seq/nresseq,nseq_ori,kcode(maxres,maxseq)	
	character*255 f3590,f3590aln
	
	DATA AA/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/
c	======================================================
c	readin  the directories"
c	&&&&&&& ENTER 3590 set:    &&&&&&&&&&
	open(unit=1,file='3590')
	read(unit=1,fmt='A255')f3590
	n3590=0
	do i=255,1,-1
	if(f3590(i:i) .eq. ' ')n3590=i
	end do
	n3590=n3590-1
	open(unit=1,file='3590aln')
	read(unit=1,fmt='A255')f3590aln
	naln90=0
	do i=255,1,-1
	if(f3590aln(i:i) .eq. ' ')naln90=i
	end do
	naln90=naln90-1
c	&&&&&&& home:    &&&&&&&&&&
	open(unit=1,file='homedir')
	read(unit=1,fmt='A255')home
	nhome=0
	do i=255,1,-1
	if(home(i:i) .eq. ' ')nhome=i
	end do
	nhome=nhome-1
c	&&&&&&& home:    &&&&&&&&&&


c	======================================================
	open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)

	open(unit=55,file='nametarg1')
        read(55,13)prodin
13      format(a6)
	OPEN(unit=4,file='LIST.targ'//prodin)
	read(4,*)nd,genome

	text='ATOM     1   CA  GLY '
	OPEN(unit=20,file=home(1:nhome)//'strlist2/BLOSUM_INT')
	read(20,*)
	do i=0,19
        read(20,*)ad,(imut(i,j),j=0,19)
	end do
	close(20)
	open(unit=1,file='listbest')
	read(1,*)nfz
	do ik=1,nfz
	read(1,1)namep5(ik)(1:nsize)
1	format(a)
	end do
	OPEN(unit=33, file='STRUCTURE_rap3orienrev')
	OPEN(unit=63, file='list_verygood')
	rav=0
	read(33,*)ndata2
	iden=0
	rmsmin=0.
	itot=0
	favo=0.
	atc=0
	ntarg=0

	do 1001 IIk=1,NDATA2
        read(33,9432)iselected,namep(1:nsize),nameselect,zselect,head
9432    format(i6,3x,a,3x,a5,1f7.3,,3x,a60)
	do jk=1,nfz
	if(namep(1:nsize).eq.namep5(jk)(1:nsize)) go to 1001
	end do
	OPEN(UNIT=31,FILE=f3590(1:n3590)//namep(1:nsize)//f3590aln(1:naln90))
	rewind(31)
	call read_msafasta(nunit)
	close(31)
	mres=nresseq
	nseq_ori=1
	nc0=0
	do i=1,mres
	do j=1,mres	
	conpred3(i,j)=0
	end do
	end do
        OPEN(unit=3,file=namep(1:nsize)//'.threadrap3orienrev')
	OPEN(unit=34,file=home(1:nhome)//genome//'pdborienrev/'//namep(1:nsize)//'rap3orienrev1.pdb')
	read(3,*)nth
	if(nth .lt.2) go to 1001
	nth=2
c	if(nth .gt.1) go to 1001
	do 2000 ithread=1,nth
	do i=1,mres
	tm1(i,ithread)=.false.
	end do
	read(34,15)name5,ic
c	read(3,25)name5,zth
25	format(a5,1x,1f7.3)
15	format(a5,1x,i4,1x,1f7.3,1x,i3)
201	continue
        rewind(30)
	OPEN(unit=30,file='/library/orien6/'//NAME5//'.SEQ')
	rewind(30)
        read(30,2)name2,NSEQ3
2       format(a5,1x,i5)
        DO J=1,nSEQ3
        map(j)=0
        end do
	if(ithread.eq.1)then
	do id=1,ic
        read(34,1033)a33,i,(xa(l,i,ithread),l=1,3),j
	tm1(i,ithread)=.true.
	map(j)=i
	n2(id)=i
	i2(i)=id
	x2(id)=xa(1,i,ithread)
	y2(id)=xa(2,i,ithread)
	z2(id)=xa(3,i,ithread)
	end do
	L2=IC
	elseif(ithread.eq.2)then
	L1=IC
	do id=1,ic
        read(34,1033)a33,i,(xa(l,i,ithread),l=1,3),j
1033    FORMAT(a20,2X,I4,1x,3X,3F8.3,1x,i4,1x,a3)
	tm1(i,ithread)=.true.
	n1(id)=i
	i1(i)=id
	x1(id)=xa(1,i,ithread)
	y1(id)=xa(2,i,ithread)
	z1(id)=xa(3,i,ithread)
	map(j)=i
	end do
	end if

	read(34,*)

        OPEN(unit=71,file='/library/orien6fit/'//name5//'.FITCONT_OR')
        rewind(71)
        do i=1,NSEQ3
        read(71,*)ii,jc(ii),(jcon(j,ii),j=1,jc(ii))
        read(71,*)
        do j=1,jc(ii)
        if(jcon(j,ii).gt. nmax)jcon(j,ii)=0
        end do
        end do
        close(71)
       Do J=1,NSEQ3
         if(map(j).gt.0)then
        i=map(j)
        if(i .gt.0)then
        if(jc(j).gt.0)then
        do jj=1,jc(j)
        ji=jcon(jj,j)
        ij=map(ji)
		if(ij .gt.0)then
        	conpred3(i,ij)=conpred3(i,ij)+5
		end if
        end do

        end if
        end if
        end if
        enddo
2000	continue
	IF(L1.eq.0) go to 1001
        call TMscore(L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,TMS,Rcomm,Lcomm,SS)
	IF(TMS .LT.0.70) go to 1001
	rrms=0.
	nn=0
	reject=d0
	do i=1,mres
	tgood(i)=.false.
	tgood2(i)=.false.
	if(tm1(i,1) .and. tm1(i,2))then
	tgood(i)=.true.
	nn=nn+1
	ii1=i1(i)
	ii2=i2(i)
	xsc(1,nn)=x1(ii1)
	xsct(1,nn)=x2(ii2)
	xsc(2,nn)=y1(ii1)
	xsct(2,nn)=y2(ii2)
	xsc(3,nn)=z1(ii1)
	xsct(3,nn)=z2(ii2)
	end if
	end do

	do i=1,mres
	do j=1,mres
	if(tgood(i).and. tgood(j))then
	go to 10
	else
	conpred3(i,j)=0
	end if
10	continue
	end do
	end do
	nt=0
	do j=1,mres
	do i=1,mres
	if(iabs(i-j).gt.4)then
	if(conpred3(i,j).gt.9)then
	nt=nt+1
	end if
	end if
	end do
	end do
	nt=nt/2
	open(unit=59,file=home(1:nhome)//genome//'goodr1aorienrev/'//namep(1:nsize)//'.goodcon')
	write(59,*)nt
	do i=1,mres
	do j=i+5,mres
	if(conpred3(i,j).gt.9)then
	write(59,*)i,j
	end if
	end do
	end do
	rms=0.
	nn=0
	nn2=0
	do i=1,mres
	if(tgood(i))then
	nn=nn+1
		dd=0.
		do jj=1,3
		dd=dd+(xsc(jj,nn)-xsct(jj,nn))**2
		end do
	if(dd .lt. reject)then
	tgood2(i)=.true.
	rms=rms+dd
	nn2=nn2+1	
	do jj=1,3
	xsc2(jj,nn2)=xsc(jj,nn)
	xsct2(jj,nn2)=xsct(jj,nn)
	end do
	end if
	end if
	end do
c	call rmscalc(nn2,xsc2,xsct2,rrmso,drmso)
	rms=sqrt(rms/nn2)
	write(9,39)namep(1:nsize),L1,L2,LCOMM,TMS, RCOMM,'comparison -<d0 region',reject,nn2,rms
39	format(a,1x,3(i3,1x),2(1f7.3,1x),1x,a22,1x,1f7.3,1x,i3,1x,1f7.3)

	open(unit=59,file=home(1:nhome)//genome//'goodr1aorienrev/'//namep(1:nsize)//'.freeze')
	ni=0
	do i=1,mres
	If(tgood2(i))ni=ni+1
	end do
	write(59,*)ni
	do i=1,mres
	If(tgood2(i))write(59,*)i
	end do
	close(59)

c	now compare to native structure
2001    continue
	ntarg=ntarg+1
	namev(ntarg)(1:nsize)=namep(1:nsize)
	itot=itot+1
	atc=atc+float(LCOMM)/mres
	TMSAV=TMSAV+TMS
	rav=rav+rcomm
	rms0=rms0+rms
	iden=iden+1


1001	continue

	write(9,*)'total number of proteins examined',itot
	write(9,*)'fraction of good cases=',float(iden)/itot
	write(9,*)'average fraction of coverage=',atc/iden
	write(9,*)'average Tm of the good selected fragments', tmsav/iden
	write(9,*)'average rms of the good selected fragments', rav/iden
	write(9,*)'average rmsd <d0= ',rms0/iden

	write(63,*)ntarg
	do ik=1,ntarg
	write(63,63)namev(ik)(1:nsize)
63	format(a)
	end do
	RETURN
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
      subroutine TMscore(L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,TM,Rcomm,Lcomm,SMAX)
      PARAMETER(nmax=3000)
      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0	
      common/align_2/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score

      dimension xa_max(nmax),ya_max(nmax),za_max(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score,score_maxsub
      double precision score,score_max
      double precision score_maxsub

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ****************
      nseqA=L1
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
      enddo
      nseqB=L2
      do i=1,L2
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
      enddo

******************************************************************
*     pickup the aligned residues:
******************************************************************
      k=0
      do i=1,nseqA
         do j=1,nseqB
            if(nresA(i).eq.nresB(j))then
               k=k+1
               iA(k)=i
               iB(k)=j
               goto 205
            endif
         enddo
 205     continue
      enddo
      n_ali=k                   !number of aligned residues
      Lcomm=n_ali
      if(n_ali.lt.1)then
c         write(*,*)'There is no common residues in the input structures'
         TM=0
         Rcomm=0
         return
      endif

************/////
*     parameters:
*****************
***   d0------------->
      d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      if(d0.lt.0.5)d0=0.5
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=5                    !maximum number of iterations
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
      score_maxsub_max=-1       !MAXSUB SCORE
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
           LL=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xx=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yy=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zz=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              xa(j)=xx
              ya(j)=yy
              za(j)=zz
           enddo
           d=d0_search-1
           call score_fun       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              do i=1,nseqA
                 xa_max(i)=xa(i)
                 ya_max(i)=ya(i)
                 za_max(i)=za(i)
              enddo
           endif
          if(score_maxsub_max.lt.score_maxsub)score_maxsub_max=score_maxsub
***   iteration for extending ---------------------------------->
           do 301 it=1,n_it
              LL=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
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
                 xx=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yy=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zz=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
                 xa(j)=xx
                 ya(j)=yy
                 za(j)=zz
              enddo
              if(it.gt.1)d=d+1
              call score_fun    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 do i=1,nseqA
                    xa_max(i)=xa(i)
                    ya_max(i)=ya(i)
                    za_max(i)=za(i)
                 enddo
              endif
 301       continue             !for iteration
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      do i=1,L1
         x1(i)=xa_max(i)
         y1(i)=ya_max(i)
         z1(i)=za_max(i)
      enddo
      TM=score_max
      SMAX=score_maxsub_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun
      PARAMETER(nmax=3000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align_2/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score,score_maxsub
      double precision score,score_max
      double precision score_maxsub

      n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      score_maxsub_sum=0        !Maxsub-score
  
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         score_sum=score_sum+1/(1+(dis/d0)**2)
         if(dis.lt.3.5)then
            score_maxsub_sum=score_maxsub_sum+1/(1+(dis/3.5)**2)
         endif
      enddo
      score=score_sum/float(nseqB) !TM-score
***   for MAXsub-score:
        score_maxsub=score_maxsub_sum/float(nseqB) !MAXsub-score
      return
      end

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc

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



      SUBROUTINE READ_MSAFASTA(NUNIT)
C
C------------------------------------------------------------------------

C------------------------------------------------------------------------
      PARAMETER     (MAXSEQ =1000)
      PARAMETER     (MAXRES=1800)
      PARAMETER	  (NRESMAX=1800)
	PARAMETER	  (KMAX=10000)
      CHARACTER     SEQ*1
      CHARACTER*1   AA1(0:19),SEQD(KMAX)
      DIMENSION SEQ(MAXSEQ,MAXRES)

      common/seq/nresseq,nseq_ori,kcode(nresmax,maxseq)
      DATA AA1/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/
	nunit=31
	read(nunit,*)nseq_ori,nresseq
	if(nresseq.gt.maxres)then 
	nresseq=maxres
	end if
	
	if(nseq_ori.gt.maxseq)nseq_ori=maxseq
	if(nresseq.le.maxres)then
	do i=1,nseq_ori
	read(nunit,1)(seq(i,j),j=1,nresseq)
1	format(50a1)
	end do
	else
	do i=1,nseq_ori
	read(nunit,1)(seqd(j),j=1,nresseq)
	do j=1,maxres
	seq(i,j)=seqd(j)
	end do
	end do
	nresseq=maxres
	end if

	
c	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	do i=1,nseq_ori
	do j=1,nresseq
		kcode(j,i)=-1
		do ires=0,19
			if(aa1(ires) .eq. seq(i,j))then
			kcode(j,i)=ires
			go to 18
			End if
			end do	
18 	continue
	end do
	end do
c	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



99	 return
      	end