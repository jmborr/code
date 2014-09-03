	OPEN(UNIT=9,FILE='stat_domains')
	CALL r14STAT()
	open(unit=1,file='finish_domains')
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
	character*1 a3
	CHARACTER*6 a6
	CHARACTER*60 head
	CHARACTER*255 NAMEP,NAMEV(0:ncases),namep5(0:ncases)
	CHARACTER*1 AA(0:19),ad
	CHARACTER*21 TEXT
	CHARACTER*5 name5
	character*4 genome
	CHARACTER*6 prodin
	DIMENSION NRES(0:NCASES)
	dimension ixm(0:maxres),idomb(50),idome(50)
	dimension idom(50),len(50)
        parameter  (maxseq =3000)	
        common/seq/nresseq,nseq_ori,kcode(maxres,maxseq)	
	character*255 f3590,f3590aln
c	=================================
c	directories
	character*255 home,templatedat
c	=================================
	
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
c	======================================================
c	&&&&&&& home:    &&&&&&&&&&
	open(unit=1,file='homedir')
	read(unit=1,fmt='A255')home
	nhome=0
	do i=255,1,-1
	if(home(i:i) .eq. ' ')nhome=i
	end do
	nhome=nhome-1
c	&&&&&&& home:    &&&&&&&&&&

	open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)

	read(5,13)prodin
13      format(a6)
	OPEN(unit=4,file='LIST.targ'//prodin)
	read(4,*)nd,genome

	text='ATOM     1   CA  GLY '

	OPEN(unit=33, file='STRUCTURE_rap3orienrev')

	rav=0
	read(33,*)ndata2
	iden=0
	rmsmin=0.
	itot=0
	favo=0.
	atc=0
	ntarg=0

	do 1001 IIk=1,NDATA2
        read(33,9432)iselected,namep(1:nsize)
9432    format(i6,3x,a)
	OPEN(UNIT=31,FILE=f3590(1:n3590)//namep(1:nsize)//f3590aln(1:naln90))
	rewind(31)
	nunit=31
	call read_msafasta(nunit)
	close(31)
	mres=nresseq
        OPEN(unit=3,file=namep(1:nsize)//'.threadrap3orienrev')
	OPEN(unit=34,file=home(1:nhome)//genome//'pdborienrev/'//namep(1:nsize)//'rap3orienrev1.pdb')
	read(3,*)nthtot
	! decide on what max threshold will be used. For convenience purposes, we will use
	nth=0
	do 2000 ith=1,nthtot
c	identify alignment regions on a per template basis that hits something significant:
	read(3,18)name5,zth
18	format(a5,1x,1f7.3)
	if(zth.lt.7)go to 2001
	nth=nth+1
	do i=1,mres
	ixm(i)=0
	end do
	read(34,15)name5,ic
15	format(a5,1x,i4,1x,1f7.3,1x,i3)
	do id=1,ic
        read(34,1033)a33,i
1033	FORMAT(A20,2X,I4)
	ixm(i)=1
	end do
	read(34,*)
	ib=0
	do i=1,mres
		if(ixm(i).gt.0)then
			ib=i
			go to 2
		end if
	end do
2       continue
	ie=mres
	do i=mres,1,-1
	if(ixm(i).gt.0)then
		ie=i
		go to 3
		end if
	end do
3   	continue
	idomb(nth)=ib
	idome(nth)=ie
2000	continue
2001	continue
	ntho=nth
!   how many domains/aligned regions are there????
	do i=1,mres
		ixm(i)=0
	end do
	! perhaps better to keep orignal ranking from most to least confident...
c   if the alignment >80% of the chain length no need to rebuild profiles
	if(nth.eq.0)then
		ndom=1
		idomb(1)=1
		idome(1)=mres
		do i=1,mres
		ixm(i)=1
		end do
		go to 59
	end if


	do ith=1,nth
	len(ith)=idome(ith)-idomb(ith)+1
	end do

	do ith=1,nth-1
	do jth=ith+1,nth
	if(len(ith).lt.len(jth))then
	itt=len(ith)
	len(ith)=len(jth)
	len(jth)=itt

	itt=idomb(ith)
	idomb(ith)=idomb(jth)
	idomb(jth)=itt

	itt=idome(ith)
	idome(ith)=idome(jth)
	idome(jth)=itt
	end if
	end do
	end do
	len(1)=idome(1)-idomb(1)+1
		ith=1
		if(float(len(ith))/mres .gt. 0.8)then
		ndom=1
		idomb(1)=1
		idome(1)=mres
		do i=1,mres
		ixm(i)=1
		end do
		go to 59
		end if

	kth=0
	do ith=1,nth
		jt=0
	kthold=kth
			do i=idomb(ith),idome(ith)
			if(ixm(i).eq.0)then
				if(jt.eq.0)then
				jt=1
				kth=kth+1
				end if
			ixm(i)=kth
			end if
			end do
	if(kth.ne. kthold)then
	kb=0
	ib=idomb(ith)
	ie=idome(ith)
	do i=ib,ie
	if(ixm(i).eq.kth)then
	if(kb.eq.0)then
	kb=1
	idomb(kth)=i
	end if
	idome(kth)=i
	end if
	end do
	end if
	end do

	nnth=kth

c	        expansion of the domains:
c   		here is the leftmost boundary
	ie=0
	do i=1,mres
	if(ixm(i).eq.0)then
	ie=i
	else
		ith=ixm(i)
		go to 32
	end if
	end do
32      continue
	IF(Ie .GT.0)THEN
	if(ie .lt.40)then
		do i=1,ie
		ixm(i)=ith
		end do
	idomb(ith)=1
	else
		nnth=nnth+1
		idom(nnth)=nnth
		idomb(nnth)=1
		idome(nnth)=ie
		do i=1,ie
		ixm(i)=nnth
		enddo
	end if
	END IF
c        here is the expansion to the right most boundary

	ib=mres+1
	do i=mres,1,-1
	if(ixm(i).eq.0)then
	ib=i
	else
		ith=ixm(i)
		go to 33
	end if
	end do
33      continue
	if(ib .lt.mres+1)then
	if(mres-ib+1 .lt.40)then
		do i=ib,mres
		ixm(i)=ith
		end do
	idome(ith)=mres
	else
		nnth=nnth+1
		idom(nnth)=nnth
		idomb(nnth)=ib
		idome(nnth)=mres
	do i=ib,mres
		ixm(i)=nnth
	enddo
	end if
	end if

	NTH=nnth
	ixm(mres+1)=1
	ixm(0)=1
 	do 500 ith=1,nth
	ie=idome(ith)
	ib2=ie+1
	if(ixm(ib2).gt.0)then
	go to 40
	end if
	do i=ib2+1,mres
	if(ixm(i).eq.0)then
	ie2=i
	else
		jth=ixm(i)
		go to 401
		end if
	end do
401	continue
	idiff=ie2-ib2+1
	IF(IDIFF.LT. 40)THEN
!	 add to the existing fragments
		idiff2=idiff/2

		idome(ith)=ie+idiff2
		do i=ie+1,ie+idiff2
			ixm(i)=ith
		end do

		idomb(jth)=ie+idiff2+1
		do i=ie+idiff2+1,ie2
			ixm(i)=jth
		end do
	ELSE
		nnth=nnth+1
		idomb(nnth)=ib2
		idome(nnth)=ie2
		do i=ib2,ie2
		ixm(i)=nnth
		end do
	END IF	

40	continue
c	 now do the left side
        ie2=idomb(ith)-1
	if(ixm(ie2).gt.0)then
	go to 41
	end if

        do i=ie2,1,-1
	if(ixm(i).eq.0)then
	ib2=i
	else
	        jth=ixm(i)
	        go to 411
        end if
        end do
411 	continue

	idiff=ie2-ib2+1
	IF(IDIFF.LT.40)THEN
	idiff2=idiff/2
	idomb(ith)=ie2-idiff2
	do i=ie2,ie2-idiff2,-1
	ixm(i)=ith
	end do
	ie2o=idome(jth)
	idome(jth)=ie2-idiff2-1

	do i=ib2,ie2-idiff2-1
	ixm(i)=jth
	end do

	ELSE
	nnth=nnth+1
	idomb(nnth)=ib2
	idome(nnth)=ie2
	do i=ib2,ie2
	ixm(i)=nnth
	end do
	END IF	
41      continue
500	continue
	ndom=nnth
59 	 continue

	open(unit=60,file=home(1:nhome)//genome//'domains/'//namep(1:nsize)//'.domains')
	nal=0
	do i=1,mres
	if(ixm(i).gt.0)nal=nal+1
	end do
	write(9,9)namep(1:nsize),mres,ntho,ndom,nal
9	format(a,1x,4(i4,1x))
	write(60,*)ndom
	do ith=1,ndom
	ib=idomb(ith)
	ie=idome(ith)
	il=ie-ib +1
	write(9,60)ith,idomb(ith),idome(ith),il
	write(60,60)ith,idomb(ith),idome(ith),il
60	format(4(i5,1x))
	end do
	write(9,*)'	 ====== '
1001	continue


	RETURN
	END 	

      SUBROUTINE READ_MSAFASTA(NUNIT)
C
C------------------------------------------------------------------------

C------------------------------------------------------------------------
      PARAMETER     (MAXSEQ =3000)
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

