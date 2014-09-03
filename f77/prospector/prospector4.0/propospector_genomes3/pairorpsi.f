c	like pmfhomgsororienbaln3590.f but uses 50a1 input of msac	like pmfhomgsor.f but includes orientational information
c	like homgsavlon2301 but uses 'OR' criterion for fragment selection
c	move to 3. selection criteria
c	movees selection to 3. std units rather than 2.5
c	 like pmfhomgslong_fixed.f but uses SCALEAV2
c	9/8 gap problem fixed
c	like pmfprodfull2.f but weights contacts by their random probability
c	like  pmfhom30avfull but uses std cutoffs to select the structures
c	
c	like pmfhomdist5wtb4.f but uses homology of 0.3
c	using different weighting for homologous sequence (accounts
c	for gaps
cc	like pmfnat_dist5wtb3 but includes homologous sequences
c	pmfnat_dist5wt.f constuctions whole contact map fragments.
c	threads all local pieces independent of size.
c	not just small to big
C	*********************************************************
c	*       builds  sequence based potential		*
C	*	ALL coordinates are in ANGSTROMS		*
C	*********************************************************
	PARAMETER(NMAX=1700)
	PARAMETER(NCASES=5700)
        PARAMETER(MAXSEQ=300)
	PARAMETER(NRESMAX=530)
	PARAMETER(IWIN=6)
	PARAMETER(IWIN2=13)
	CHARACTER*3 AA3(0:21),AD3
	CHARACTER*1 AD
	CHARACTER*6 TARGET
	CHARACTER*4 GENOME
	CHARACTER*5 NAME2,NAME,NAMET
	CHARACTER*5 NAMESELECT
	CHARACTER*255 F3590,F3590ALN,NAMES(200)

	LOGICAL HOMMAP
	LOGICAL *1 TOUCH(-1:NRESMAX,-1:NRESMAX)

	COMMON/NAMES/NAME(NCASES)
	COMMON/C/NDATA
	COMMON/CONT/NC(NMAX,NCASES),ICON(15,NMAX,NCASES)
	COMMON/CONT2/LOR(15,NMAX,NCASES)	
        COMMON/SEQ/NRESSEQ,NSEQ_ORI,KCODE(NRESMAX,MAXSEQ)
	COMMON/SEQ2/IMUT(0:19,0:19)
	COMMON/SEQ3/NSUC(NRESMAX)
	COMMON/SEQ4/NEX(0:19),IPOS(NRESMAX,0:19)
	COMMON/SCALE/AEXP(NRESMAX,NRESMAX),AOBS(-4:NRESMAX,NRESMAX,3)
	COMMON/SEQHOM/HOMMAP(NCASES)
	COMMON/SCALEPAIR/APABL(-1:19,-1:19,3)
	COMMON/SCALEPAIR2/AN(NRESMAX,NRESMAX,3)

	DIMENSION JCODE(NMAX)
	DIMENSION AMAP(NRESMAX,NRESMAX)
	DIMENSION ICODE(NMAX,NCASES),NRES(NCASES)
c	=================================
c	directories
	common/namesize/nsize		
	character*255 home,templatedat
	common/home/nhome,home
c	=================================
c	======================================================
	DATA AA3/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &		     'MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP','CYX','HIX'/

c	======================================================

c	======================================================
c      			 PSIBLAST INPUTS
        logical*1 tncbi(0:22)
        common/profile/xm(0:19,0:nmax,0:ncases)
        dimension ix(0:22),ncbi(22)
	character*1 aa(0:19),ancbi(1:22)
        DATA AA/ 'G','A','S','C','V','T','I','P',
     &               'M','D','N','L',
     &               'K','E','Q','R',
     &               'H','F','Y','W'/
        data ancbi/'A','B','C','D','E','F','G','H','I','K','L','M','N','P',
     &  'Q','R','S','T','V','W','X','Y'/

        do jres=1,22
        tncbi(jres)=.false.
        do ires=0,19
        if(aa(ires).eq.ancbi(jres))then
        tncbi(jres)=.true.
                ncbi(jres)=ires
                go to 22
        end if
        end do
22      continue
        end do
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

c	&&&&&&& templatedir:    &&&&&&&&&&
	open(unit=1,file='templatedir')
	read(unit=1,fmt='A255')templatedat
	nadr=0
	do i=255,1,-1
	if(templatedat(i:i) .eq. ' ')nadr=i
	end do
	nadr=nadr-1
c	&&&&&&& templatedir:    &&&&&&&&&&

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

	read(5,1112)target
1112	format(a6)
	Open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)

	open(unit=3,file='completedpairorpsi'//target)

	do i=-1,19	
	do j=-1,19
	do it=1,3
	apabl(i,j,it)=0.
	end do
	end do
	end do
	OPEN(unit=26,file=home(1:nhome)//'strlist//SCALET4')
	do it=1,3
	read(26,*)
	read(26,*)
	do i=0,19
	read(26,17) ad3, (apabl(i,j,it),j=0,19)
	enddo
	end do
17	FORMAT(A3,20F5.1)
	close(26)

	OPEN(unit=20,file=home(1:nhome)//'strlist//BLOSUM_INT')
	read(20,*)
	do i=0,19
	read(20,*)ad,(imut(i,j),j=0,19)
	end do
	OPEN(unit=50,file='LIST.jan2005good35')
	read(50,*)NDATA
	NTOT=0
	DO IK=1,NDATA
	read(50,221)name(ik)
221	format(a5)
	OPEN(unit=30,file='/library/orien6/'//NAME(ik)//'.SEQ')
	rewind(30)
	READ(30,2)NAME2,Nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
	read(30,*)(icode(i,ik),i=1,nres(ik))
	do i=1,nres(ik)
	if(icode(i,ik).gt.20)then
	icode(i,ik)=0
	stop
	end if
	end do
	close(30)
2	format(a5,1x,i5)
	do i=1,nres(ik)
	do ires=0,19
	xm(ires,i,ik)=0.
	end do
	end do
        open(unit=31,file=templatedat(1:nadr)//'pdb_jan05/profile/'//name(ik)//'.mtx')
        read(31,*)nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
        do id=1,12
        read(31,*)
        end do
        do i=1,nres(ik)
        read(31,*)(ix(ires),ires=0,22)
                do jres=1,22
                if(tncbi(jres))then
                ires=ncbi(jres)
                xm(ires,i,ik)=float(ix(jres))/100.
                end if
        end do
        end do
	END DO
	do ik=1,ndata
	NTOT=NTOT+NRES(ik)
	end do
	close(50)
	call contb(NRES)

	OPEN(unit=60,file='LIST.targ'//target)
	read(60,*)ndataseq,genome
	OPEN(unit=25,file='problems')
	do 10000  iik=1,ndataseq
	READ(60,33)NAMES(IIK)
33	format(a)
	OPEN(UNIT=6,FILE='outhomgsororien.'//NAMES(IIK)(1:nsize))
	nunit=31
	OPEN(UNIT=31,FILE=f3590(1:n3590)//names(iik)(1:nsize)//f3590aln(1:naln90))
	call read_msa(nunit)
	itest=0
	do i=1,nresseq
	jcode(i)=kcode(i,1)
	end do
	close(31)
c	=======================================
c	read in of homologous sequence homologous sequences
  
	ntot=0
	nz=0
	do ik=1,ndata
	NTOT=NTOT+NRES(ik)
	nz=nz+1
	end do
	do jk=1,ndata
	hommap(jk)=.true.
	end do

	mres=nresseq
	do iseq=0,19
	nex(iseq)=0
	end do
	do i=1,mres
	do iseq=0,19
	if(jcode(i) .eq. iseq)then
	nex(iseq)=nex(iseq)+1
	ipos(nex(iseq),iseq)=i
	end if
	end do
	end do
	do iseq=0,19

	end do
	CALL EHMSTAT(NAME,NRES,ICODE)

c	SEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQSEQ
	mres=nresseq
	OPEN (UNIT=7,FILE=home(1:nhome)//genome//'PAIROR/'//NAMES(iik)(1:nsize)//'.PAIRORPSI')
	rewind(7)
	icorrect=0
	ineg=0
	do it=1,3
	do i=1,mres
	iseq=jcode(i)
	write(7,*)i,aa3(iseq)
	WRITE(7,171) (an(i,j,it),j=1,mres)
	write(7,*)'================================='
171	FORMAT(25F5.1)
	end do
	end do
	close(7)
c	----------------------
	write(3,3)names(iik)(1:nsize)
3	format(a)
10000 	continue
	STOP
	END
	SUBROUTINE EHMSTAT(NAME,NRES,ICODE)
	PARAMETER(NMAX=1700)
        PARAMETER(MAXSEQ=300)
	PARAMETER(NCASES=5700)
	PARAMETER(NRESMAX=530)
	PARAMETER(IWIN=6)
	PARAMETER(IWIN2=13)

	LOGICAL HOMMAP,NOGAP

	CHARACTER*5 NAME(NCASES)

	DIMENSION ICODE(NMAX,NCASES),NRES(NCASES)
	DIMENSION IMAP(-5:NMAX,-5:NMAX,3)
	DIMENSION IHIST(NRESMAX)	
	DIMENSION JCODE(NRESMAX)
        DIMENSION AMUT(-1:19,0:NRESMAX),ANAT(NRESMAX)
        DIMENSION AHIST(-IWIN2:NRESMAX),AHIST2(-IWIN2:NRESMAX),ALEFT(-IWIN2:NRESMAX)

	COMMON/C/NDATA
	COMMON/CONT/NC(NMAX,NCASES),ICON(15,NMAX,NCASES)
	COMMON/CONT2/LOR(15,NMAX,NCASES)	

	COMMON/SEQHOM/HOMMAP(NCASES)
        COMMON/SEQ/NRESSEQ,NSEQ_ORI,KCODE(NRESMAX,MAXSEQ)
	COMMON/SCALEPAIR/APABL(-1:19,-1:19,3)
	COMMON/SCALEPAIR2/AN(NRESMAX,NRESMAX,3)
	COMMON/SEQ2/IMUT(0:19,0:19)
	COMMON/SEQ3/NSUC(NRESMAX)
	COMMON/SEQ4/NEX(0:19),IPOS(NRESMAX,0:19)
	COMMON/SCALE/AEXP(NRESMAX,NRESMAX),AOBS(-4:NRESMAX,NRESMAX,3)
	COMMON/SEQ5/NOGAP(NRESMAX,1)
        COMMON/PROFILE/XM(0:19,0:NMAX,0:NCASES)


	mres=nresseq
	do i=1,mres
        AHIST(I)=0.
        AHIST2(i)=0.
	IHIST(I)=0
	JCODE(I)=KCODE(I,1)
		do j=1,mres
			do it=1,3
			an(j,i,it)=0.
			aobs(j,i,it)=0.
			end do
		end do
	end do
        do i=1,mres
                do iseq=0,19
                amut(iseq,i)=0.
                end do
        nden=0
        do ic=1,nseq_ori
                ires=kcode(i,ic)
        nden=nden+1
                if(ires.gt.-1)then
                amut(ires,i)=amut(ires,i)+1
                end if
        end do  
                do iseq=0,19
                amut(iseq,i)=amut(iseq,i)/nden
                end do
        end do

	IHOM=1
	do 150 i=1,mres
	anat(i)=0.
	nsuc(i)=0
	NOGAP(I,IHOM)=.TRUE.
	ii=0
		do ii=-iwin,iwin
		iv=i+ii
			if(iv .gt. 0 .and. iv.le.mres)then
    	jres=jcode(iv)
		if(jcode(iv).eq. -1)then
				NOGAP(I,IHOM)=.FALSE.
			go to 150
		else
			do ires=0,19
			anat(i)=anat(i)+amut(ires,iv)*amut(ires,iv)
			end do
		end if
		end if
		end do
	if(anat(i).lt.0)NOGAP(I,IHOM)=.false.
150	continue


	do 700 jk=1,NDATA
	if( .not.hommap(jk))go to 700
c	sequence is ik
c	structure is jk
	DO 161 Ij=1,nres(jk)
	do i=1,mres
	if(nogap(i,ihom))then
		aSIM=0.
		is1=jcode(i)
		js1=icode(ij,jk)
	if(is1 .eq.js1)then
		do ii=-iwin,iwin
		iv=i+ii
		jv=ij+ii
		if(iv .gt. 0 .and. iv.le.mres)then
		if(jv .gt. 0 .and. jv.le.nres(jk))then
c		is1=jcode(iv)
c		js1=icode(jv,jk)
		do jres=0,19
		asim=asim+amut(jres,iv)*xm(jres,jv,jk)
		end do
		end if
		end if
		end do
	ahist(i)=ahist(i)+asim
	ahist2(i)=ahist2(i)+asim*asim
	IHIST(I)=IHIST(I)+1
	end if
	end if
	end do
161	CONTINUE
700	CONTINUE

	DO I=1,MRES
	if(nogap(i,ihom))then
	if(ihist(i).gt.0)then
	nd=ihist(i)
	ah=ahist(i)/nd
	ah2=ahist2(i)/nd
	std=sqrt(ah2-ah*ah)
	aleft(i)=ah+3*std
	end if
	END IF
	END DO



	do 1000 jk=1,NDATA
	if(.not.hommap(jk))go to 1000
	isimt=0
	do i=1,nres(jk)+1
	do j=i,nres(jk)+1
	do it=1,3
	imap(i,j,it)=0
	imap(j,i,it)=0
	end do
	end do
	end do
	do ij=1,nres(jk)
	IF(NC(IJ,JK).GT.0)THEN
		do kk=1,nc(ij,jk)
		jj=icon(kk,ij,jk)
		it=lor(kk,ij,jk)
		imap(ij,jj,it)=1
		end do
	end if
	end do
c	sequence is ik
c	structure is jk
	DO 171 Ij=1,nres(jk)
	iseq1=icode(ij,jk)
	IF(NC(IJ,JK).GT.0)THEN
	do 170 kk=1,nc(ij,jk)
		jj=icon(kk,ij,jk)
		it=lor(kk,ij,jk)
		imap(ij,jj,it)=1
	if(iabs(jj-ij) .gt.iwin)then
		iseq2=icode(jj,jk)
	if(nex(iseq1) .gt. 0. and. nex(iseq2).gt.0)then
	DO 169  I1=1,NEX(ISEQ1)
		i=ipos(i1,iseq1)
	DO 168 J1=1,NEX(ISEQ2)
		j=ipos(j1,iseq2)
	IF(NOGAP(i,ihom) .and. nogap(j,ihom))then
	IF(IABS(I-J).GT. IWIN)THEN
			ASIM1=0.
			do ii=-iwin,iwin
			iv=i+II
			jv=ij+ii
			if(iv .gt.0 .and. iv.le.mres)then
			if(jv .gt. 0 .and. jv.le.nres(jk))then
c			is1=jcode(iv)
c			js1=icode(jv,jk)
c			ISIM1=ISIM1+IMUT(IS1,JS1)
                        do jres=0,19
                        asim1=asim1+amut(jres,iv)*xm(jres,jv,jk)
                        end do
  
			end if
			end if
			end do
			ASIM2=0.
			do ii=-iwin,iwin
			iv=j+II
			jv=jj+ii
			if(iv .gt.0 .and. iv.le.mres)then
			if(jv .gt.0 .and. jv.le.nres(jk))then
c			is1=jcode(iv)
c			js1=icode(jv,jk)
c			ISIM2=ISIM2+imut(is1,js1)
                        do jres=0,19
                        asim2=asim2+amut(jres,iv)*xm(jres,jv,jk)
                        end do
  
			end if
			end if
			end do
			is1=jcode(i)
			js1=jcode(j)
	if(asim1 .ge. aleft(i) .or. asim2 .ge.aleft(j))then
			Do IV=-5,5
			iz=i+iv
			ijz=ij+iv
			if(iz .gt. 0 .and. ijz .gt.0)then
			if(iz .le. mres .and. ijz .le.nres(jk))then
			iz1=jcode(iz)
			ijz1=icode(ijz,jk)
	if(imut(iz1,ijz1) .gt. 0 )then
			do jv=-5,5
			jz=j+jv
			jjz=jj+jv
	if(jz .gt.0 .and. jz.le.mres)then
	if(jjz .gt.0 .and. jjz.le.nres(jk))then
			jz1=jcode(jz)
			jjz1=icode(jjz,jk)
		js1=kcode(jz,1)
	if(asim1 .gt. 0 .and. asim2 .gt. 0)then
	if(imut(jz1,jjz1) .gt. 0)then
	sf=asim1*asim2/(anat(i)*anat(j))
	sf=0.2+sqrt(sf)
	do it=1,3
	aobs(iz,jz,it)=aobs(iz,jz,it)+sf*imap(ijz,jjz,it)
	aobs(jz,iz,it)=aobs(jz,iz,it)+sf*imap(ijz,jjz,it)
	end do
			end if
			end if
			end if
			end if
			end do
			end if
			end if
			end if
			END DO
		nsuc(i)=nsuc(i)+1
		nsuc(j)=nsuc(j)+1
		end if
	END IF
	END IF
168 	continue
169	continue
	end if
	end if
170	CONTINUE
	end if
171	CONTINUE
1000	continue

1001	continue

	atot=0.
	do i=1,mres
	do j=1,mres
	do it=1,3
		atot=atot+aobs(i,j,it)
	end do
	end do
	end do
	anorm=atot/(3*mres*mres)
	write(6,*)'number of observations',anorm,atot
	isum=0
	do i=1,mres
		do j=1,mres
	do it=1,3	
	if(aobs(i,j,it) .gt. 0.)then
	an(i,j,it)=an(i,j,it)-alog(aobs(i,j,it)/anorm)
	an(j,i,it)=an(j,i,it)-alog(aobs(i,j,it)/anorm)
		isum=isum+1
		else
	iseq=jcode(i)
	jseq=jcode(j)
		an(i,j,it)=an(i,j,it)+apabl(iseq,jseq,it)
		an(j,i,it)=an(j,i,it)+apabl(iseq,jseq,it)
		end if
	end do
	end do
	end do
2000	continue
	do i=1,mres
	do j=1,mres
	do it=1,3
	an(i,j,it)=an(I,j,it)/(2)
	end do
	end do
	end do
	RETURN
	END 	
	
	subroutine contb(NRES)
	PARAMETER(nmax=1700)
	PARAMETER(NCASES=5700)
	CHARACTER*5 NAME
	common/names/name(ncases)
	common/c/ndata
	common/cont/nc(nmax,ncases),icon(15,nmax,ncases)
	common/cont2/lor(15,nmax,ncases)
	DIMENSION nres(ncases)
C	*********************************************************
c	*        constructs the contact map library		*
c	* obtained for NEWDATA base and is the mean position	*
C	*	ALL coordinates are in ANGSTROMS		*
C	*********************************************************


	DO IK=1,NDATA
	write(71,*)ik,name(ik)
	OPEN(unit=71,file='/library/orien6fit/'//NAME(ik)//'.FITCONT_OR')
	rewind(71)
	do i=1,nres(ik)
	read(71,*)ii,nc(ii,ik),(icon(j,ii,ik),j=1,nc(ii,ik))
	nt=0
	do j=1,nc(ii,ik)
        if(icon(j,ii,ik).gt. nmax)icon(j,ii,ik)=0
        end do


	read(71,*)ii,nct,(lor(j,ii,ik),j=1,nct)
	end do
	CLOSE(71)

	END DO

	RETURN
	end 
c	##########################################################
c========================================================================
      subroutine read_msa(nunit)
C------------------------------------------------------------------------
      PARAMETER     (MAXSEQ=300)
      PARAMETER     (MAXRES=530)
      PARAMETER	    (NRESMAX=530)
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

