	PARAMETER(nmax=1000)
        parameter (maxseq=300)
        parameter (maxres=500)
	PARAMETER(ncases=8700)
	CHARACTER*5 NAME(0:ncases),name2,name5
	character*1 ad
	character*3 amino
	character*20 a20
	DIMENSION ICODE(0:NMAX,0:NCASES),nres(0:ncases)
	common/profile/xm(0:59,0:nmax,0:ncases)
	common/c/ndata,ndata2
	COMMON/S19/jbin(nmax+1,0:ncases)
        common/seq/nresseq,nseq_ori,kcode(nmax,maxseq)
	common/cont/nc(0:nmax,0:ncases),icon(15,0:nmax,0:ncases)
	common/cont2/ior(15,0:nmax,0:ncases)
	COMMON/SEQALL/imut(-1:19,-1:19)	
c	=================================
c	directories
	common/namesize/nsize		
	character*255 home,templatedat
	common/home/nhome,home
c	=================================
c       PSIBLAST
        logical*1 tncbi
	common/psi/tncbi(0:22),ncbi(22)
        dimension ix(0:22),nseq2(ncases)
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
c	=================================

	Open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)
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

	OPEN(unit=20,file=home(1:nhome)//'strlist/BLOSUM_INT')
	read(20,*)
	do i=0,19
	read(20,*)ad,(imut(i,j),j=0,19)
	end do
	close(20)
	OPEN(UNIT=20,FILE='LIST.jul2007good35')


	READ(20,*)NDATA
	NTOT=0
	DO 1100 IK=1,NDATA
	READ(20,1)NAME(ik)
1	format(a5)
1100	continue
	DO IK=1,NDATA
	OPEN(unit=30,file='/local/images/templates-2007081400/orien6/'//NAME(ik)//'.SEQ')
	read(30,2)name2,nres(ik)
	nseq2(ik)=nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
	read(30,*)(icode(i,ik),i=1,nres(ik))
2	format(a5,1x,i5)
	do i=1,nres(ik)
	if(icode(i,ik).gt.19)icode(i,ik)=0.
	end do
	close(30)
	END DO
	close(20)

	DO IK=1,NDATA
	OPEN(unit=71,file='/local/images/templates-2007081400/orien6fitgly/'//NAME(ik)//'.FITCONT_ORGLY')
	do i=1,nres(ik)
	read(71,*)ii,nc(ii,ik),(icon(j,ii,ik),j=1,nc(ii,ik))
	read(71,*)ii,nc(ii,ik),(ior(j,ii,ik),j=1,nc(ii,ik))
        do j=1,nc(ii,ik)
        if(icon(j,ii,ik).gt. nmax)icon(j,ii,ik)=0
        end do
	end do
	close(71)
	END DO


	DO IK=1,NDATA
c	open(unit=14,file='/local/images/templates-2007081400/templatesseqjul07/'//name(ik)//'.sec')
       open(unit=14,file='/local/images/templates-2007090400/templatesseqjul07/'//name(ik)//'.sec')

	read(14,*)
	read(14,*)(jbin(i,ik),i=1,nres(ik))
	close(14)
	END DO

	DO 1101 IK=1,NDATA
        open(unit=59,file='/local/images/templates-2007082700/struprofilejul07/'//name(ik)//'.profile2')
        read(59,*)nseq1
	if(nseq1.gt.nmax)nseq1=nmax
        do i=1,nseq1
        read(59,59)ii,(xm(ires,i,ik),ires=0,59)
59      format(i4,1x,20(1f7.3,1x),/,2(5x,20(1f7.3,1x),/))
        end do
        close(59)
1101	continue
	CALL r14STAT(name,NRES,ICODE)
	
c	END DO

	STOP
	END


	
	SUBROUTINE r14STAT(name,NRES,ICODE)
	PARAMETER(nmax=1000)
	PARAMETER(ncases=8700)
	parameter(maxseq=300)
        parameter (maxres=500)

	CHARACTER*1 AA(0:19)
	CHARACTER*3 AD3,AMINO
	CHARACTER*4 GENOME
	CHARACTER*5 NAME(0:NCASES),NAME2,NAMESELECT
	CHARACTER*6 PRODIN	
	CHARACTER*21 TEXT
	CHARACTER*255 HOME,NAMEP
	CHARACTER*255 F3590,F3590ALN,FE10,FE10ALN,PRED,PREDSEQ
	character*255 fmtx

	logical*1 touch2(-4:nmax,-4:nmax)
	logical*1 tpsec(30,4)
	LOGICAL*1 HOMOLOGY(0:NCASES)

	DIMENSION ICODE(0:NMAX,0:NCASES)
	DIMENSION NRES(0:NCASES)
	DIMENSION JSEQ(0:MAXRES),ESEQ(0:3,0:3)
	DIMENSION CONPRED3(0:NMAX,0:NMAX)
	DIMENSION EN(-2:7),EN2S(-2:7),SD(-2:7)
	DIMENSION E(-2:7),JCODE(-1:MAXRES)
	DIMENSION APABL(-1:19,-1:19)
	DIMENSION ENERG(NCASES,4),ITEM(NCASES,4)
	DIMENSION AN(-1:MAXRES,-1:MAXRES,3)
	DIMENSION AMUT(-1:39,0:MAXRES),AMUT2(-1:39,0:MAXRES)

	logical*1 tncbi
	common/psi/tncbi(0:22),ncbi(22)
	dimension ix(0:22)

	dimension an2(-1:maxres,-1:maxres,3),conpred(-1:maxres,-1:maxres)
	dimension an01(-1:maxres,-1:maxres,3)

	common/home/nhome,home
        common/seq/nresseq,nseq_ori,kcode(nmax,maxseq)
	common/profile/xm(0:59,0:nmax,0:ncases)
	common/mappings2/jk,lst
	common/cont/nc(0:nmax,0:ncases),icon(15,0:nmax,0:ncases)
	common/cont2/ior(15,0:nmax,0:ncases)
	COMMON/S19/jbin(nmax+1,0:nCASES)
	COMMON/C/NDATA,NDATA2
        COMMON /ALIGNPARAM/GAP_OPEN, GAP_EXTN,SCORE(0:MAXRES,0:NMAX)
	common/alignment/invmap(0:nmax,0:ncases,4)	
	COMMON/SEQALL/imut(-1:19,-1:19)	
	common/namesize/nsize		
	DATA AA/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/

	text='ATOM     1   CA  GLY '          
c	======================================================
c	readin  the directories"
c	&&&&&&&& ENTER MTX set	   &&&&&&&&&&
	open(unit=1,file='mtx')
	read(unit=1,fmt='A255')fmtx
	nmtx=0
	do i=255,1,-1
	if(fmtx(i:i) .eq. ' ')nmtx=i
	end do
	nmtx=nmtx-1

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
c	&&&&&&& ENTER e10 set:    &&&&&&&&&&
	open(unit=1,file='e10')
	read(unit=1,fmt='A255')fe10
	ne10=0
	do i=255,1,-1
	if(fe10(i:i) .eq. ' ')ne10=i
	end do
	ne10=ne10-1
	open(unit=1,file='e10aln')
	read(unit=1,fmt='A255')fe10aln
	naln10=0
	do i=255,1,-1
	if(fe10aln(i:i) .eq. ' ')naln10=i
	end do
	naln10=naln10-1
c	&&&&&&& ENTER secondary structure prediction set:    &&&&&&&&&&

	open(unit=1,file='secon')
	read(unit=1,fmt='A255')pred
	npred=0
	do i=255,1,-1
	if(pred(i:i) .eq. ' ')npred=i
	end do
	npred=npred-1
	open(unit=1,file='seconend')
	read(unit=1,fmt='A255')predseq
	npredseq=0
	do i=255,1,-1
	if(predseq(i:i) .eq. ' ')npredseq=i
	end do
	npredseq=npredseq-1
c	======================================================
	OPEN(UNIT=27,file='/local/images/templates-2007081400/orien6fitgly/PAIRGLY',status='old')
	read(27,*)
	read(27,*)
	DO I=0,19
	read(27,17) Ad3,(apabl(I,J),J=0,19)
	end do
17	FORMAT(A3,23F5.1)
	close(27)

        read(5,13)prodin
13      format(a6)
	close(11)
	OPEN(UNIT=9,FILE='stat_r5saorienrev'//prodin)
	OPEN(UNIT=69,FILE='completedr5saorienrev'//prodin)	
	OPEN(unit=33,file='LIST.targ'//prodin)
	read(33,*)ndata2,genome
	do 1001 Ik=1,NDATA2
	read(33,11)namep(1:nsize)
11	format(a)
	do jk=1,ndata
	homology(jk)=.true.
	end do
        open(unit=8,file=home(1:nhome)//'/pdbhomoljul07/'//namep(1:nsize)//'.homol')
        read(8,*)nhom
        do id=1,nhom
        read(8,818)nameselect
818     format(a5)
                do jk=1,ndata
                if(name(jk).eq.nameselect)then
                homology(jk)=.false.
                go to 218
                end if
                end do
218     continue
        end do
	nunit=31
	OPEN(UNIT=31,FILE=fe10(1:ne10)//namep(1:nsize)//fe10aln(1:naln10))
	rewind(31)
	call read_msa(nunit)
	close(31)
	mres=nresseq
	OPEN(UNIT=14,FILE=pred(1:npred)//namep(1:nsize)//predseq(1:npredseq))
	rewind(14)
	do i=1,nresseq
	jcode(i)=kcode(i,1)
	read(14,*)id,amino,iseq
	if(iseq.eq.2)then
	jseq(i)=3
	elseif(iseq .eq.4)then
	jseq(i)=2
	elseif(iseq.eq.3)then
	jseq(i)=1
	else
	jseq(i)=0
	end if
	end do
	close(14)

	do ib=0,3
	do jb=0,3
	eseq(ib,jb)=0
	end do
	end do
        eseq(3,3)=1.
        eseq(2,2)=1.
        eseq(2,3)=-1.
        eseq(3,2)=-1.
	do j=-1,mres
	do i=-1,mres
	do ils=1,3
	an(j,i,ils)=0.
	end do
	conpred(i,j)=0.
	end do
	end do
       OPEN(unit=44,file=home(1:nhome)//genome//'PAIROR/'//namep(1:nsize)//
     &  '.PAIRORPSIGLY3')
        do ils=1,3
       do i=1,nresseq
       read(44,*)
       read(44,167) (an01(i,j,ils),j=1,nresseq)
       read(44,*)
       end do
        end do
       close(44)
c
!	these have to be updated
c	==============================================
	do i=1,mres
		do iseq=0,19
		amut(iseq,i)=0.
		end do
	nden=0
	do ic=1,nseq_ori
		ires=kcode(i,ic)
		if(ires.gt.-1)then
		nden=nden+1
		amut(ires,i)=amut(ires,i)+1
		end if
	end do	
		do iseq=0,19
		amut(iseq,i)=0.5*amut(iseq,i)/nden
		end do
	end do

        open(unit=31,file=fmtx(1:nmtx)//namep(1:nsize)//'.mtx')
	rewind(31)
	read(31,*)mres1
	do id=1,13
	read(31,*)
	end do
	do i=1,mres1
	read(31,*)(ix(ires),ires=0,22)
		do jres=1,22
		if(tncbi(jres))then
		ires=ncbi(jres)+20
		amut(ires,i)=float(ix(jres))/200.
		end if
	end do
	end do
	close(31)
	CWT=1.8
	ap1=0.5*CTW
	do j=1,mres
	do i=1,mres
	do ic=1,nseq_ori
	ires=kcode(i,ic)
	jres=kcode(j,ic)
	do ils=1,3
	an(i,j,ils)=an(i,j,ils)-apabl(ires,jres)
	end do
	end do
	do ils=1,3
	an(i,j,ils)=5*an(i,j,ils)/nseq_ori
c        an01(i,j,ils)=(an0(i,j)+0.25*an01(i,j,ils))/1.25
	an(i,j,ils)=0.9*(an(i,j,ils)-an01(i,j,ils))
	end do
	end do
	end do

	OPEN(UNIT=31,FILE=f3590(1:n3590)//namep(1:nsize)//f3590aln(1:naln90))
	rewind(31)
	call read_msa(nunit)
	rewind(31)	
	close(31)

	do j=-1,mres
	do i=-1,mres
	do ils=1,3
	an2(j,i,ils)=0.
	end do
	end do
	end do

	do i=1,mres
		do iseq=0,19
			amut2(iseq,i)=0.
		end do
		nden=0
		do ic=1,nseq_ori
			ires=kcode(i,ic)
		if(ires.gt.-1)then
		nden=nden+1
			amut2(ires,i)=amut2(ires,i)+1
		end if
		end do	
        	do iseq=0,19
		amut2(iseq,i)=0.5*amut2(iseq,i)/nden
		end do
	end do
c	==============================================

	CWT=0.9

	do j=1,mres
	do i=1,mres
	do ic=1,nseq_ori
	jres=kcode(j,ic)
	ires=kcode(i,ic)
	do ils=1,3
	an2(i,j,ils)=an2(i,j,ils)-apabl(ires,jres)
	end do
	end do
	do ils=1,3
	an2(i,j,ils)=5.*an2(i,j,ils)/nseq_ori
	an2(i,j,ils)=0.45*(an2(i,j,ils)-an01(i,j,ils))
	end do
	end do
	end do

	itarg=0
	name(0)='NULL'
	do jk=1,ndata	
	do k=1,4
	energ(jk,k)=-1000.
	end do
	end do

	
	write(9,*)'==========================='
	write(9,*)namep(1:nsize),' ',mres
	write(9,*)'==========================='
	NSEQ1=mres

	icount=0
	do k=1,4
		en(k)=0.
		en2s(k)=0.
	end do
	do 1000 jk=1,NDATA
	if(homology(jk))then
	NSEQ2=nres(jk)
	icount=icount+1

c	structure is jk
c	sequence is ik

	ils=0
	DO J=0,NSEQ2	
	DO I=0,NSEQ1
	SCORE(I,J)=0.
	END DO
	END DO
c Prepare a SCORE table
c       GAP_OPEN=-5.47
c      GAP_EXTN=-0.16
c      AWT=1.52
c      BWT=0.98
	GAP_OPEN=-7.0
	GAP_EXTN=-0.05
	AWT=1.5
	BWT=0.7
        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1
	ss=0.
	do jres=0,39
	ss=ss+amut(jres,i)*xm(jres,j,jk)
	end do
 	score(ii,j)=ss+eseq(ibinseq,ibin)*bwt+awt
	end do
	end do
C Prepare a SCORE table
	lst=1
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(1)=result
C Prepare a SCORE table
c	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
c	GAP_OPEN=-10.5
c	GAP_EXTN=-0.05
c	GAP_OPEN=-11.0
c	GAP_EXTN=-0.05
c	AWT=0.1
c	BWT=1.0
	GAP_OPEN=-14.0
	GAP_EXTN=-1.05
	BWT=0.4
c	AWT=0.1
	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1
 	score(ii,j)=score(ii,j)+eseq(ibinseq,ibin)*bwt
        if(nn.gt. 0)then
                do jj=1,nc(j,jk)
                ll=icon(jj,j,jk)
	        ils=ior(jj,j,jk)		
                ll=invmap(ll,jk,1)
        score(ii,j)=score(ii,j)+an(i,ll,ils)
                end do
        end if
        ENDDO
	ENDDO
	lst=2	
      	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(2)=result

c Prepare a SCORE table for the original profile
	DO I=0,NSEQ1
	DO J=0,NSEQ2
	SCORE(I,J)=0.
	END DO
	END DO
c        GAP_OPEN=-10.47
c       GAP_EXTN=-0.12
c       AWT=1.52
c       BWT=1.33
	GAP_OPEN=-10.
	GAP_EXTN=-0.1
	AWT=1.3
	BWT=0.8

        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1	
	ss=0.
	do jres=0,19
	ss=ss+amut2(jres,i)*(xm(jres,j,jk))
	end do
	do jres=20,39
	ss=ss+amut(jres,i)*(xm(jres+20,j,jk))
	end do

	score(ii,j)=ss+eseq(ibinseq,ibin)*bwt+awt
	end do
	end do

	lst=3	
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(3)=result		
	GAP_OPEN=-14.5
	GAP_EXTN=-0.75
c	AWT=1.1
	BWT=0.2

c	GAP_OPEN=-14.
c	GAP_EXTN=-0.7
c	AWT=1.1
c	BWT=0.6
	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1	
 	score(ii,j)=score(ii,j)+eseq(ibinseq,ibin)*bwt
        if(nn.gt. 0)then
                do jj=1,nc(j,jk)
                ll=icon(jj,j,jk)
	        ils=ior(jj,j,jk)				
                ll=invmap(ll,jk,3)
        score(ii,j)=score(ii,j)+an2(i,ll,ils)
                end do
        end if
        ENDDO
	ENDDO


	lst=4		
      	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(4)=result
777	continue
c        GAP_OPEN=-5.47
c       GAP_EXTN=-0.16
c       AWT=1.52
c       BWT=0.98
	GAP_OPEN=-7.0
	GAP_EXTN=-0.05
	AWT=1.5
	BWT=0.7
        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	ss=0.
	do jres=0,39
	ss=ss+amut(jres,i)*xm(jres,j,jk)
	end do
 	score(i,j)=ss+eseq(ibinseq,ibin)*bwt+awt
	end do
	end do


C Prepare a SCORE table
	lst=1
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(1)=result-e(1)

c	GAP_OPEN=-11.0
c	GAP_EXTN=-0.05
	GAP_OPEN=-14.0
	GAP_EXTN=-1.05
	BWT=0.4

	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
 	score(i,j)=score(i,j)+eseq(ibinseq,ibin)*bwt
        if(nn.gt. 0)then
                do jj=1,nc(j,jk)
                ll=icon(jj,j,jk)
	        ils=ior(jj,j,jk)		
                ll=invmap(ll,jk,1)
        score(i,j)=score(i,j)+an(i,ll,ils)
                end do
        end if
        ENDDO
	ENDDO
	lst=2	
c	GAP_OPEN=-10.5
c	GAP_EXTN=-0.05
      	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(2)=result-e(2)
c Prepare a SCORE table for the original profile
	DO I=0,NSEQ1
	DO J=0,NSEQ2
	SCORE(I,J)=0.
	END DO
	END DO
c        GAP_OPEN=-10.47
c       GAP_EXTN=-0.12
c       AWT=1.52
c       BWT=1.33
	GAP_OPEN=-10.
	GAP_EXTN=-0.1
	AWT=1.3
	BWT=0.8



        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	ss=0.
	do jres=0,19
	ss=ss+amut2(jres,i)*xm(jres,j,jk)
	end do
	do jres=20,39
	ss=ss+amut(jres,i)*(xm(jres+20,j,jk))
	end do

 	score(i,j)=ss+eseq(ibinseq,ibin)*bwt+awt
	end do
	end do


	lst=3	
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(3)=result-e(3)
c	GAP_OPEN=-14.
c	GAP_EXTN=-0.7
c	AWT=1.1
c	BWT=0.6

	GAP_OPEN=-14.5
	GAP_EXTN=-0.75
c	AWT=1.1

	BWT=0.2
	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
 	score(i,j)=score(i,j)+eseq(ibinseq,ibin)*bwt
        if(nn.gt. 0)then
                do jj=1,nc(j,jk)
                ll=icon(jj,j,jk)
	        ils=ior(jj,j,jk)				
                ll=invmap(ll,jk,3)
        score(i,j)=score(i,j)+an2(i,ll,ils)
                end do
        end if
        ENDDO
	ENDDO
	lst=4		
      	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(4)=result-e(4)
	do k=1,4
	energ(jk,k)=e(k)
	en(k)=en(k)+e(k)
	en2s(k)=en2s(k)+e(k)**2
	end do
	END IF
1000	continue
	write(9,*)'number of cases= ',icount
	Aden=icount


	do k=1,4
	en(k)=en(k)/aden
	en2s(k)=en2s(k)/aden
	sd(k)=sqrt(en2s(k)-en(k)**2)
64	format(i2,1x,2(1f7.3,1x))
	do jk=1,ndata
	item(jk,k)=jk
	end do
	do jk=1,ndata
	do kk=jk+1,ndata
	if(energ(jk,k).lt.energ(kk,k))then
	itt=item(jk,k)
	item(jk,k)=item(kk,k)
	item(kk,k)=itt
	etem=energ(jk,k)
	energ(jk,k)=energ(kk,k)
	energ(kk,k)=etem
	end if
	end do
	end do

	imv=0
	do jk=1,30
	jjm=item(jk,k)
	zz=(energ(jk,k)-en(k))/sd(k)
        if(zz .gt. 1.3 )then
	imv=imv+1
	if(imv .gt.5)go to 1111
	tpsec(jk,k)=.true.
c       attempt to extractcontacts
	
	nn=0
        do j=1,nres(jjm)
        if(invmap(j,jjm,k).gt.-1)then
	nn=nn+1
        i=invmap(j,jjm,k)
                do jj=1,nc(j,jjm)
                ll=icon(jj,j,jjm)
                ll=invmap(ll,jjm,k)
        if(ll .gt. -1)conpred(i,ll)=conpred(i,ll)+1
                end do
        end if
        end do
	write(9,*)jjm,' ',name(jjm),energ(jk,k),zz,nn

        end if
1111	continue	
	end do
	write(9,*)'++++++++++'

	END DO
	do j=1,mres
	do i=1,mres
	conpred3(i,j)=0
	end do
	end do
	do k=1,4
	do jk=1,30
	if(.not.tpsec(jk,k))go to 6543
	jjm=item(jk,k)

	do i=-4,mres+4
	do j=i,mres+4
	touch2(i,j)=.false.
	touch2(j,i)=.false.
	end do
	end do
        do j=1,nres(jjm)
        if(invmap(j,jjm,k).gt.-1)then
	i=invmap(j,jjm,k)
	ires=jcode(i)
	jres=icode(j,jjm)
	if(imut(ires,jres).gt.0)then
                do jj=1,nc(j,jjm)
                ll=icon(jj,j,jjm)
                li=invmap(ll,jjm,k)
		if(li .gt. -1)then
		ires2=jcode(li)
		jres2=icode(ll,jjm)
		if(imut(ires2,jres2).gt.0)then
		touch2(i,li)=.true.
		touch2(li,i)=.true.		
		end if
		end if
		end do
	end if
	end if
	end do
	do i=1,mres
	do j=i+3,mres
	amx=conpred(i,j)
	if(amx.gt.conpred(j,i))then
	conpred(j,i)=amx
	else
	conpred(i,j)=conpred(j,i)
	end if
	end do
	end do
	do j=1,mres
	do i=1,mres
	if(conpred(i,j) .gt.3)then
	conpred3(i,j)=conpred3(i,j)+1
	if(touch2(i,j))then
	do id=-3,3
	iid=i+id
	do jd=-3,3	
	jjd=j+jd
		if(touch2(iid,jjd))then
		conpred3(iid,jjd)=conpred3(iid,jjd)+1
		end if
	end do
	end do
	end if
	end if
	end do
	end do
6543	continue	
	END DO
	END DO

	nt=0
	do j=1,mres	
	do i=1,mres
	if(iabs(i-j).gt.4)then
	if(conpred3(i,j).gt.3)then
	nt=nt+1
	end if
	end if
	end do
	end do
	
	if(nt .gt.0)then
	ares=mres
	aden=float(nt)/(ares**2)
	end if
	do j=1,mres	
	jres=jcode(j)
	do i=1,mres
	ires=jcode(i)
	if(iabs(i-j).gt.4)then
	if(conpred3(i,j).gt.3)then
	do ils=1,3
 	an2(i,j,ils)=-alog(conpred3(i,j)/aden)
	end do
	end if
	end if
	end do
	end do
	write(9,*)' number of contacts predict=',nt
	write(9,*)'============================='
	
	OPEN(unit=21, file=home(1:nhome)//genome//'zpot4a/'//namep(1:nsize)//'.predictedraorienrev')
	write(21,*)nt/2
	do j=1,mres	
	do i=1,mres
	if(iabs(i-j).gt.4)then
	if(conpred3(i,j).gt.4)then
	do ils=1,3
 	an2(i,j,ils)=-alog(conpred3(i,j)/aden)
	end do
	if(i .lt.j)then
	write(21,*)i,j
c871 	format(a1,3x,a1,2x,i4,2x,i4)
	end if
	end if
	end if
	end do
	end do
	close(21)

 	OPEN(unit=44,file=home(1:nhome)//genome//'zpot4a/'//namep(1:nsize)//'.potcontactexpr5saorienrev')
	do ils=1,3
        do i=1,mres
        iseq=jcode(i)
        write(44,*)i,aa(iseq)
        WRITE(44,167) (an2(i,j,ils),j=1,mres)
 167     FORMAT(25F5.1)
        write(44,*)'================='
	end do
	end do


        write(69,9991)ik,namep(1:nsize)
9991    format(i4,2x,a)	

1001	continue
	rewind(33)
	CLOSE(33)
	RETURN
	END 	

C
      subroutine read_msa(nunit)
c
c------------------------------------------------------------------------

c------------------------------------------------------------------------
       parameter     (maxseq =300)
       parameter     (nmax=1000)
      character     seq*1
      character*1   aa1(0:19)
      dimension seq(maxseq,2000)
      common/seq/nresseq,nseq_ori,kcode(nmax,maxseq)
      DATA AA1/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/

	read(nunit,*)nseq_ori,nresseq
	if(nseq_ori .gt.maxseq)nseq_ori=maxseq
	do i=1,nseq_ori
	read(nunit,1)(seq(i,j),j=1,nresseq)
1	format(50a1)
	end do
	if(nresseq.gt.nmax)nresseq=nmax
c	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	do i=1,nseq_ori
	do j=1,nresseq
		kcode(j,i)=-1
		do ires=0,19
			if(aa1(ires) .eq. seq(i,j))then
			kcode(j,i)=ires
			go to 18
			end if
			end do	
18 	continue
	end do
	end do
c	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

99	 return
      	end
c       *********************************************************!
c
	SUBROUTINE DP(NSEQ1,NSEQ2,RESULT,ils,ilen)
	PARAMETER(nmax=1000)
	PARAMETER(ncases=8700)
        parameter (maxres=500)       
        COMMON /ALIGNPARAM/GAP_OPEN, GAP_EXTN,SCORE(0:MAXRES,0:NMAX)
	common/mappings2/jk,lst
	common/alignment/invmap(0:nmax,0:ncases,4)		
	double precision prev, preh, pred, val,D,H,V,best
        dimension idir(0:maxres,0:nmax),val(0:maxres,0:nmax)
        dimension preV(0:maxres,0:nmax),preH(0:maxres,0:nmax),preD(0:maxres,0:nmax)
        dimension idirH(0:maxres,0:nmax),idirV(0:maxres,0:nmax)
	logical*1 tval(0:maxres,0:nmax)

c	DO DYNAMIC PROGRAMMING FIRST
c	===========================================================
       val(0,0)=0.0
       do i=1,nseq1
	tval(i,0)=.false.
         val(i,0)=0
         idir(i,0)=0
         preD(i,0)=0.0
         preH(i,0)=-1000.0
         preV(i,0)=-1000.0
      enddo
      do j=1,nseq2
	invmap(j,jk,lst)=-1
       	tval(0,j)=.false.
         val(0,j)=0
         idir(0,j)=0
         preD(0,j)=0.0
         preH(0,j)=-1000.0
         preV(0,j)=-1000.0
      enddo
         do 222 j=1,nseq2
      do 111 i=1,nseq1
	tval(i,J)=.false.
            preD(i,j)=val(i-1,j-1)+score(i,j)
            D=preD(i-1,j)+gap_open
            H=preH(i-1,j)+gap_extn
            V=preV(i-1,j)+gap_extn
            if(D.gt.H.and.D.gt.V)then
               preH(i,j)=D
               idirH(i-1,j)=1
            elseif(H.gt.V)then
               preH(i,j)=H
               idirH(i-1,j)=2
            else
               preH(i,j)=V
               idirH(i-1,j)=3
            endif
            D=preD(i,j-1)+gap_open
            H=preH(i,j-1)+gap_extn
            V=preV(i,j-1)+gap_extn
            if(D.gt.H.and.D.gt.V)then
               preV(i,j)=D
               idirV(i,j-1)=1
            elseif(H.gt.V)then
               preV(i,j)=H
               idirV(i,j-1)=2
            else
               preV(i,j)=V
               idirV(i,j-1)=3
            endif
            
            if(preD(i,j).gt.preH(i,j).and.preD(i,j).gt.preV(i,j))then
               idir(i,j)=1
               val(i,j)=preD(i,j)
            elseif(preH(i,j).gt.preV(i,j))then
               idir(i,j)=2
               val(i,j)=preH(i,j)
            else
               idir(i,j)=3
               val(i,j)=preV(i,j)
            endif
		IF(val(i,j).gt.0.)then
		tval(i,j)=.true.
		else
		val(i,j)=0.
		end if
		if(preH(i,j).lt.0.)preH(i,j)=0.
		if(preV(i,j).lt.0.)preV(i,j)=0.

 111  continue
 222     continue


c	===========================================================


	ib=0
	jb=0
	best=0.
	do j=1,nseq2
	do i=1,nseq1
	if(tval(i,j))then
	if(val(i,j).gt.best)then
	ib=i
	jb=j
	best=val(i,j)
	end if
	end if
	end do
	end do
	i=ib
	j=jb
        DO WHILE ((i.GT.0).AND.(j.GT.0))
	IF(IDIR(I,J).EQ.1)THEN
  		if(tval(i,j))then
		invmap(j,jk,lst)=i
        	     i=i-1
	             j=j-1
                 else
			go to 10	
		end if
        ELSEif(idir(i,j).eq.2)then
		if(tval(i-1,j))then
                i=i-1
  	       idir(i,j)=idirH(i,j)
		else
		go to 10
		end if
	ELSE 	
		if(tval(i,j-1))then
	        j=j-1
              idir(i,j)=idirV(i,j)
		else
		go to 10
		end if
             ENDIF
         ENDDO
10	result=best
	return
        end
