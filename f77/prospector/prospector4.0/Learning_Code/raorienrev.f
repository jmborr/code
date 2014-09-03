c	like alignb.f but uses specific penalties
c	uses r14 type of potential
c	match
	PARAMETER(nmax=1000)
        parameter (maxseq=300)
        parameter (maxres=500)
	PARAMETER(ncases=5700)
	CHARACTER*5 NAME(0:ncases),name2,name5
	character*20 a20
	character*1 ad
	character*3 amino
	LOGICAL*1 DIR
	DIMENSION ICODE(0:NMAX,0:NCASES),nres(0:ncases),x3(3)
	common/profile/xm(0:19,0:nmax,0:ncases)
	common/c/ndata,ndata2
	COMMON/S19/jbin(nmax+1,0:ncases)
	common/cont/nc(0:nmax,0:ncases),icon(15,0:nmax,0:ncases)
	common/cont2/ior(15,0:nmax,0:ncases)
	COMMON/SEQALL/imut(-1:19,-1:19)	
         COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX), 
     *    GAP_OPEN, GAP_EXTN
	common/dir2/DIR(0:MAXRES,0:NMAX)
c	=================================
c	directories
	common/namesize/nsize		
	character*255 home,templatedat
	common/home/nhome,home
c	=================================
c       PSIBLAST
        logical*1 tncbi(0:22)
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
	OPEN(UNIT=20,FILE='LIST.jan2005good35')
	READ(20,*)NDATA
	NTOT=0
	DO 1100 IK=1,NDATA
	READ(20,1)NAME(ik)
1	format(a5)
	OPEN(unit=30,file='/library/orien6/'//NAME(ik)//'.SEQ')
	rewind(30)
	read(30,2)name2,nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
        read(30,*)(icode(i,ik),i=1,nres(ik))
2	format(a5,1x,i5)
	do i=1,nres(ik)
	if(icode(i,ik).gt.19)icode(i,ik)=0.
	end do
	rewind(30)
	close(30)
	
	OPEN(unit=71,file='/library/orien6fit/'//NAME(ik)//'.FITCONT_OR')
	rewind(71)
	do i=1,nres(ik)
	read(71,*)ii,nc(ii,ik),(icon(j,ii,ik),j=1,nc(ii,ik))
	read(71,*)ii,nc(ii,ik),(ior(j,ii,ik),j=1,nc(ii,ik))
        do j=1,nc(ii,ik)
        if(icon(j,ii,ik).gt. nmax)icon(j,ii,ik)=0
        end do

	end do
	close(71)
	OPEN(UNIT=14,FILE=templatedat(1:nadr)//'pdb_jan05/seq/'//name(ik)//'.SEQ')
	do i=1,nres(ik)
	read(14,*)id,amino,iseq
	if(iseq.eq.2)then
	jbin(i,ik)=3
	elseif(iseq.eq.4)then
	jbin(i,ik)=2
	elseif(iseq.eq.3)then
	jbin(i,ik)=1
	else
	jbin(i,ik)=0
	end if
	end do
	rewind(14)
	close(14)
1100	continue
	CLOSE(20)
	DO 1101 IK=1,NDATA
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
1101	continue

	CALL r14STAT(name,NRES,ICODE)
	
c	END DO

	STOP
	END


	
	SUBROUTINE r14STAT(name,NRES,ICODE)
	PARAMETER(nmax=1000)
	PARAMETER(ncases=5700)
	parameter     (maxseq=300)
        parameter     (maxres=500)
	LOGICAL*1 HOMMAP(0:ncases)
	LOGICAL*1 DIR
	character*5 name(0:ncases),name2,nameselect
	character*255 namep

	character*4 genome
	character*6 prodin	
	character*1 seq1,seq2,aa(0:19)
	character*3 amino
	character*21 text
	dimension jseq(0:maxres),eseq(0:3,0:3)
	character*255 home
	common/home/nhome,home
	CHARACTER*3 ad3
	COMMON/S19/jbin(nmax+1,0:nCASES)
	COMMON/C/NDATA,NDATA2
        COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX), 
     *    GAP_OPEN, GAP_EXTN
	common/dir2/DIR(0:MAXRES,0:NMAX)    
	COMMON/SEQALL/imut(-1:19,-1:19)	
	DIMENSION ICODE(0:NMAX,0:NCASES)
	DIMENSION NRES(0:ncaseS)
	common/alignment/invmap(0:nmax,0:ncases,4)	
	logical*1 touch2(-4:nmax,-4:nmax)
	logical*1 tpsec(30,4)
	dimension conpred3(0:nmax,0:nmax)
	dimension en(-2:7),en2s(-2:7),sd(-2:7)
	dimension e(-2:7)
	common/profile/xm(0:19,0:nmax,0:ncases)
	common/mappings/jcode(-1:maxres)
	common/mappings2/jk,lst
	common/cont/nc(0:nmax,0:ncases),icon(15,0:nmax,0:ncases)
	common/cont2/ior(15,0:nmax,0:ncases)
	dimension apabl(-1:19,-1:19)
	dimension energ(ncases,4),item(ncases,4)
	dimension an(-1:maxres,-1:maxres,3)
	dimension amut(-1:19,0:maxres),amut2(-1:19,0:maxres)
	dimension an2(-1:maxres,-1:maxres,3),conpred(-1:maxres,-1:maxres)
	dimension an0(-1:maxres,-1:maxres)
	dimension an01(-1:maxres,-1:maxres,3)
        common/seq/nresseq,nseq_ori,kcode(maxres,maxseq)
	LOGICAL*1 HOMOLOGY(0:NCASES)
	character*255 f3590,f3590aln,fe10,fe10aln,pred,predSEQ
	common/namesize/nsize		
	common/gv/go,gv
	DATA AA/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/

	text='ATOM     1   CA  GLY '          
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

	VAL(0,0)=0.
	dir(0,0)=.false.
	do i=1,nmax
	dir(0,i)=.true.	
	end do
	do i=1,maxres
	dir(i,0)=.true.
	end do
	OPEN(UNIT=27,file=home(1:nhome)//'strlist/SCALEAV2',status='old')
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
	OPEN(UNIT=9,FILE='stat_raorienrev'//prodin)
	OPEN(UNIT=69,FILE='completedraorienrev'//prodin)	
	OPEN(unit=33,file='LIST.targ'//prodin)
	read(33,*)ndata2,genome
	do 1001 Ik=1,NDATA2
	read(33,11)namep(1:nsize)
11	format(a)
	do jk=1,ndata
	homology(jk)=.true.
	end do
        open(unit=8,file=home(1:nhome)//'/pdbhomoljan05/'//namep(1:nsize)//'.homol')
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
	call read_msa(nunit)
	rewind(31)	
	close(31)
	mres=nresseq
	OPEN(UNIT=14,FILE=pred(1:npred)//namep(1:nsize)//predseq(1:npredseq))
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
	rewind(14)
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
        OPEN(unit=44,file=home(1:nhome)//genome//'PAIR/'//namep(1:nsize)//
     &	'.PAIRPSI')
       do i=1,nresseq
       read(44,*)
       read(44,167) (an0(i,j),j=1,nresseq)
       read(44,*)
       end do
       close(44)
      OPEN(unit=44,file=home(1:nhome)//genome//'PAIROR/'//namep(1:nsize)//
     &  '.PAIRORPSI')
        do ils=1,3
       do i=1,nresseq
       read(44,*)
       read(44,167) (an01(i,j,ils),j=1,nresseq)
       read(44,*)
       end do
        end do
       close(44)

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

       
	nold=nseq_ori
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
c	an01(i,j,ils)=an0(i,j)
        an01(i,j,ils)=(an0(i,j)+0.25*an01(i,j,ils))/1.25

	an(i,j,ils)=0.55*(an(i,j,ils)-an01(i,j,ils))
	end do
	end do
	end do

	OPEN(UNIT=31,FILE=f3590(1:n3590)//namep(1:nsize)//f3590aln(1:naln90))
	call read_msa(nunit)
	rewind(31)	
	close(31)


	do i=-1,mres
	do j=-1,mres
	do ils=1,3
	an2(j,i,ils)=0.
	end do
	end do
	end do
	do j=1,mres
	jres=kcode(j,nseq_ori+1)
	do i=1,mres
	do ic=1,nseq_ori
	ires=kcode(i,ic)
	do ils=1,3
	an2(i,j,ils)=an2(i,j,ils)-apabl(ires,jres)
	end do
	end do
	do ils=1,3
	an2(i,j,ils)=5.*an2(i,j,ils)/nseq_ori
	an2(i,j,ils)=0.6*(an2(i,j,ils)-an01(i,j,ils))
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
	nden=nden+1
		if(ires.gt.-1)then
		amut2(ires,i)=amut2(ires,i)+1
		end if
	end do	
		do iseq=0,19
		amut2(iseq,i)=amut2(iseq,i)/nden
		end do
	end do


	

	close (15)
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
	
	icount=icount+1
c	structure is jk
c	sequence is ik

	ils=0
	NSEQ2=nres(jk)
	DO J=0,NSEQ2	
	DO I=0,NSEQ1
	SCORE(I,J)=0.
	END DO
	END DO
c Prepare a SCORE table
        GAP_OPEN=-5.47
        GAP_EXTN=-0.16
        AWT=1.52
        BWT=0.98
        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1
	ss=0.
	do jres=0,19
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
	GAP_OPEN=-10.5
	GAP_EXTN=-0.05
	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1
 	score(ii,j)=score(ii,j)+eseq(ibinseq,ibin)*1.4
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
        GAP_OPEN=-10.47
        GAP_EXTN=-0.12
        AWT=1.52
        BWT=1.33

        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1	
	ss=0.
	do jres=0,19
	ss=ss+amut2(jres,i)*xm(jres,j,jk)
	end do
	score(ii,j)=ss+eseq(ibinseq,ibin)*bwt+awt
	end do
	end do

	lst=3	
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(3)=result		
	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	II=NSEQ1-I+1	
 	score(ii,j)=score(ii,j)+eseq(ibinseq,ibin)*1.3
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
	GAP_OPEN=-14.0
	GAP_EXTN=-0.7
	lst=4		
      	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(4)=result
        GAP_OPEN=-5.47
        GAP_EXTN=-0.16
        AWT=1.52
        BWT=0.98
        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	ss=0.
	do jres=0,19
	ss=ss+amut(jres,i)*xm(jres,j,jk)
	end do
 	score(i,j)=ss+eseq(ibinseq,ibin)*bwt+awt
	end do
	end do


C Prepare a SCORE table
	lst=1
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(1)=result-e(1)
C Prepare a SCORE table
c	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
 	score(i,j)=score(i,j)+eseq(ibinseq,ibin)*1.4
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
	GAP_OPEN=-10.5
	GAP_EXTN=-0.05
      	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(2)=result-e(2)

c Prepare a SCORE table for the original profile
	DO I=0,NSEQ1
	DO J=0,NSEQ2
	SCORE(I,J)=0.
	END DO
	END DO
        GAP_OPEN=-10.47
        GAP_EXTN=-0.12
        AWT=1.52
        BWT=1.33


        DO j=1,NSEQ2
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
	ss=0.
	do jres=0,19
	ss=ss+amut2(jres,i)*xm(jres,j,jk)
	end do
 	score(i,j)=ss+eseq(ibinseq,ibin)*bwt+awt
	end do
	end do


	lst=3	
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	e(3)=result-e(3)
	ivs=1
        DO j=1,NSEQ2
	nn=nc(j,jk)	
	ibin=jbin(j,jk)
        DO i=1,NSEQ1	
        ibinseq=jseq(i)
 	score(i,j)=score(i,j)+eseq(ibinseq,ibin)*1.3
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
	GAP_OPEN=-14.
	GAP_EXTN=-0.7
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
	do i=1,mres
	do j=1,mres
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

	do i=1,mres
	do j=1,mres
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
 	OPEN(unit=21,file=home(1:nhome)//genome//'zpot4a/'//namep(1:nsize)//'.predictedr1aorienrev')
	write(21,*)nt/2
	do i=1,mres
	ires=jcode(i)
	do j=1,mres	
	jres=jcode(j)
	if(iabs(i-j).gt.4)then
	if(conpred3(i,j).gt.3)then
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
	write(9,*)' number of contacts predict=',nt
	write(9,*)'============================='
	

 	OPEN(unit=44,file=home(1:nhome)//genome//'zpot4a/'//namep(1:nsize)//'.potcontactexpr1aorienrev')
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

C  Dynamic Programming for Sequence Alignment, Fortran version
C
C
	SUBROUTINE DP(NSEQ1,NSEQ2,RESULT,ils,ilen)
C
C  Dynamic Programming for Sequence Alignment, Fortran version
C
	PARAMETER(nmax=1000)
	PARAMETER(ncases=5700)
        parameter (maxres=500)       
	integer pos1,pos2
	CHARACTER*1 seq1,SEQ2
	LOGICAL*1 DIR
       COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *    GAP_OPEN, GAP_EXTN
	common/dir2/DIR(0:MAXRES,0:NMAX)         
	common/mappings/jcode(-1:maxres)
	common/mappings2/jk,lst
	common/cont/nc(0:nmax,0:ncases),icon(15,0:nmax,0:ncases)
	common/alignment/invmap(0:nmax,0:ncases,4)		
c       WRITE(*,*)'Aligning ',SEQ1,' to ',SEQ2
       result=ALIGN(NSEQ1,NSEQ2)

c       WRITE(9,*)'Alignment score: ',result
C Extract alignment
	do j=1,nseq2
	invmap(j,jk,lst)=-1
	end do
	ib=0
	jb=0
         pos=0
	best=0.
	do i=1,nseq1
	do j=1,nseq2
	if(val(i,j).gt.best)then
	ib=i
	jb=j
	best=val(i,j)
	end if
	end do
	end do
	i=ib
	j=jb
         DO WHILE ((i.GT.0).AND.(j.GT.0))
           IF (.not.DIR(i,j)) THEN
  		if(val(i,j).gt.0)then
		invmap(j,jk,lst)=i
        	     i=i-1
	             j=j-1
                 else
			go to 10	
		end if
           ELSE
             IF (VAL(i-1,j).GT.VAL(i,j-1)) THEN
		if(val(i-1,j).gt.0)then
                i=i-1
		else
		go to 10
		end if
             ELSE
		if(val(i,j-1).gt.0)then
               j=j-1
		else
		go to 10
		end if
             ENDIF
           ENDIF
         ENDDO
10	result=best
	return
       END

C       
C  Main DP routine
C  The SCORE array should be filled in
C
       FUNCTION ALIGN(NSEQ1,NSEQ2)
	PARAMETER(nmax=1000)
	PARAMETER(ncases=5700)
        parameter (maxres=500)      
	LOGICAL*1 DIR
        COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *    GAP_OPEN, GAP_EXTN
	common/dir2/DIR(0:MAXRES,0:NMAX)         
	common/mappings2/jk,lst
	common/mappings/jcode(-1:maxres)
	common/cont/nc(0:nmax,0:ncases),icon(15,0:nmax,0:ncases)
	common/cont2/ior(15,0:nmax,0:ncases)
	common/gv/go,gv
       REAL H, V
C Init matrices

        do j=1,nseq2
         do i=1,nseq1	 
	dir(i,j)=.false.
	val(i,j)=0.
	end do
	end do
        DO j=1,NSEQ2	
        DO i=1,NSEQ1
             H=VAL(i-1,j)
             if (DIR(i-1,j)) THEN
               H=H+GAP_EXTN
             ELSE
               H=H+GAP_OPEN
	     END IF
             V=VAL(i,j-1)
             if (DIR(i,j-1)) THEN
               V=V+GAP_EXTN
             ELSE
               V=V+GAP_OPEN
             ENDIF
             D=VAL(i-1,j-1)+SCORE(i,j)
             IF ((D.GE.H).AND.(D.GE.V)) THEN
               VAL(i,j)=D
               DIR(i,j)=.false.
             ELSE
               IF ((H.GE.D).AND.(H.GE.V)) THEN
                 VAL(i,j)=H
                 DIR(i,j)=.true.
               ELSE
                 VAL(i,j)=V
                 DIR(i,j)=.true.
               ENDIF
             ENDIF
		IF(val(i,j).lt.0)then
		val(i,j)=0.
		dir(i,j)=.false.
		end if
           ENDDO
	END DO
c	update the score matrix for residues above
c	residue j is updated
	align=VAL(NSEQ1,NSEQ2)
	return
	end 

c
      subroutine read_msa(nunit)
c
c------------------------------------------------------------------------

c------------------------------------------------------------------------
       parameter     (maxseq =300)
       parameter     (maxres=500)
      parameter	    (nresmax=500)
      character     seq*1
      character*1   aa1(0:19)
      dimension seq(maxseq,maxres)
      common/seq/nresseq,nseq_ori,kcode(nresmax,maxseq)
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


