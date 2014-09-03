	PARAMETER(NMAX=1700)
	PARAMETER(MAXSEQ=300)
	PARAMETER(MAXRES=500)
	PARAMETER(NCASES=1000)

	LOGICAL*1 TNCBI(0:22)

	CHARACTER*1 AD
	CHARACTER*1 AA(0:19),ANCBI(1:22)
	CHARACTER*3 AMINO	
	CHARACTER*4 GENOME
	CHARACTER*5 NAME5,NAME(0:NCASES),NAME2
	CHARACTER*6 PRODIN
	CHARACTER*20 HOME,A20
	CHARACTER*255 NAMEP
	
	DIMENSION ICODE(0:NMAX,0:NCASES),nres(0:ncases)
	DIMENSION EN(4),SD(4)
	DIMENSION IX(0:22),NCBI(22)

	COMMON/HOME/HOME	
	COMMON/PROFILE/XM(0:19,0:NMAX,0:NCASES)
	COMMON/C/NDATA,NDATA2
	COMMON/S19/JBIN(NMAX+1,0:NCASES)
	COMMON/CONT/NC(0:NMAX,0:NCASES),ICON(15,0:NMAX,0:NCASES)
	COMMON/CONT2/IOR(15,0:NMAX,0:NCASES)
	COMMON/NAMESIZE/NSIZE		
	COMMON/PAIRS/APABL(-1:19,-1:19,3)
	COMMON/BACKBONE/XA(3,NMAX,NCASES)
c	=================================
c       PSIBLAST

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
        read(5,13)prodin
13      format(a6)

	open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)

	home='/nfs/users/skolnick/'		
	
	OPEN(UNIT=9,FILE='stat_rap4orienrev'//prodin)
	OPEN(UNIT=69,FILE='completedrap4orienrev'//prodin)	
	OPEN(unit=33,file='LIST.targ'//prodin) ! note the conversion from STRUCTURE_rap3orienrev in
c	pap4orienrev now we split it into chunks
	
	read(33,*)ndata2,genome
        DO 4000 IIK=1,NDATA2
        READ(33,9432)namep(1:nsize)
9432    format(a,3x,a5,1f7.3,3x,a60)

        OPEN(unit=3,file=namep(1:nsize)//'.threadrap4orienrev')
        OPEN(unit=34,file=home//genome//'domainjan05/'//namep(1:nsize)//'rap4orienrev.pdb')
        OPEN(unit=53,file=namep(1:nsize)//'.threadrap3orienrev')	 
	read(53,93,END=4000, ERR=4000)nth,(en(k),sd(k),k=1,4)
	if(nth.gt.5)nth=5
93	format(i5, 8(1f9.3,1x))
	write(3,*)nth
	do 3000 lth=1,nth
        read(53,18)name5,zd,ksel
18      format(a5,1x,1f7.3,1x,i4)
	write(9,*)'==========================='
	write(9,*)namep(1:nsize),' name of family member: ',name5
	write(9,*)'==========================='
c        open(unit=20,file='/nfs/users/piotr/pdblist/jun2004/clusters/'//name5//'.cluster')
c        read(20,*)ndata
        OPEN(unit=27,file=home//'domainsjan05/'//name5//'/'//NAME5//'.SCALEFAMOR')
	rewind(27)
	do it=1,3
	read(27,*)
	read(27,*)
	DO I=0,19
	read(27,17) Ad3,(apabl(I,J,it),J=0,19)
	end do
	end do
17	FORMAT(A3,23F5.1)
	close(27)

	NDATA=1
	IF(NDATA .GT. 100)NDATA=100
	DO 1100 IK=1,NDATA
c	READ(20,1)NAME(ik)
	name(ik)=NAME5
1	format(a5)

        OPEN(unit=30,file=home//'domainsjan05/'//name5//'/seq/'//NAME(ik)//'.SEQ')
	rewind(30)
	read(30,2)name2,nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
2	format(a5,1x,i5)
	read(30,*)(icode(i,ik),i=1,nres(ik))
 	do i=1,nres(ik)
 	if(icode(i,ik).gt.19)icode(i,ik)=0.
 	end do
	rewind(30)
	close(30)

        OPEN(unit=71,file=home//'domainsjan05/'//name5//'/fit/'//NAME(ik)//'.FITCONT_OR')
	do i=1,nres(ik)
	read(71,*)ii,nc(ii,ik),(icon(j,ii,ik),j=1,nc(ii,ik))
	read(71,*)ii,nc(ii,ik),(ior(j,ii,ik),j=1,nc(ii,ik))
        do j=1,nc(ii,ik)
        if(icon(j,ii,ik).gt. nmax)icon(j,ii,ik)=0
        end do
	end do
	close(71)

	OPEN(unit=45,file=home//'domainsjan05/'//name5//'/real/'//NAME(ik)//'.real')
	rewind(45)
	read(45,*)nres(ik)
	if(nres(ik).gt.nmax)nres(ik)=nmax
	do i=1,nres(ik)
	read(45,*)id,(xa(j,i,ik),j=1,3)
	end do
	close(45)
	
	
! 	fix this input file for the most general case.
c   OPEN(UNIT=14,FILE='/nfs/users/adrian/dat/pdblib/'//name5//'/'//name(ik)//'/'//name(ik)//'.SEQ')
	OPEN(UNIT=14,FILE='/nfs/users/adrian/dat/pdb_jan05/seq/'//name(ik)//'.SEQ')
        do i=1,nres(ik)
        read(14,*)id,amino,iseq
        if(iseq.eq.2)then
        jbin(i,ik)=3
        elseif(iseq.eq.4)then
        jbin(i,ik)=2
        elseif(iseq .eq.3)then
        jbin(i,ik)=1
        else
        jbin(i,ik)=0
        end if
        end do
        rewind(14)
        close(14)
	do i=1,nres(ik)
	do ires=0,19
	xm(ires,i,ik)=0.
	end do
	end do

! fix this..... input file
        open(unit=31,file='/nfs/restore/DDPSC2/adrian_dat/pdb_jan05/profile/'//name(ik)//'.mtx')
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

1100	continue
	close(20)
        CALL r14STAT(name,NRES,genome,namep,name5)
3000	CONTINUE	
121        write(69,9991)ik,namep(1:nsize)
9991    format(i4,2x,a)	

4000	CONTINUE

	open(unit=1,file='finish_rap4orienrev'//prodin)
	rewind(1)
	write(1,*)'1'

	STOP
	END


	
	SUBROUTINE r14STAT(name,NRES,genome,namep,name5)

	PARAMETER(NMAX=1700)
	PARAMETER(NCASES=1000)
	PARAMETER(MAXSEQ=300)
	PARAMETER (MAXRES=500)
	PARAMETER (ZMAX3=4.)
	LOGICAL*1 THREADG(0:NCASES)
	LOGICAL*1 DIR,homology(ncases)

	CHARACTER*1 AA(0:19)
	CHARACTER*3 AD3,AMINO,AA3(0:19)
	CHARACTER*4 GENOME	
	CHARACTER*5 NAME(0:NCASES)
	CHARACTER*6 PRODIN
	CHARACTER*17 TEXT
	CHARACTER*20 HOME
	CHARACTER*255 NAMEP
	CHARACTER*255 F3590,F3590ALN,FE10,FE10ALN,PRED,PREDSEQ

	DIMENSION ICODE(0:NMAX,0:NCASES)
	DIMENSION NRES(0:NCASES),ZU(0:NCASES),JUD(0:NCASES)
	DIMENSION ZMAX(0:NCASES,4)
	DIMENSION JSEQ(0:MAXRES),ESEQ(0:3,0:3)
	DIMENSION EN(0:4),EN2S(0:4),SD(0:4)
	DIMENSION E(0:4)
	DIMENSION ENERG(NCASES,4),ITEM(NCASES,4)
	DIMENSION AN(-1:MAXRES,-1:MAXRES,3)
	DIMENSION AMUT(-1:19,0:MAXRES),AMUT2(-1:19,0:MAXRES)
	DIMENSION AN2(-1:MAXRES,-1:MAXRES,3)
	DIMENSION AN01(-1:MAXRES,-1:MAXRES,3),AN02(-1:MAXRES,3)
	DIMENSION KMX(0:NCASES),JCODE(-1:MAXRES)

	COMMON/HOME/HOME
	COMMON/S19/JBIN(NMAX+1,0:NCASES)
	COMMON/C/NDATA,NDATA2
	COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX), 
     *    GAP_OPEN, GAP_EXTN
	COMMON/DIR2/DIR(0:MAXRES,0:NMAX)    
	COMMON/ALIGNMENT/INVMAP(0:NMAX,0:NCASES,4)	
	COMMON/PROFILE/XM(0:19,0:NMAX,0:NCASES)
	COMMON/BACKBONE/XA(3,NMAX,NCASES)
	COMMON/MAPPINGS2/JK,LST
	COMMON/CONT/NC(0:NMAX,0:NCASES),ICON(15,0:NMAX,0:NCASES)
	COMMON/CONT2/IOR(15,0:NMAX,0:NCASES)
	COMMON/PAIRS/ APABL(-1:19,-1:19,3)
	
	COMMON/SEQ/NRESSEQ,NSEQ_ORI,KCODE(MAXRES,MAXSEQ)
	COMMON/NAMESIZE/NSIZE		
	     DATA AA3/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &               'MET','ASP','ASN','LEU',
     &               'LYS','GLU','GLN','ARG',
     &               'HIS','PHE','TYR','TRP'/




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

	
        text='ATOM     1   CA  '	

	VAL(0,0)=0.
	dir(0,0)=.false.
	do i=1,nmax
	dir(0,i)=.true.	
	end do
	do i=1,maxres
	dir(i,0)=.true.
	end do

	do 1001 Ik=1,NDATA
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
	elseif(iseq.eq.4)then
	jseq(i)=2
	elseif(iseq .eq.3)then
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
        eseq(2,2)=1.
        eseq(3,3)=1.
        eseq(2,3)=-1.
        eseq(3,2)=-1.

	do j=-1,mres
	do i=-1,mres
	do ils=1,3
	an(j,i,ils)=0.
	end do
	end do
	end do


	OPEN(unit=44,file=home//genome//'zpot4a/'//namep(1:nsize)//'.potcontactexpr1aorienrev')
       do ils=1,3
       do i=1,mres
       read(44,*)
       read(44,*) (an01(i,j,ils),j=1,mres)
       read(44,*)
       end do
       end do
       close(44)

	OPEN(unit=44,file=home//genome//'zpot4a/'//namep(1:nsize)//'.potcontactexpr1ap2orienrev')
       do ils=1,3
       do i=1,mres
       read(44,*)
       read(44,*) (an02(j,ils),j=1,mres)
	do j=1,mres
        an01(i,j,ils)=(an02(j,ils)+0.25*an01(i,j,ils))/1.25
	end do
       read(44,*)
       end do
       end do
       close(44)


	do j=1,mres
	do i=1,mres
	do ic=1,nseq_ori
	ires=kcode(i,ic)
	jres=kcode(j,ic)
	do ils=1,3
	an(i,j,ils)=an(i,j,ils)-apabl(ires,jres,ils)
	end do
	end do
	do ils=1,3
	an(i,j,ils)=5*an(i,j,ils)/nseq_ori
c        an01(i,j,ils)=(an02(i,j,ils)+0.25*an01(i,j,ils))/1.25
	an(i,j,ils)=0.55*(an(i,j,ils)-an01(i,j,ils))
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
	an2(i,j,ils)=an2(i,j,ils)-apabl(ires,jres,ils)
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
        do jk=1,ndata
        threadg(jk)=.false.
	do k=1,4
	zmax(jk,k)=0.
	end do
        end do
	zselect=-5000
	do jk=1,ndata
	homology(jk)=.true.
	end do
	ndata=1 ! fix this later
	do 1000 jk=1,NDATA
	if(homology(jk))then
	icount=icount+1
c	structure is jk
c	sequence is ik
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
	GAP_OPEN=-14.
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
	sd(k)=sqrt(en2s(k)-en(k)**2+0.000001)

	do jk=1,ndata
	item(jk,k)=jk
	write(9,*)energ(jk,k)
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

	do jk=1,1
	jjm=item(jk,k)
c	zz=(energ(jk,k)-en(k))/sd(k)
	zz=energ(jk,k)
	zmax(jjm,k)=zz/10.
c        if(zz .gt. -30)then
c		if(zz .gt. zmax3)then

		threadg(jjm)=.true.
c		end if
	       imv=imv+1
        if(imv .gt.5)go to 1111
	write(9,87)jjm,name(jjm),energ(jk,k),zz
87 	format(i4,1x,a5,1x,2(1f9.3,1x),1x,9(a9,1x))
c	end if
	write(9,*)'++++++++++'
1111	continue
	end do
	END DO

	nn=0
	do jk=1,ndata
	if(threadg(jk))then
	nn=nn+1
	zu(nn)=0
	zm=-3000.
	kx=0
	do k=1,4
	if(zmax(jk,k).gt.zm)then
	zm=zmax(jk,k)
	kx=k
	endif
	end do
	kmx(nn)=kx
	zu(nn)=zm
	write(9,*)'nn,zu(nn)',nn,zu(nn)
	jud(nn)=jk
	end if
	end do
93	format(i5, 8(1f9.3,1x))
	
	k=1
	do jk=1,nn
	item(jk,k)=jk
	end do
	do jk=1,nn
	do kk=jk+1,nn
	if(zu(jk).lt.zu(kk))then
	itt=item(jk,k)
	item(jk,k)=item(kk,k)
	item(kk,k)=itt
	etem=zu(jk)
	zu(jk)=zu(kk)
	zu(kk)=etem
	ii=kmx(jk)
	kmx(jk)=kmx(kk)
	kmx(kk)=ii
	end if
	end do
	end do
	do jk=1,nn
c	k=kmx(jk)	
	k=4
	jjm=item(jk,1)
	iik=jud(jjm)	
	nn=0
	NSEQ2=nres(iik)
        do j=1,nseq2
        i=invmap(j,iik,k)
        if(i .gt.0)then
	nn=nn+1
	end if
	end do
	write(34,665)name(iik),nn,zu(jk),k
665	format(a5,1x,i4,1x,1f7.3,1x,i3)
	pid=0
        do j=1,nseq2
        i=invmap(j,iik,k)
        if(i .gt.0)then
	ires=jcode(i)
	jres=icode(j,iik)
	if(ires.eq.jres)pid=pid+1
        write(34,1033)text,aa3(ires),i,(xa(l,j,iik),l=1,3),j,aa3(jres)
1033	FORMAT(A17,a3,2X,I4,1x,3X,3F8.3,1x,i4,1x,a3)			
        end if
        end do
	write(34,906)'TER',pid/nn
906	format(a3,1x,1f7.3)
	write(3,18)name(iik),zu(jk),kmx(jk)
18	format(a5,1x,1f7.3,1x,i4)
	end do

	
1001	continue

	RETURN
	END 	

C  Dynamic Programming for Sequence Alignment, Fortran version
C
	SUBROUTINE DP(NSEQ1,NSEQ2,RESULT,ils,ilen)
C  Dynamic Programming for Sequence Alignment, Fortran version
C
	PARAMETER(nmax=1700)
	PARAMETER(ncases=1000)
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
	PARAMETER(nmax=1700)
	PARAMETER(ncases=1000)
        parameter (maxres=500)      
	LOGICAL*1 DIR
        COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *    GAP_OPEN, GAP_EXTN
	common/dir2/DIR(0:MAXRES,0:NMAX)         
	common/mappings2/jk,lst
	common/mappings/jcode(-1:maxres)
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

C
      subroutine read_msa(nunit)
c
c------------------------------------------------------------------------

c------------------------------------------------------------------------
       parameter     (maxseq =300)
       parameter     (maxres=500)
       character     seq*1
       character*1   aa1(0:19)
      dimension seq(maxseq,maxres)
      common/seq/nresseq,nseq_ori,kcode(maxres,maxseq)
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
