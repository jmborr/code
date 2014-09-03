c	based on c10b6c7b7.f
c	4/7/01  Jeffrey Skolnick
	PARAMETER(NMAX=3000)
	PARAMETER(NCASES=4200)
	CHARACTER*5 NAME(0:ncases),name2
	CHARACTER*3 nametarg,name3
	DIMENSION ICODE(NMAX,0:NCASES),NRES(0:NCASES)
	COMMON/nam/name3
	OPEN(UNIT=9,FILE='stat_tmeasyfrag_5')
	CALL r14STAT(name,NRES,ICODE)

	STOP
	END


	
	SUBROUTINE r14STAT(name,NRES,ICODE)
	PARAMETER(NMAX=3000)
	PARAMETER(NCASES=4200)
	PARAMETER(MAXRES=1800)
	PARAMETER(REJECT=33.)	
	CHARACTER*5 NAME(0:NCASES),NAMESELECT2(10000),NAMESELECT
	CHARACTER*5 NAMETARG
	CHARACTER*3 NAME3,TER
	character*20 a33
	character*1 a3
	CHARACTER*6 a6
	CHARACTER*60 head
	CHARACTER*5 NAMEP
	CHARACTER*1 AA(0:19),ad
	CHARACTER*21 TEXT
	CHARACTER*5 name5
	character*4 genome
	CHARACTER*6 prodin
	double precision u,t,umin(3,3),tm(3)
	common/matrix/u(3,3),t(3)
	DIMENSION XSC(3,NMAX),XSCT(3,NMAX)
	DIMENSION XSC2(3,NMAX),XSCT2(3,NMAX)	
	COMMON/C/NDATA,NDATA2
	DIMENSION XA(3,NMAX,0:NCASES)
        COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *   LIGNMENT(0:2*MAXRES), GAP_OPEN, GAP_EXTN
	COMMON/SEQALL/IMUT(-1:19,-1:19)	
	COMMON/INVERSE/INVMAP(0:MAXRES,NCASES,4)
	COMMON/MAPPINGS2/JK,LST
	COMMON/nam/name3
	COMMON/MAPPINGS/JCODE(-1:MAXRES),IMAP(0:MAXRES)     	
	dimension isz(3000)	
	DIMENSION ICODE(NMAX,0:NCASES)
	DIMENSION NRES(0:NCASES)
	DIMENSION AMUT(-1:19,0:NMAX)
	dimension ist(0:2),ist2(0:2)
	common/rotate/mres,xs(3,NMAX)
        dimension energ(0:ncases)
	dimension ires1(0:nmax),jres1(0:nmax)
	dimension rg(0:30),iav(0:30),cav(0:30)
        dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
        dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax),x3(3,nmax)
	dimension ithist(0:20),i2map(nmax),jmap(nmax)
	dimension ib(nmax),iend(nmax)
        parameter  (maxseq =3000)	
        common/seq/nresseq,nseq_ori,kcode(maxres,maxseq)	
	
	DATA AA/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/
	open(unit=55,file='nametarg1')
        read(55,13)prodin
13      format(a6)
	OPEN(unit=4,file='LIST.targ'//prodin)
	read(4,*)nd,genome

	text='ATOM     1   CA  GLY '
	OPEN(unit=20,file='/nfs/users/skolnick/strlist/BLOSUM_INT')
	read(20,*)
	do i=0,19
        read(20,*)ad,(imut(i,j),j=0,19)
	end do
	close(20)
	do i=0,19
	imut(-1,i)=-10
	imut(i,-1)=-10
	end do

	OPEN(unit=33, file='list.easy1')

	rav=0
	read(33,*)ndata2
	iden=0
	fav=0.
	favo=0.
	do imin=1,30
	rg(imin)=0.
	iav(imin)=0
	cav(imin)=0.
	end do
	
	do 1001 IIk=1,NDATA2
        read(33,9432)namep
9432    format(a5)
c	3590 % id
	write(9,*)namep,'   ****'
	OPEN(unit=45,file='/nfs/users/skolnick/cafiles4/'//NAMEp//'.real')
	rewind(45)
	read(45,*)mres
 	L2=MRES
        do i=1,L2
        read(45,*)n2(i),x2(i),y2(i),z2(i)
	end do
	close(45)
 103  format(A22,i4,A4,3F8.3)

	nunit=31
	OPEN(UNIT=31,FILE='/nfs/users/adrian/dat/old_benchmark/3590nr/'//namep//'_ini.aln')
	rewind(31)
	call read_msafasta(nunit)
	close(31)
	n3590=nseq_ori
	mres=nresseq
	nseq_ori=1
	do i=1,mres
		do iseq=-1,19
		amut(iseq,i)=0.
			do ic=1,nseq_ori
				ires=kcode(i,ic)
				amut(iseq,i)=amut(iseq,i)+imut(iseq,ires)
			end do	
		amut(iseq,i)=amut(iseq,i)/nseq_ori
		end do
	END Do
	nunit=31


	NSEQ1=MRES
	jk=1
	OPEN(unit=30,file='/nfs/users/skolnick/orien4/'//NAMEp//'.SEQ')
		rewind(30)
		read(30,2)name2,nres(jk)
		if(nres(ik).gt. nmax)nres(jk)=nmax
2		format(a5,1x,i5)
		read(30,*)(icode(i,jk),i=1,nres(jk))
	rewind(30)
	close(30)
	NSEQ2=nres(jk)
c	first pass
	DO I=0,NSEQ1
	DO J=0,NSEQ2
	SCORE(I,J)=0.
	END DO
	END DO
c Prepare a SCORE table
        DO j=1,NSEQ2
	jres=icode(j,jk)	
        DO i=1,NSEQ1
	score(i,j)=amut(jres,i)
        ENDDO
	ENDDO
C Prepare a SCORE table
	GAP_OPEN=-12.0
	GAP_EXTN=-1.2

	ks=1
	lst=1
	ils=1
	call  DP(NSEQ1,NSEQ2,RESULT,ils,length)
	pid=0.
	nn=0
	do j=1,nseq2
	if(invmap(j,jk,1).gt.-1)then
	i=invmap(j,1,1)
	ires=kcode(i,1)
	jres=icode(j,1)
	if(ires.eq.jres)pid=pid+1
	nn=nn+1
	end if
	end do
	pid=pid/nn
	write(9,*)namep,nn,'identity to pdb file=',pid
	nrank=0
        OPEN(unit=3,file='/nfs/users/skolnick/'//genome//'qzsummary/'//namep//'.thread')
        OPEN(unit=34,file='/nfs/users/skolnick/'//genome//'qztemplates/'//namep//'.pdb')
	read(3,*)nth
	if(nth.gt.5)nth=5
	rmsmin=9000
	rmsmino=9000
	TMSmax=0.
	nnmin=0
	nnmin=0
	nnmino=0
c	if(nth .gt. 1 )then
c	go to 1001
c	end if
	do 2000 ithread=1,nth
	read(34,15,END=1001, ERR=1001)name5,ic
	read(3,25)name5,zth
25	format(a5,1x,1f7.3)
15	format(a5,1x,i4,1x,1f7.3,1x,i3)
201	continue	
	nn=0
	pid=0.
********* read  structure1--------->      
	do id=1,ic
        read(34,1033)a33,i,(x3(l,id),l=1,3),j
1033    FORMAT(a20,2X,I4,1x,3X,3F8.3,1x,i4,1x,a3)
	ii=imap(i)
	i2map(id)=ii
	jmap(id)=j
	end do
	read(34,*)
c	now determine the breaks
	ib(1)=1
	isec=1
	nd=0
	do id=1,ic-1
	nd=nd+1
	if(nd.gt.15)then
	iend(isec)=id
	nd=0
	isec=isec+1
	ib(isec)=id+1
	go to 134
	end if	

	if(jmap(id)+1.ne. jmap(id+1))then
	iend(isec)=id
	if(id .lt.ic)then
	isec=isec+1
	ib(isec)=id+1
	end if
	end if
134	continue
	end do
	iend(isec)=ic
	do jsec=1,isec
	nn=0
	DO ID=IB(JSEC),IEND(JSEC)
	Ii=i2map(id)
	if(ii .gt. -1)then
		nn=nn+1			
		n1(nn)=ii
		x1(nn)=x3(1,id)
		y1(nn)=x3(2,id)
		z1(nn)=x3(3,id)
	end if
	END DO
	L1=NN
	IF(NN.GT.4)THEN	
        call TMscore(L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,TMS,Rcomm,Lcomm)
	write(9,39)namep,name5,mres,ib(jsec),iend(jsec),LCOMM, Rcomm,TMS
39	format(2(a5,1x),4(i4,1x),2(1f7.3,1x))
	imin=int(rcomm)+1
			if(imin.gt.30)imin=30
		do ir=imin,30
		rg(ir)=rg(ir)+rcomm
		iav(ir)=iav(ir)+1
		end do
	iden=iden+1
	END IF
	end do
	if(TMS .gt.TMSmax)then
	nrank=ithread
	TMSmax=TMS
	nnmin=Lcomm
	pidmin=pid
	rmsmin=RCOMM
	zmin=zth	
	end if
2000	continue
2001    continue
c	IF(TMSmax .lt.0.17)go to 1001
	nav=nav+nrank
	ntmgood=ntmgood+1
	TMav=tMav+TMSmax
20	format(a5,1x,1f7.3,1x,3(i4,1x),1f7.3,1x,i3)
22	continue
	rav=rav+rmsmin
	fav=fav+float(nnmin)/mres		
	pidav=pidav+pidmin
		ihist=20*tmsmax+1
	do iv=0,ihist
	ithist(iv)=ithist(iv)+1
	end do

	if(rmsmin .lt.6.5)then
 		ngood=ngood+1
         	ravg=ravg+rmsmin
		favg=favg+float(nnmin)/mres
	end if

	write(22,228)namep,rmsmin,nnmin
228	format(a5,1x,1f7.3,1x,i4)
1001	continue
	write(9,*)'number of cases=',iden,' <rms>=',rav/iden
	write(9,*)'number of cases=',iden,' <%Id>=',pidav/iden
	write(9,*)'number of proteins<6.5 A =',ngood
	do ihist=0,20
	write(9,50)float(ihist)/20,ithist(ihist)
50	format(1f7.4,1x,i4)
	end do
	itot=iden
	do imin=1,30
	iden=iav(imin)
	if(iden .gt.0)then
	write(9,21)'rms< ',imin,iden,rg(imin)/iden,float(iden)/itot
	end if
	end do
21	format(a5,1x,i3,1x,i6,1x,3(1f7.3,1x))

	return	
	end

	SUBROUTINE DP(NSEQ1,NSEQ2,RESULT,ils,ilen)
C
C  Dynamic Programming for Sequence Alignment, Fortran version
C
	PARAMETER(NMAX=3000)
	PARAMETER(NCASES=4200)
        PARAMETER(MAXRES=1800)

	INTEGER POS1,POS2

	LOGICAL*1 DIR
		
	COMMON/ALN2/DIR(0:MAXRES,0:NMAX)
        COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *   LIGNMENT(0:2*MAXRES), GAP_OPEN, GAP_EXTN
	COMMON/MAPPINGS/JCODE(-1:MAXRES),IMAP(0:MAXRES)     
	COMMON/MAPPINGS2/JK,KS
	COMMON/INVERSE/INVMAP(0:MAXRES,NCASES,4)

C Define penalty gaps and sequence lengths
       

C Prepare a SCORE table
      result=ALIGN(NSEQ1,NSEQ2)

C Extract alignment


       length=IEXTRACT(NSEQ1,NSEQ2)
      	ilen=0
C Write alignment
       ilipos=0
       pos1=1
       pos2=1
	do i=1,nseq1
	imap(i)=-1
	end do
	do j=1,nseq2
	invmap(j,jk,ks)=-1
	end do
       DO i=0,length-1
         IF (LIGNMENT(ilipos).EQ.0) THEN
           imap(pos1)=pos2
           invmap(pos2,jk,ks)=pos1
           ilen=ilen+1
           pos1=pos1+1
           pos2=pos2+1
         ELSEIF (LIGNMENT(ilipos).EQ.-1) THEN
           pos2=pos2+1
         ELSEIF(LIGNMENT(ilipos).EQ.1)THEN
           pos1=pos1+1
         ENDIF
         ilipos=ilipos+1	   	 
       ENDDO
	return
       END

C       
C  Main DP routine
C  The SCORE array should be filled in
C
      FUNCTION ALIGN(NSEQ1,NSEQ2)
	PARAMETER(NMAX=3000)
	PARAMETER(NCASES=4200)
	PARAMETER(MAXRES=1800)

	REAL H, V

	LOGICAL ALIGRES(NMAX)
	LOGICAL*1 DIR

	COMMON/ALN2/DIR(0:MAXRES,0:NMAX)
      COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *   LIGNMENT(0:2*MAXRES), GAP_OPEN, GAP_EXTN
	COMMON/MAPPINGS2/JK,KS
	COMMON/MAPPINGS/JCODE(-1:MAXRES),IMAP(0:MAXRES)     	
	
        
        do j=1,nseq2
        do i=1,nseq1	
	dir(i,j)=.false.
	val(i,j)=0.
	end do
	end do
	go=GAP_OPEN*0.6
	ge=GAP_EXTN*0.6

	do i=1,nseq1
	if(i.lt.nseq1)then
	score(i,nseq2)=-20
	end if
 	val(i,0)=(go+ge*(i-1))
	end do
	do j=1,nseq2
	if(j.lt.nseq2)then
	score(nseq1,j)=-20
	end if
 	val(0,j)=(go+ge*(j-1))
	end do


C Perform DP computation
	
        DO j=1,NSEQ2	
        DO i=1,NSEQ1
             H=VAL(i-1,j)
	if(j .lt.nseq2)Then
             if (DIR(i-1,j)) THEN
               H=H+GAP_EXTN
             ELSE
               H=H+GAP_OPEN
             ENDIF
	else
		if (DIR(i-1,j)) THEN
               H=H+ge
             ELSE
               H=H+go
             ENDIF
	END IF
             V=VAL(i,j-1)
	IF(I.lt.Nseq1)then
             if (DIR(i,j-1)) THEN
               V=V+GAP_EXTN
             ELSE
               V=V+GAP_OPEN
             ENDIF
	else
		if (DIR(i,j-1)) THEN
               V=V+ge
             ELSE
               V=V+go
             ENDIF
	end if
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
           ENDDO
	END DO
c	update the score matrix for residues above
c	residue j is updated
	align=VAL(NSEQ1,NSEQ2)
         RETURN
       END

C
C  Extract alignment from the distance matrix
C
      FUNCTION IEXTRACT(NSEQ1,NSEQ2)
	PARAMETER(NMAX=3000)
      PARAMETER(MAXRES=1800)
	       
	INTEGER POS
	LOGICAL*1 DIR

	COMMON/ALN2/DIR(0:MAXRES,0:NMAX)		
       COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX), 
     *   LIGNMENT(0:2*MAXRES), GAP_OPEN, GAP_EXTN


C Traceback alignment

         i=NSEQ1
         j=NSEQ2
         pos=0
         
         DO WHILE ((i.GT.0).AND.(j.GT.0))
           IF (.not.DIR(i,j)) THEN
             i=i-1
             j=j-1
             LIGNMENT(pos)=0
           ELSE
             IF (VAL(i-1,j).GT.VAL(i,j-1)) THEN
               i=i-1
               LIGNMENT(pos)=1
             ELSE
               j=j-1
               LIGNMENT(pos)=-1
             ENDIF
           ENDIF
           pos=pos+1
         ENDDO
    
         DO WHILE (i.GT.0)
           i=i-1
           LIGNMENT(pos)=1
           pos=pos+1
         ENDDO

         DO WHILE (j.GT.0)
           j=j-1
           LIGNMENT(pos)=-1
           pos=pos+1
         ENDDO

         DO i=0,pos/2-1
           itmp=LIGNMENT(i)
           LIGNMENT(i)=LIGNMENT(pos-1-i)
           LIGNMENT(pos-1-i)=itmp
         ENDDO         

         iextract=pos

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
      subroutine TMscore(L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,TM,Rcomm,Lcomm)
      PARAMETER(nmax=3000)
      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align_2/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score

      dimension xa_max(nmax),ya_max(nmax),za_max(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max

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
      common/scores/score
      double precision score,score_max

      n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         score_sum=score_sum+1/(1+(dis/d0)**2)
      enddo
      score=score_sum/float(nseqB) !TM-score

      return
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


