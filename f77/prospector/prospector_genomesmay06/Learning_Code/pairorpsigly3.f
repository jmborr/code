c       This program generates for each query sequence one file
c       containing the orientation dependent contact potential
c
c	like pmfhomgsororienbaln3590.f but uses 50a1 input of msac like
c	pmfhomgsor.f but includes orientational information
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
	PARAMETER(NMAX=1700)    !maximum sequence length for template
	PARAMETER(NCASES=7100)  !maximum number of read templates
        PARAMETER(MAXSEQ=300)   !maximum number of aligned sequences in the msa
	PARAMETER(NRESMAX=530)  !maximum sequence length for query
	PARAMETER(IWIN=6)       !a sequence fragment size
	PARAMETER(IWIN2=13)     !1+2*iwin
	CHARACTER*3 AA3(0:21),AD3
	CHARACTER*1 AD
	CHARACTER*6 TARGET !prefix for one run, of the form "pdbbX" or "pdbbXX" or "pdbXXX"
	CHARACTER*4 GENOME
	CHARACTER*5 NAME2,NAMET
	CHARACTER*5 NAMESELECT
	CHARACTER*255 F3590 !full path for file containing 3590 alignments for queries
	integer n3590 !number of non-blank characters in "F3590" variable
        CHARACTER*255 F3590ALN !full path for file containing file extension to 3590 alignment
	integer naln90 !number of non-blank characters in "F3590ALN" variable
	CHARACTER*255 NAMES(200) !query identifiers
	LOGICAL *1 TOUCH(-1:NRESMAX,-1:NRESMAX)
	COMMON/NAMES/NAME(NCASES) !name(i) contain template identifier for "i"th template
	CHARACTER*5 NAME
	COMMON/C/NDATA !number of templates
	integer nc !nc(i,j) is number of contacts that residue "i" makes in template "j"
	integer icon !icon(i,j,k) is sequence-position of the residue making the "i"th contact
	             !with residue "j" in template "k"
	COMMON/CONT/NC(NMAX,NCASES),ICON(15,NMAX,NCASES)
	COMMON/CONT2/LOR(15,NMAX,NCASES) !lor(i,j,k) is the orientation of the "i"th contact
				!made by residue icon(i,j,k) with residue "j" in template "k"	
	integer nresseq !query seq. length for purposes of dealing with multiple seq. alignments
	integer nseq_ori !read number of aligned sequences (including query sequence)
	integer kcode    !kcode(i,j) amino acid or gap numeric code query-sequence-position "i"
			 !in aligned sequence "j". Numeric indexes are those of data AA/
        COMMON/SEQ/NRESSEQ,NSEQ_ORI,KCODE(NRESMAX,MAXSEQ)
	COMMON/SEQ2/IMUT(0:19,0:19) !BLOSUM scoring matrix
	COMMON/SEQ3/NSUC(NRESMAX) !never used
	integer nex !there are nex(i) residues of type "i" in query sequence
	integer ipos !ipos(i,j) is the query-sequence position for the "i"th amino acid of type
	             !"j". There are a total of nex(j) amino acids of type "j" in query sequence
	COMMON/SEQ4/NEX(0:19),IPOS(NRESMAX,0:19)
	real aobs !aobs(i,j,k) weighted number of times we see a pseudo-contact
	          !between "i" and "j" with orientation "k"
	real aexp !never used
	COMMON/SCALE/AEXP(NRESMAX,NRESMAX),AOBS(-4:NRESMAX,NRESMAX,3)
	LOGICAL HOMMAP
	COMMON/SEQHOM/HOMMAP(NCASES) !hommap(i).eq. .true. means template "i" is non-homolog
	                             !to query-sequence
	COMMON/SCALEPAIR/APABL(-1:19,-1:19,3) !orientation dependent potential. Three orientations
	                                      !(parallel,perpendicular,anti-parallel)
	COMMON/SCALEPAIR2/AN(NRESMAX,NRESMAX,3)!orientation dependent potential as apabl, but
	                                       !parametrized for the query-sequence
	DIMENSION JCODE(NMAX)!jcode(i) numeric code for residue at position "i" of query sequence
	DIMENSION AMAP(NRESMAX,NRESMAX)
	DIMENSION ICODE(NMAX,NCASES) !icode(i,ik) is amino acid numeric index for residue position
                                     !"i" of "ik"th template. Indexes are those of data AA/
        DIMENSION NRES(NCASES) !nres(i) is sequence length for "i"th template
c       imap(i,j,k).eq.1 means residue "i" and "j" make a contact with orientation "k"

c	directories
	common/namesize/nsize
	integer nsize !nsize is length of word identifying each template
	character*255 home !full path of working directory
	character*255 templatedat !full path for directory of templates
	common/home/nhome,home
	integer nhome !number of non-blank characters in "home" variable
	integer nadr  !number of non-blank characters in "templatedat" variable
c	=================================
c	======================================================
	DATA AA3/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &		     'MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP','CYX','HIX'/

c	======================================================

c	======================================================
c      			 PSIBLAST INPUTS
        logical*1 tncbi(0:22) !tncbi(i).eq..true if "i" in NCBI index of one of the 20 amino acids
c       xm(0:19,j,k) is position-specific scoring matrix at residue-position "j" for "k"th template
        common/profile/xm(0:19,0:nmax,0:ncases)
        dimension ix(0:22)
        dimension ncbi(22) !relate NCBI and Jeff's amino acid codes numbering
	character*1 aa(0:19),ancbi(1:22)
c       Jeff's amino acid codes
        DATA AA/ 'G','A','S','C','V','T','I','P',
     &               'M','D','N','L',
     &               'K','E','Q','R',
     &               'H','F','Y','W'/
c       ancbi amino acid codes from the NCBI database
        data ancbi/'A','B','C','D','E','F','G','H','I','K','L','M','N',
     &  'P','Q','R','S','T','V','W','X','Y'/


c       initialize "tncbi" and "ncbi"
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

c       initialize "home" and "nhome"
	open(unit=1,file='homedir') !file containing full path of working directory
	read(unit=1,fmt='A255')home
	nhome=0
	do i=255,1,-1
	   if(home(i:i) .eq. ' ')nhome=i
	end do
	nhome=nhome-1
c	write(*,*) 'home=',home(1:nhome),'nhome=',nhome

c       initialize "templatedat" and "nadr"
	open(unit=1,file='templatedir') !file containing full path for directory of templates
	read(unit=1,fmt='A255')templatedat
	nadr=0
	do i=255,1,-1
	   if(templatedat(i:i) .eq. ' ')nadr=i
	end do
	nadr=nadr-1
c	write(*,*) 'templatedat=',templatedat(1:nadr),'nadr=',nadr

c	init "f3590" and "n3590"
c	'3590' is file containing full path for file containing 3590 alignments for queries
	open(unit=1,file='3590') 
	read(unit=1,fmt='A255')f3590
	n3590=0
	do i=255,1,-1
	   if(f3590(i:i) .eq. ' ')n3590=i
	end do
	n3590=n3590-1
c	write(*,*) '3590=',f3590(1:n3590),'n3590=',n3590

c	init "3590aln" and "naln90"
	open(unit=1,file='3590aln')!file containing file extension for 3590 alignments
	read(unit=1,fmt='A255')f3590aln
	naln90=0
	do i=255,1,-1
	   if(f3590aln(i:i) .eq. ' ')naln90=i
	end do
	naln90=naln90-1
c	write(*,*) '3590aln=',f3590aln(1:naln90),'naln90=',naln90	
	
c       we invoque this file like "pairorpsigly3 < nametargXXX"
	read(5,1112) target !prefix for this run, of the form pdbbX or pdbbXX or pdbXXX
1112	format(a6)
c       file namesize contains length of word identifying each template
	Open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)

c       this file will store identifiers of processed queries
	open(unit=3,file='completedpairorpsigly'//target)

c       init apabl(-1:19,-1:19)
	do i=-1,19	
	   do j=-1,19
	      do it=1,3
		 apabl(i,j,it)=0.
	      end do
	   end do
	end do

c       init apabl(0:19,0:19)
c       !file PAIRORGLY3 contains info for orientation dependent potential
	OPEN(unit=26,file=home(1:nhome)//'strlist/PAIRORGLY3')
	do it=1,3
	   read(26,*)
	   read(26,*)
	   do i=0,19
	      read(26,17) ad3, (apabl(i,j,it),j=0,19)
	   enddo
	end do
17	FORMAT(A3,20F5.1)
	close(26)

c       init "imut"
	OPEN(unit=20,file=home(1:nhome)//'strlist//BLOSUM_INT')!file contains BLOSUM scoring matrix
	read(20,*)
	do i=0,19
	   read(20,*)ad,(imut(i,j),j=0,19)
	end do

c       init "name", "nres", "icode", "xm"
	OPEN(unit=50,file='LIST.may2006good35')!file contains template identifiers
	read(50,*)NDATA
	NTOT=0 !number of residues of all templates
	DO IK=1,NDATA
	   read(50,221)name(ik) !template identifier
221	   format(a5)
c       file /library/orien6/101m_.SEQ contains amino acid sequence for
c       101m_ template, but not in the form of one-letter or
c       three-letter code, but in numeric form. For example, "4" means
c       the amino acid at index 4 in data block AA, which corresponds to
c       "V" (remember AA begins at index 0, not 1)
	   OPEN(unit=30,file='/library/orien6/'//NAME(ik)//'.SEQ')
	   rewind(30)
	   READ(30,2)NAME2,Nres(ik) !template identifier and sequence length
	   if(nres(ik).gt.nmax)nres(ik)=nmax !restrict sequence length to "nmax"
	   read(30,*)(icode(i,ik),i=1,nres(ik)) !store sequence in numeric code
	   do i=1,nres(ik)
	      if(icode(i,ik).gt.20)then !strange codes are mutated to GLY
		 icode(i,ik)=0
		 stop
	      end if
	   end do
	   close(30)
2	   format(a5,1x,i5)
	   do i=1,nres(ik)
	      do ires=0,19
		 xm(ires,i,ik)=0. !init all profiles
	      end do
	   end do
c       profile/101m_.mtx contains the PSI-BLAST-derived position specific scoring matrix
	   open(unit=31,file=templatedat(1:nadr)//'profile/'//name(ik)//
     &          '.mtx')
	   read(31,*)nres(ik) !first line contains template sequence length
	   if(nres(ik).gt.nmax)nres(ik)=nmax !restrict to first nmax amino acids
	   do id=1,13 !obviate next 13 lines
	      read(31,*)
	   end do
	   do i=1,nres(ik)
c       there are 22 "amino acids" corresponding to the NCBI amino acid codes of data ancbi/
	      read(31,*)(ix(ires),ires=0,22)
	      do jres=1,22
		 if(tncbi(jres))then !if the code corresponds to one of the 20 amino acids
		    ires=ncbi(jres) !translate from NCBI code to Jeff's code
		    xm(ires,i,ik)=float(ix(jres))/100. !PSI-BLAST multiplies scores by 100
		 end if
	      end do
	   end do
	END DO
	do ik=1,ndata
	   NTOT=NTOT+NRES(ik)
	end do
	close(50)!close the file containing template identifiers (LIST.may2006good35)

	call contb(NRES) !init "nc" "icon" "lor"

c	init "names"
	OPEN(unit=60,file='LIST.targ'//target)
	read(60,*)ndataseq,genome !number of query sequences and prefix
	OPEN(unit=25,file='problems')
	do 10000  iik=1,ndataseq !cycle over all query sequences
	   READ(60,33)NAMES(IIK) !read query identifiers
33	   format(a)
	   OPEN(UNIT=6,FILE='outhomgsororien.'//NAMES(IIK)(1:nsize))
	   nunit=31
c       file 101m_.aln (in this case 101m_ would be the query) contains
c       the multiple sequence alignment compiled from PSI-BLAST results
	   OPEN(UNIT=31,FILE=f3590(1:n3590)//names(iik)(1:nsize)//f3590aln(1:naln90))
	   call read_msa(nunit) !init "nresseq" "nseq_ori" "kcode"
	   itest=0
	   do i=1,nresseq !cycle over query seq. length and init "jcode"
	      jcode(i)=kcode(i,1) !kcode(i,1) is numeric code for amino acid at position "i" of
	   end do                 !aligned sequence 1, which is the query sequence
	   close(31)

	   ntot=0
	   nz=0
	   do ik=1,ndata
	      NTOT=NTOT+NRES(ik) !total number of residues among all templates
	      nz=nz+1 !number of templates
	   end do
	   do jk=1,ndata !init all templates as non-homologs to the query sequence
	      hommap(jk)=.true.
	   end do

c       init "nex" and "ipos"
	   mres=nresseq
	   do iseq=0,19
	      nex(iseq)=0
	   end do
	   do i=1,mres !cycle over query-sequence length
	      do iseq=0,19
		 if(jcode(i) .eq. iseq)then
		    nex(iseq)=nex(iseq)+1  !update the query-sequence amino acid composition
		    ipos(nex(iseq),iseq)=i !mark the position of this amino acid type
		 end if
	      end do
	   end do

	   CALL EHMSTAT(NAME,NRES,ICODE)

	   mres=nresseq !query-sequence length
c       homedir/genomePAIROR10ml_.PAIRORPSIGLY3 contains orientation dependent potential for 101m_
	   OPEN(UNIT=7,FILE=home(1:nhome)//genome//'PAIROR/'//NAMES(iik)(1:nsize)//'.PAIRORPSIGLY3')
	   rewind(7)
	   icorrect=0
	   ineg=0
	   do it=1,3 !cycle over all orientations
	      do i=1,mres !cycle over all query sequence
		 iseq=jcode(i) !amino acid type at position "i"
		 write(7,*)i,aa3(iseq) !aa3(iseq) is amino acid type in the three letter code
		 WRITE(7,171) (an(i,j,it),j=1,mres) !orientation specific potential
		 write(7,*)'================================='
171		 FORMAT(25F5.1)
	      end do
	   end do
	   close(7)
c	----------------------
	   write(3,3)names(iik)(1:nsize)
3	   format(a)
10000 	continue
	STOP !end of program pairorpsigly3
	END


c       cccccccccccccccccccccccccccccccccccccccccccccccccccc
c       arguments to the subroutine:
c       name(i) contain template identifier for "i"th template
c       nres(i) is sequence length for "i"th template
c       icode(i,ik) is amino acid numeric index for residue position "i" of
c       "ik"th template. Indexes are those of data AA/
c       cccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE EHMSTAT(NAME,NRES,ICODE)
	PARAMETER(NMAX=1700)   !maximum sequence length for template
        PARAMETER(MAXSEQ=300)  !maximum number of aligned sequences in the msa
	PARAMETER(NCASES=7100) !maximum number of read templates
	PARAMETER(NRESMAX=530) !maximum sequence length for query
	PARAMETER(IWIN=6)
	PARAMETER(IWIN2=13)
	LOGICAL HOMMAP
	CHARACTER*5 NAME(NCASES) !name(i) contain template identifier for "i"th template
	DIMENSION ICODE(NMAX,NCASES) !icode(i,ik) is amino acid numeric index for residue position
	                             !"i" of "ik"th template. Indexes are those of data AA/
	DIMENSION NRES(NCASES) !nres(i) is sequence length for "i"th template
c       imap(i,j,k).eq.1 means residue "i" and "j" make a contact with orientation "k"
	DIMENSION IMAP(-5:NMAX,-5:NMAX,3)
c       ihist(i) is the number of instances we calculate a smoothed
c       score at position "i". This will tipically be the number of
c       identical residues to "i" in the template library.
	DIMENSION IHIST(NRESMAX)	
	DIMENSION JCODE(NRESMAX)
        DIMENSION AMUT(-1:19,0:NRESMAX)!amut(-1:19,j) amino acid profile at query-sequence-position
	                               !"j". -1 is reserved for non-amino acid (gap or other)
c       anat(i) is a measure of how disperse is the profile in the [i-iwin,i+iwin] window of the
c       query sequence, assuming that this window does not contain gaps or some atypical amino
c       acid type in the query sequence, that is, nogap(i,1) .eq. .true.
c       anat(i)=sum_{iv=i-iwin}^{i+iwin} sum_{ires=0}^19 amut(ires,iv)*amut(ires,iv)
	DIMENSION ANAT(NRESMAX)
c       ahist(i) is the cummulative of the smoothed scores at position
c       "i". This cummulative step is summing over all template
c       positions "ij" satisfying icode(ij,jk).eq.jcode(i) and over all
c       templates "jk". Note that if the [i-iwin,i+iwin] contains a gap
c       or a non-conventional amino acid, then ahist(i).eq.0
        DIMENSION AHIST(-IWIN2:NRESMAX)
c       ahist2(i) is the cummulative of the square of the smoothed scores
        DIMENSION AHIST2(-IWIN2:NRESMAX)
c       aleft(i) is a smoothed score which is three standard deviations above the average of
c       the distribution of smooted scores at query-sequence position "i"
	DIMENSION ALEFT(-IWIN2:NRESMAX)
	COMMON/C/NDATA
	COMMON/CONT/NC(NMAX,NCASES),ICON(15,NMAX,NCASES)
	COMMON/CONT2/LOR(15,NMAX,NCASES)	
	COMMON/SEQHOM/HOMMAP(NCASES) !hommap(i).eq. .true. means template "i" is non-homolog
	                             !to query-sequence
	integer nresseq !query seq. length for purposes of dealing with multiple seq. alignments
	integer kcode    !kcode(i,j) amino acid or gap numeric code query-sequence-position "i"
			 !in aligned sequence "j". Numeric indexes are those of data AA/
        COMMON/SEQ/NRESSEQ,NSEQ_ORI,KCODE(NRESMAX,MAXSEQ)
	COMMON/SCALEPAIR/APABL(-1:19,-1:19,3)!orientation dependent potential. Three orientations
                                             !(parallel,perpendicular,anti-parallel)
	COMMON/SCALEPAIR2/AN(NRESMAX,NRESMAX,3)!orientation dependent potential as apabl, but
                                               !parametrized for the query-sequence
	COMMON/SEQ2/IMUT(0:19,0:19)
	COMMON/SEQ3/NSUC(NRESMAX) !never used
	COMMON/SEQ4/NEX(0:19),IPOS(NRESMAX,0:19)
	real aobs !aobs(i,j,k) weighted number of times we see a pseudo-contact
	          !between "i" and "j" with orientation "k"
	real aexp !never used
	COMMON/SCALE/AEXP(NRESMAX,NRESMAX),AOBS(-4:NRESMAX,NRESMAX,3)
	logical NOGAP
c       nogap(i,1).eq.true none of the amino acids in the window
c       [i-iwin,i+iwin] of the query sequence is a gap or a
c       non-conventional amino acid
	COMMON/SEQ5/NOGAP(NRESMAX,1)
        COMMON/PROFILE/XM(0:19,0:NMAX,0:NCASES)


c       init "ahist" "ahist2" "ihist" "jcode" "an" "aobs"
	mres=nresseq
	do i=1,mres !cycle over query-sequence length
	   ahist(i)=0.
	   ahist2(i)=0.
	   ihist(i)=0
	   jcode(i)=kcode(i,1)!numeric code for residue at position "i" of query sequence
	   do j=1,mres
	      do it=1,3
		 an(j,i,it)=0.
		 aobs(j,i,it)=0.
	      end do
	   end do
	end do

c       init "amut"
        do i=1,mres !cycle over query-sequence length
	   do iseq=0,19
	      amut(iseq,i)=0.
	   end do
	   nden=0
	   do ic=1,nseq_ori !cycle over aligned sequences in the 3590 msa
	      ires=kcode(i,ic) !amino acid type (numeric code)
	      nden=nden+1
	      if(ires.gt.-1)then !ires.eq.-1 means a gap, or just non-amino acid type
		 amut(ires,i)=amut(ires,i)+1 !update amino acid profile
	      end if
	   end do  
	   do iseq=0,19 !normalize profile
	      amut(iseq,i)=amut(iseq,i)/nden
	   end do
        end do

c       init "nogap" "nsuc" "anat"
	IHOM=1 !denotes the query sequence for the set of aligned sequences in the 3590 msa
	do 150 i=1,mres !cycle over query-sequence length
	   anat(i)=0.
	   nsuc(i)=0
	   NOGAP(I,IHOM)=.TRUE.
	   ii=0
	   do ii=-iwin,iwin
	      iv=i+ii
	      if(iv .gt. 0 .and. iv.le.mres)then
		 jres=jcode(iv) !amino acid type at query-sequence position "iv"
		 if(jcode(iv).eq. -1)then
		    NOGAP(I,IHOM)=.FALSE. !the [i-iwin,i+iwin] fragment of the query sequence
		    go to 150             !contains a gap. Stop and go to next position "i"
		 else
		    do ires=0,19
		       anat(i)=anat(i)+amut(ires,iv)*amut(ires,iv)
		    end do
		 end if
	      end if
	   end do
	   if(anat(i).lt.0)NOGAP(I,IHOM)=.false.
150	continue


c       init "ahist"  "ahist2" "ihist"
	do 700 jk=1,NDATA !cycle over all templates
	   if( .not.hommap(jk))go to 700 !skip if template is homolog to query sequence
	   DO 161 ij=1,nres(jk) !cycle over the template sequence
	      do i=1,mres !cycle over the query sequence
		 if(nogap(i,ihom))then !if [i-iwin,i+iwin] of query seq. does not contain gaps or
		                       !a non-conventional amino acid
		    aSIM=0. !average score
		    is1=jcode(i) !amino acid type at position "i" of query sequence
		    js1=icode(ij,jk) !amino acid type at position "ij" of template "jk"
		    if(is1 .eq.js1)then !if types coincide
c       find smoothed score between query-sequence position "i" and
c       template-position "ij" as the sum of the scores resulting from
c       the gapless alignment of window [i-iwin,i+iwin] in query
c       sequence to window [ij-iwin,ij+iwin] in template sequence. You
c       can also think of this score as the cummulative score over the
c       window.
		       do ii=-iwin,iwin !open corresponding windows in query and template seqs.
			  iv=i+ii  !running index over query
			  jv=ij+ii !running index over template
			  if(iv .gt. 0 .and. iv.le.mres)then
			     if(jv .gt. 0 .and. jv.le.nres(jk))then
c       The score of query-sequence position "iv" when aligned to
c       template-sequence position "jv" is given by the "matrix product"
c       of the query profile at "iv" by the position-specific scoring
c       matrix of template "jk" at position "jv"
				do jres=0,19 
				   asim=asim+amut(jres,iv)*xm(jres,jv,jk)
				end do
			     end if
			  end if
		       end do
c       ahist(i) is the cummulative of the previous smoothed scores. The
c       second cummulative step is summing over all template positions
c       "ij" satisfying icode(ij,jk).eq.jcode(i) and over all templates "jk"
		       ahist(i)=ahist(i)+asim
c       ahist2(i) is the cummulative of the square of the previous smoothed scores
		       ahist2(i)=ahist2(i)+asim*asim
c       ihist(i) is the number of instances we calculate a smoothed
c       score at position "i". This will tipically be the number of
c       identical residues to "i" in the template library.
		       ihist(i)=ihist(i)+1
		    end if
		 end if
	      end do
161	   CONTINUE
700	CONTINUE


c       init "aleft"
	do i=1,mres !cycle over the query sequence
	   if( nogap(i,ihom) )then ![i-iwin,i+iwin] only contains conventional amino acids
	      if( ihist(i).gt.0 )then !we found residues identical to "i" in the template library
		 nd=ihist(i) !number of such residues
		 ah=ahist(i)/nd !average of the smoothed scores
		 ah2=ahist2(i)/nd !normalize cummulative of square of smoothed scores
		 std=sqrt(ah2-ah*ah) !standard deviation of the distribution of smoothed scores
c       aleft(i) is a smoothed score which is three standard deviations above the average of
c       the distribution of smooted scores at query-sequence position "i"
		 aleft(i)=ah+3*std 
	      end if
	   end if
	end do

	do 1000 jk=1,NDATA !cycle over all templates
	   if(.not.hommap(jk))go to 1000 !skip homologs templates to query sequence
	   isimt=0
c       init "imap" for current template
	   do i=1,nres(jk)+1 !cycle over template residues
	      do j=i,nres(jk)+1 !cycle over template residues
		 do it=1,3 !cycle over all orientations
		    imap(i,j,it)=0
		    imap(j,i,it)=0
		 end do
	      end do
	   end do
	   do ij=1,nres(jk) !cycle over template residues
	      if( nc(ij,jk).gt.0 )then !if "ij" makes contacts with other residues in the template
		 do kk=1,nc(ij,jk) !cycle over contacts betwee "ij" and other residues in the templ
		    jj=icon(kk,ij,jk) !residue making contact number kk with "ij"
		    it=lor(kk,ij,jk)  !orientation of the contact between "ij" and "jj"
		    imap(ij,jj,it)=1  !mark there's a contact between "ij" and "jj" and its orient
		 end do
	      end if
	   end do
c	sequence is ik
c	structure is jk
	   do 171 ij=1,nres(jk) !cycle over template residues
	      iseq1=icode(ij,jk) !amino acid type at position "ij" in template "jk"
	      if( nc(ij,jk).gt.0 )then !if "ij" makes contacts with other residues in the template
		 do 170 kk=1,nc(ij,jk)
		    jj=icon(kk,ij,jk)!residue making contact number kk with "ij"
		    it=lor(kk,ij,jk) !orientation of the contact between "ij" and "jj"
		    imap(ij,jj,it)=1 !mark there's a contact between "ij" and "jj" and its orient
		    if(iabs(jj-ij) .gt.iwin)then !not a local contact
		       iseq2=icode(jj,jk) !amino acid type at position "jj" in template "jk"
c       check that query sequence contains amino acids of type iseq1 and iseq2
		       if(nex(iseq1) .gt. 0. and. nex(iseq2).gt.0)then
			  do 169  i1=1,nex(iseq1)!cycle over query-seq amino acids of type iseq1
			     i=ipos(i1,iseq1) !position in the query sequence
			     do 168 j1=1,nex(iseq2) !cycle over query-seq amino acids of type iseq2
				j=ipos(j1,iseq2) !position in the query sequence
c       check that positions "i" and "j" are centers of respective
c       windows containing only conventional residues
				if(nogap(i,ihom) .and. nogap(j,ihom))then 
				   if(iabs(i-j).gt. iwin)then !contact order not low
c       calculate smoothed score between position "i" in query sequence
c       and position "ij" in template "jk". Remember "i" and "ij" are amino acids of same type
				      asim1=0.
				      do ii=-iwin,iwin !cycle over respective windows
					 iv=i+ii
					 jv=ij+ii
					 if(iv .gt.0 .and. iv.le.mres)then
					    if(jv .gt. 0 .and. jv.le.nres(jk))then
					       do jres=0,19 
						  asim1=asim1+amut(jres,iv)*xm(jres,jv,jk)
					       end do					       
					    end if
					 end if
				      end do
c       calculate smoothed score between position "j" in query sequence
c       and position "jj" in template "jk". Remember "j" and "jj" are amino acids of same type
				      asim2=0.
				      do ii=-iwin,iwin
					 iv=j+ii
					 jv=jj+ii
					 if(iv .gt.0 .and. iv.le.mres)then
					    if(jv .gt.0 .and. jv.le.nres(jk))then
					       do jres=0,19
						  asim2=asim2+amut(jres,iv)*xm(jres,jv,jk)
					       end do
					    end if
					 end if
				      end do
				      is1=jcode(i) !type of amino acid at position "i" (iseq1)
				      js1=jcode(j) !type of amino acid at position "j" (iseq2)
c       we require both smoothed scores are significantly bigger than
c       the average smoothed scores ahist(i) and ahist(j), that is,
c       bigger than aleft(i) and aleft(j) which are three standard
c       deviations away from the averages. This says that the stretch of
c       amino acids in "i" for query sequence is "homolog" to the
c       stretch of amino acids in "ij" in the template. The same goes
c       for "j" and "jj". Thus we are justified in tranferring in some
c       way the contact between "ij" and "jj" (if there is one) as
c       contact between "i" and "j". However, instead of this simplified
c       scheme, we acknowledge that there may be some noise in the
c       alignments and that we could transfer contacts between neighbor
c       amino acids of "ij" and "jj" as contacts between neighboring
c       amino acids of "i" and "j". Thus, we open an 11-residue window
c       around "ij" and "jj". Then, we gather all contacts between
c       residues in these two windows. If say, there's a contact between
c       positions ijz and jjz in the template (ijz is neighbor of ij and
c       jjz is neighbor of jj), then we look at corresponding positions
c       in the query sequence iz and jz, and if these positions are
c       evolutary related to the positions in the template, then we can
c       transfer the (ijz,jjz) contact as (iz,jz) contact
				      if(asim1 .ge. aleft(i) .or. asim2 .ge.aleft(j))then
					 do iv=-5,5 !open a window of size 11. It could have been 2*iwin+1 instead
					    iz=i+iv !running index in query sequence centered in "i"
					    ijz=ij+iv !running index in template centered in "ij"
					    if(iz .gt. 0 .and. ijz .gt.0)then
					       if(iz .le. mres .and. ijz .le.nres(jk))then
						  iz1=jcode(iz) !amino acid type in query
						  ijz1=icode(ijz,jk) !amino acid type in template
c       if favorable blosum matrix score (naturally occurring mutation) between query residue and template residue
						  if(imut(iz1,ijz1) .gt. 0 )then 
						     do jv=-5,5 !open a window of size 11
							jz=j+jv !running index in query sequence centered in "i"
							jjz=jj+jv !running index in template centered in "jj"
							if(jz .gt.0 .and. jz.le.mres)then
							   if(jjz .gt.0 .and. jjz.le.nres(jk))then
							      jz1=jcode(jz) !amino acid type in query
							      jjz1=icode(jjz,jk) !amino acid type in template
							      js1=kcode(jz,1)
c       we require not only smoothed scores much better that the average, but also positive. This condition should be moved right before the line "if(asim1 .ge. aleft(i) .or. asim2 .ge.aleft(j))then"
							      if(asim1 .gt. 0 .and. asim2 .gt. 0)then
c       if favorable blosum matrix score (naturally occurring mutation) between query residue and template residue
								 if(imut(jz1,jjz1) .gt. 0)then
c       weight of the contact is the product of the two smoothed scores at positions "i" and "j" divided by the product of the dispersion of each profile. 
								    sf=asim1*asim2/(anat(i)*anat(j))
								    sf=0.2+sqrt(sf)
								    do it=1,3
c       accumulate number of observations where we see a pseudo-contact between iz and jz in aobs
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
					 end do
					 nsuc(i)=nsuc(i)+1 !never used
					 nsuc(j)=nsuc(j)+1 !never used
				      end if
				   end if
				end if
168			     continue
169			  continue
		       end if
		    end if
170		 continue
	      end if
171	   continue
1000	continue

1001	continue

c       calculate anorm, normalization factor
	atot=0.
	do i=1,mres !cycle over query sequence
	   do j=1,mres !cycle over query sequence
	      do it=1,3 !cycle over all orientations
		 atot=atot+aobs(i,j,it)
	      end do
	   end do
	end do
	anorm=atot/(3*mres*mres)
	write(6,*)'number of observations',anorm,atot
	isum=0
	do i=1,mres !cycle over query sequence
	   do j=1,mres !cycle over query sequence
	      do it=1,3	 !cycle over all orientations
		 if(aobs(i,j,it) .gt. 0.)then
		    an(i,j,it)=an(i,j,it)-alog(aobs(i,j,it)/anorm) !specific orientation potential
		    an(j,i,it)=an(j,i,it)-alog(aobs(i,j,it)/anorm)
		    isum=isum+1
		 else !use generic orientation potential if we don't find any pseudo-contact
		    iseq=jcode(i) !amino acid type at query-sequence position "i"
		    jseq=jcode(j)
		    an(i,j,it)=an(i,j,it)+apabl(iseq,jseq,it) !generic orientation potential
		    an(j,i,it)=an(j,i,it)+apabl(iseq,jseq,it)
		 end if
	      end do
	   end do
	end do
2000	continue
	do i=1,mres
	   do j=1,mres
	      do it=1,3
		 an(i,j,it)=an(i,j,it)/(2)
	      end do
	   end do
	end do
	RETURN !end of ehmstat subroutine
	END 	
	

c	*********************************************************
c	*        constructs the contact map library		*
c	* obtained for NEWDATA base and is the mean position	*
c	*	ALL coordinates are in ANGSTROMS		*
c	*********************************************************
	subroutine contb(NRES)
	PARAMETER(nmax=1700) !maximum sequence length for template
	PARAMETER(NCASES=7100) !maximum number of read templates
	CHARACTER*5 NAME
	common/names/name(ncases) !name(i) contain template identifier for "i"th template
	common/c/ndata !number of templates
	integer nc !nc(i,j) is number of contacts that residue "i" makes in template "j"
	integer icon !icon(i,j,k) is sequence-position of the residue making the "i"th contact
	             !with residue "j" in template "k"
	common/cont/nc(nmax,ncases),icon(15,nmax,ncases)
	common/cont2/lor(15,nmax,ncases) !lor(i,j,k) is the orientation of the "i"th contact
				!made by residue icon(i,j,k) with residue "j" in template "k"
	DIMENSION nres(ncases) !nres(i) is sequence length for "i"th template

c       init "nc" "icon" "lor"
	DO IK=1,NDATA
c       file /library/orien6fitgly/101m_.FITCONT_ORGLY contains contact map of 101m_. For example,
c       3    8    1    7    8   11  131  134  135  138
c       3    8    3    2    2    2    2    1    2    1
c       means residue 3 makes 8 contacts. The first line shows
c       sequence-position of contacting residues, and second line
c       contains the orientation (parallel,perpendicular,anti-parallel)
c       of each contact
	OPEN(unit=71,file='/library/orien6fitgly/'//NAME(ik)//'.FITCONT_ORGLY')
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
c       read multiple sequence alignment compiled from PSI-BLAST results
c       nunit: file descriptor for query.aln file
c       Initializes nresseq, nseq_ori, and kcode
c========================================================================
	subroutine read_msa(nunit)
	PARAMETER     (MAXSEQ=300) !maximum number of aligned sequences in the msa
	PARAMETER     (MAXRES=530) !maximum sequence length for query
	PARAMETER     (NRESMAX=530)!maximum sequence length for query
	PARAMETER     (KMAX=10000)
	CHARACTER     SEQ*1
	CHARACTER*1   AA1(0:19),SEQD(KMAX)
	DIMENSION SEQ(MAXSEQ,MAXRES) !seq(i,j) amino acid code or gap symbol for
                                     !query-sequence-position "j" in aligned sequence "i"
	integer nresseq  !query seq. length for purposes of dealing with multiple seq. alignments
	integer nseq_ori !read number of aligned sequences (including query sequence)
	integer kcode    !kcode(i,j) amino acid or gap numeric code query-sequence-position "i"
			 !in aligned sequence "j". Numeric indexes are those of data AA/
	common/seq/nresseq,nseq_ori,kcode(nresmax,maxseq)
c       Jeff's amino acid codes
	DATA AA1/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/

	read(nunit,*)nseq_ori,nresseq !read number of aligned sequences and query-sequence length
	if(nresseq.gt.maxres)then 
	   nresseq=maxres !restrict alignment length to maxres
	end if
	
	if(nseq_ori.gt.maxseq)nseq_ori=maxseq !restrict number of aligned sequences to maxseq
	if(nresseq.le.maxres)then
	   do i=1,nseq_ori !cycle through the aligned sequences
	      read(nunit,1)(seq(i,j),j=1,nresseq) !read aligned sequence "i"
 1	      format(50a1) !every row assumed to contain 50 characters
	   end do
	else !query sequence length bigger than maxres
	   do i=1,nseq_ori !cycle through the aligned sequences
	      read(nunit,1)(seqd(j),j=1,nresseq) !read aligned sequence to temporary array
	      do j=1,maxres
		 seq(i,j)=seqd(j) !restrict length to maxres
	      end do
	   end do
	   nresseq=maxres !update nresseq to maxres
	end if

	
c	&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	do i=1,nseq_ori
	do j=1,nresseq
		kcode(j,i)=-1 !initialize as non-amino acid
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

