ccccc 
c
c     'pair.f' is a program modified from Jeff's program of extracting
c     pair1 and pair3 potential.
c
c     Beginning with a multiple sequence alignment of homologues to a
c     target sequence, and a list of pdb structures.
c
c     pair1(i,j) is a LchxLch position potential contact potential,
c     independent of secondary structure propensities and orientation of
c     side-chains. For a given (i,j), the potential is an average over
c     all homologues of matrix1.comm(S_i,S_j), where S_i,S_j are the
c     residue types at positions i,j for a given homologue.
c
c     pair3(i,j,it) is a LchxLchx3 amino acid contact potential, which
c     is constructed from two cases:
c
c     (1) for a given (i,j), there may be one (or more) PDB structures
c     containing two sequence fragments centered at positions (ij,jj)
c     such that:
c     (a) residue type at ij same as residue type at i
c     (b) residue type at jj same as resideu type at j
c     (c) secondary structure at ij same as predicted sec.str. at i
c     (d) secondary structure at jj same as predicted sec.str. at j
c     Now we open a sequence fragment (of size 11) centered at ij
c     (fragment A) and another fragment centered at jj (fragment
c     B). Let's imagine there's a contact between position ijz in the
c     fragment A and position jjz in the fragment B,
c     imap(ijz,jjz,it)==1. Is there a contact between corresponding
c     positions in the target (iz,jz)? We assign a contact probability
c     aobs(iz,jz,it) if:
c     (e) There exists one (or more) template such
c     that a fragment in the template centered at i (fragment A') is
c     homolog to fragment A, and a fragment in the template centered at j
c     (fragment B') is homolog to fragment B.     
c     (f) residue type at iz in the previous homolog is evolutionary
c     related to residue type at position ijz in the template.
c     (g) residue type at jz in the previous homolog is evolutionary
c     related to residue type at position jjz in the template.
c     If all of the above is true, then
c     aobs(iz,jz,it)=aobs(iz,jz,it)+sf*imap(ijz,jjz,it), where sf is a
c     measure of the combined homology of fragent A to fragment A', and
c     fragment B to fragment B'
c
c     (2) for a given (i,j), if there are no template PDB structures
c     obeying the previous point, then the potential is an average over
c     all homologues of the probability matrix3.comm(S_i,S_j,it), where
c     S_i,S_j are the residue types at positions i,j for a given
c     homologue.
c
ccccc
      
      program pair3
      PARAMETER(Lch_str_max=2000) !maximum length of PDB structure
      PARAMETER(N_STR_MAX=20000) !maximum number of PDB structure
      PARAMETER(Lch_max=2000)   !maximum length of query sequence
      parameter(N_hom_max=300)  !maximum number of sequences in MSA
      COMMON/SEQ/Lch,N_MSA,KCODE(LCH_MAX,N_HOM_MAX)
      COMMON/SEQ2/IMUT(0:19,0:19)
      COMMON/SEQ4/NEX(0:19),IPOS(LCH_MAX,0:19)
      CHARACTER*3 AA(0:21),ad3
      character*1 ad,a1(Lch_str_max)
      character*20 template,dum
      common/names/template(N_str_max)
      common/c/N_stru
      dimension jcode(Lch_str_max)
      common/scalepair/apabl(-1:19,-1:19,3)
      COMMON/Scalepair2/AN(Lch_max,Lch_max,3)
      common/independent/ap(-1:19,-1:19),am(Lch_max,Lch_max)
      common/second/isec1(Lch_max),isec2(Lch_str_max,N_str_max)
      common/stru/ICODE(LCH_STR_MAX,N_STR_MAX),Lch_str(N_str_max)
      character*100 fnam1,fnam3
      character*300 fnam2
      common/libname/fnam2,ilib
      
      CHARACTER*1 SEQ(N_HOM_MAX,Lch_max)
      CHARACTER*1 seq0(Lch_max)
      CHARACTER*1 AA1(0:19)

      DATA AA/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     $     'MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG',
     $     'HIS','PHE','TYR','TRP','CYX','HIX'/

      DATA AA1/ 'G','A','S','C','V','T','I','P','M','D','N',
     $     'L','K','E','Q','R','H','F','Y','W'/

c       icode(i,str) --- AA of template structure
c       kcode(i,tem) --- AA of template in MSA
c       jcode(i)= kcode(i,1) -- AA of target sequence

ccccc easy/medium: asim_cut=0.3; hard: asim_cut=0.25
      call getarg(1,fnam1)      !fnam1: sequence identity cutoff
      call getarg(2,fnam2) !fnam2: library path
      ilib=0
      do i=1,300
c     char(32) is a blank space, that is, the end of the library path
         if(fnam2(i:i).eq.char(32))goto 52
         ilib=i
      enddo
 52   continue
      if(fnam1.eq.' ')then
         asim_cut=0.3
      else
         read(fnam1,*)asim_cut
      endif

ccccc input files
c     sec.struct. independent matrix
      OPEN(unit=1,file='matrix1.comm',status='old')
c     sec.struct. specific matrix
      OPEN(unit=2,file='matrix3.comm',status='old')
c     mutation matrix
      OPEN(unit=3,file='blosum.comm',status='old')
c     list of PDB templates
      OPEN(unit=9,file='list',status='old')
c     multiple sequence alignment for the target
      open(unit=11,file='msa.aln',status='old')
c     predicted secondary structure for the target
      open(unit=12,file='seq.dat',status='old')
c     sec.struct. and side-chain orientation independend
      OPEN(unit=22,FILE='pair.1',status='unknown')
c     sec.struct. and side-chain orientation specific
      OPEN(unit=23,FILE='pair.3',status='unknown')

ccccc apabl() probability for a contact between two residue-types
c     dependent on the side-chain orientation

c     initialize apabl()
      do i=-1,19
         do j=-1,19
            do it=1,3
               apabl(i,j,it)=0.
            end do
         end do
      end do

      do it=1,3
c     unit=2 <--> matrix3.comm
         read(2,*)
         read(2,*)
         do i=0,19
            read(2,17) ad3, (apabl(i,j,it),j=0,19)
         enddo
      end do
 17   FORMAT(A3,20F5.1)

ccccc read orientation independent probability, apa(i,j)
c     initialize apa()
      do i=-1,19
         do j=-1,19
            ap(i,j)=0.
         end do
      end do

c     unit=1 <--> matrix1.comm
      read(1,*)
      read(1,*)
      do i=0,19
         read(1,17) ad3, (ap(i,j),j=0,19)
      enddo

ccccc read BLOSUM mutation matrix
c     unit=3 <--> blosum.comm
      read(3,*)
      do i=0,19
         read(3,*)ad,(imut(i,j),j=0,19)
      end do

ccccc read MSA
c     unit=11 <--> msa.aln
      read(11,*)N_MSA0,Lch      !number of sequences in MSA, length
      if(Lch.gt.Lch_max)then
         write(*,*)'Query sequence length is too long!'
         stop
      endif
c     N_MSA current number of sequences in the multiple sequence
c     alignment with sequence identity above cutoff asim_cut
      N_MSA=1
c     seq(1,..) stores target sequence
      read(11,1)(seq(1,j),j=1,Lch)
 1    format(50a1)
      do i=2,N_MSA0
         read(11,1)(seq0(j),j=1,Lch) !target sequence
         isim=0
         do j=1,Lch
            if(seq0(j).eq.seq(1,j))isim=isim+1
         enddo
         asim=float(isim)/float(Lch)
         if(asim.gt.asim_cut)then
c     keep only those sequences above cutoff sequence identity
            N_MSA=N_MSA+1
            do j=1,Lch
               seq(N_MSA,j)=seq0(j)
            enddo
         endif
         if(N_MSA.ge.N_hom_max)goto 3
      end do
 3    continue

c     kcode(i,j) same as seq(i,j) but stores numbers (0-19)
c     instead of letters
      do i=1,N_MSA
         do j=1,Lch
            kcode(j,i)=-1
            do ires=0,19
               if(aa1(ires) .eq. seq(i,j))then
                  kcode(j,i)=ires
                  go to 18
               end if
            end do
 18         continue
         end do
      end do

c     jcode(i) stores the target sequence in numbers (0-19)
      do i=1,Lch
         jcode(i)=kcode(i,1)
      end do

c     nex(i) tells me there are nex(i) amino acids on type "i"
c     in the target sequence --> sequence composition.
c     You can access to the positions of all amino acids of
c     a given type "j" in the target sequence with ipos(i,j),
c     where i runs from 0 to next(i)
      mres=Lch
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
ccccc read MSA done

c     isec1(i) stores secondary structure codes for target sequence
      do i=1,Lch
         read(12,*)j,dum,isec1(i)
      enddo

ccccc read secondary structure codes of the templates,
c     and also their amino acid sequences

c     unit=9 <--> list
      read(9,*)N_stru
      DO 101 ik=1,N_stru
c     template(ik) stores name of a particular template
         read(9,*) template(ik)
c     im stores number characters in template name
         im=0
         do i=1,20
            if(template(ik)(i:i).eq.char(32))goto 51
            im=i
         enddo
 51      continue

         OPEN(unit=30,file=fnam2(1:ilib)//template(ik)(1:im)//'.cnt2')
c     Lch_str(ik) is seqence length of template
         read(30,*)Lch_str(ik)
         do i=1,Lch_str(ik)
c     isec2(i,ik) stores secondary structure codes for template ik
            read(30,*)ii,a1(ii),isec2(i,ik) !template sequences, sec
            read(30,*)
         enddo
         close(30)

         do i=1,Lch_str(ik)
            do j=0,19
               if(aa1(j).eq.a1(i))then
c     icode(i,ik) stores residue type (in number code, 0-19)
c     at position "i" of template "ik"
                  icode(i,ik)=j
                  goto 108
               endif
            enddo
c     assign GLY if a1(i) not matched with any of the 20
c     standar residue-types
            icode(i,ik)=0
 108        continue
         enddo
 101  continue

ccccc calculate pair1 ana pair3
      CALL EHMSTAT

ccccc output pair1
      mres=Lch
      do i=1,mres
         iseq=jcode(i)
         write(22,*)i,aa(iseq)
         WRITE(22,171) (am(i,j),j=1,mres)
         write(22,*)'================================='
      end do

ccccc output pair3
      do it=1,3
         do i=1,mres
            iseq=jcode(i)
            write(23,*)i,aa(iseq)
            WRITE(23,171) (an(i,j,it),j=1,mres)
            write(23,*)'================================='
         end do
      end do

 171  FORMAT(25F5.1)

      STOP
      END


      SUBROUTINE EHMSTAT
c     maximum length of PDB structure
      PARAMETER(Lch_str_max=2000)
c     maximum number of PDB structure
      PARAMETER(N_STR_MAX=20000)
c     maximum length of query sequence
      PARAMETER(Lch_max=2000)
c     maximum number of sequences in MSA
      parameter(N_hom_max=300)
      PARAMETER(IWIN=6)
      PARAMETER(IWIN2=13)

      character*20 template
      character*1 atmp
      common/names/template(N_str_max)
      common/c/N_stru
      dimension imap(-5:Lch_str_max,-5:Lch_str_max,3)
      common/seq/Lch,N_MSA,kcode(Lch_max,N_hom_max)
      common/scalepair/apabl(-1:19,-1:19,3)
      COMMON/Scalepair2/AN(Lch_max,Lch_max,3)
      dimension ileft(Lch_max,N_hom_max)
      common/seq2/imut(0:19,0:19)
      common/seq4/nex(0:19),ipos(Lch_max,0:19)
      dimension aobs(-4:Lch_max,Lch_max,3)
      DIMENSION IHIST(-8*IWIN2:16*IWIN2,Lch_max)
      dimension inat(Lch_max,N_hom_max)
      logical nogap(Lch_max,N_hom_max)
      common/independent/ap(-1:19,-1:19),am(Lch_max,Lch_max)
      common/second/isec1(Lch_max),isec2(Lch_str_max,N_str_max)
      common/stru/ICODE(LCH_STR_MAX,N_STR_MAX),Lch_str(N_str_max)
      character*300 fnam2
      common/libname/fnam2,ilib

      dimension icon(15,Lch_str_max),lor(15,Lch_str_max),nc(Lch_str_max)
      dimension an0(Lch_max,Lch_max,3)

ccccc pair.1=am(i,j) orientation independent
c     am(i,j) is a probability for residues at positions i and j to
c     interact, independent of the orientation of their side-chains

c     initialize am(i,j)
      do i=1,Lch                !Lch--Lch--length of sequence at MSA
         do j=1,Lch
            am(i,j)=0           !pair.1
         enddo
      enddo
  

c     N_MSA is number of homologues in the sequence-alignment
      do ihom=1,N_MSA
         do i=1,Lch
            do j=1,Lch
c     iseq is residue type at position "i" of homologue "ihom"
               iseq=kcode(i,ihom)
               jseq=kcode(j,ihom)
c     ap(i,j) probability from matrix1.comm
               am(i,j)=am(i,j)+ap(iseq,jseq) !20x20 contact matrix
               am(j,i)=am(j,i)+ap(iseq,jseq) !20x20 contact matrix
            enddo
         enddo
      enddo

c     normalize by number of possible homologue pairs
      do i=1,Lch
         do j=1,Lch
            am(i,j)=am(i,j)/(2*N_MSA)
         enddo
      enddo
ccccc pair.1 done


ccccc for pair.3=an0(i,j,it)  no fragments, global homology
c     an0(i,j,it) probability of residues at positions "i" and "j"
c     to interact if both are in the same orientation

c     initialize an0(i,j,it)
      do i=1,Lch
         do j=1,Lch
            do it=1,3
               an0(i,j,it)=0
            enddo
         enddo
      enddo

      do ihom=1,N_MSA
         do i=1,Lch
            do j=1,Lch
               do it=1,3
                  iseq=kcode(i,ihom)
                  jseq=kcode(j,ihom)
                  an0(i,j,it)=an0(i,j,it)+apabl(iseq,jseq,it)
                  an0(j,i,it)=an0(j,i,it)+apabl(iseq,jseq,it)
               enddo
            enddo
         enddo
      enddo
  
c     normalize by number of homologue pairs
      do i=1,Lch
         do j=1,Lch
            do it=1,3
               an0(i,j,it)=an0(i,j,it)/(2*N_MSA)
            enddo
         enddo
      enddo


ccccc pair.3=an(i,j,it) orientation dependent

      do i=1,Lch
         do j=1,Lch
            do it=1,3
               aobs(j,i,it)=0.  !pair3 element from one seqience in MSA
            enddo
         enddo
      enddo
 
ccccc pre-calculate   nogap(i,ihom), inat(i,ihom), and  ileft(i,ihom)
c     N_MSA is the number of homologous sequences.
      do 103 ihom=1,N_MSA

c     nogap(i,ihom) is false if there is a gap in any position within 
c     the range of positions [i-iwin,i+iwin] for homologue ihom
         do 150 i=1,Lch     
            nogap(i,ihom)=.true.
            do ii=-iwin,iwin    !iwin=6
               iv=i+ii
               if(iv .gt. 0 .and. iv.le.Lch)then
                  if(kcode(iv,ihom).eq. -1)then
                     nogap(i,ihom)=.false.
                     go to 150
                  endif
               endif
            enddo
 150     continue

ccccc calculate ileft(i,ihom). For a given homologue ihom and a given
c     position "i", we scan the library of PDB structures in search of
c     fragments of length iwin2 whose center residue has the same type
c     that residue at position "i" in homologue ihom. Then we align each
c     of these fragments to the fragment in ihom of length iwin2
c     centered at "i", and store all the scores. ileft(i,ihom) is a
c     lower bound score that guarantees that one of the selected
c     fragments in the PDB structures is homologous to the fragment in
c     ihom. (Note: the fragments have no gaps)
         do i=1,Lch
            if(nogap(i,ihom))then
               n_sim=0
               sim_a=0
               sim2_a=0
c     N_STRU is number of pdb structures
               do jk0=1,N_STRU
c     Lch_str(jk0) is length of structure jk0
                  do ij=1,Lch_str(jk0)
c     is1 is residue-type (0-19) at position "i" of homologue ihom
                     is1=kcode(i,ihom)
                     js1=icode(ij,jk0)
                     if(is1.eq.js1)then
                        isim=0
                        do ii=-iwin,iwin
                           iv=i+ii
                           jv=ij+ii
                           if(iv .gt. 0 .and. iv.le.Lch)then
                              if(jv .gt. 0 .and. jv.le.Lch_str(jk0))then
                                 is1=kcode(iv,ihom)
                                 js1=icode(jv,jk0) !structure
                                 isim=isim+imut(is1,js1) !BLOSUM matrix
                              endif
                           endif
                        enddo
                        n_sim=n_sim+1
c     sim_a is the average BLOSUM score between homologue and template in the [-win,iwin] window
                        sim_a=sim_a+isim
c     sim2_a is the second moment of the similarity measure
                        sim2_a=sim2_a+isim*isim
                     endif
                  enddo
               enddo
               if(n_sim.gt.0)then
c     finally compute average and standard deviation
                  sim_a=sim_a/float(n_sim)
                  sim2_a=sim2_a/float(n_sim)
                  dev=sqrt(sim2_a-sim_a*sim_a)
               else
                  sim_a=0
                  dev=0
               endif
c     ileft(i,ihom) is what we consider the lower cutoff for
c     similarity values that are above the similarity values
c     obtained from randomly pairing two sequence segments.
c     Different sequence positions "i" have different similarity
c     values.
               ileft(i,ihom)=int(sim_a+0.5)+3.5*int(dev+0.5)
c     -8*iwin2 is the absolute lower cutoff
               if(ileft(i,ihom).lt.-8*iwin2) ileft(i,ihom)=-8*iwin2
            endif
         enddo
ccccc finished calculating ileft


ccccc inat(i) is the score of the segment of length iwin2 centered
c     around position "i", that is, the log probability to of an
c     alignment with 100% identity. We normalize the probability (that is, we
c     substract the log probability, or score) by the lower cutoff of similarity
c     between two sequence segments ileft(i,ihom)
         do i=1,Lch
            if(nogap(i,ihom))then
               inat(i,ihom)=0
               do ii=-iwin,iwin
                  iv=i+ii
                  if(iv .gt. 0 .and. iv.le.Lch)then
                     jres=kcode(iv,ihom)
                     inat(i,ihom)=inat(i,ihom)+imut(jres,jres)
                  endif
               enddo
               inat(i,ihom)=inat(i,ihom)-ileft(i,ihom)
c     inat(i,ihom) less than one means a segment that likes to mutate
c     a lot, and we put an absolute value of 1 in this case, to
c     avoid negative values
               if(inat(i,ihom).lt.1)inat(i,ihom)=1
            endif
         enddo
c     go to next homologue
 103  continue
ccccc nogap(i,ihom), inat(i,ihom), and ileft(i,ihom) done


ccccc circle for each PDB structure
      do 1111 jk=1,N_STRU

c     im stores number characters in template name
         im=0
         do i=1,20
            if(template(jk)(i:i).eq.char(32))goto 51
            im=i
         enddo
 51      continue

         OPEN(unit=30,file=fnam2(1:ilib)//template(jk)(1:im)//'.cnt2')
c     nn is lenght of template
         read(30,*)nn
         do i=1,nn
c     nc(ii) is number of contacts that residue at position "i"
c     makes with rest of resides in the template.
c     icon(j,ii), where j=1,..,nc(ii) are the position of the
c     residues contacting 'ii'
            read(30,*)ii,atmp,itmp,nc(ii),(icon(j,ii),j=1,nc(ii)) !contact
c     lor(j,ii), where j=1,..,nct==nc(ii) are the orientations (1-3)
c     of each of the nct residues contacting "ii"
            read(30,*)nct,(lor(j,ii),j=1,nct)
         enddo
         close(30)

ccccc imap(i,j,it) is the contact map but with an extra paramenter,
c     the orientation "it". imap(i,j,it) is 1 if there's a contact.

c     initialize imap(i,j,it)
         do i=1,Lch_str(jk)+1
            do j=i,Lch_str(jk)+1
               do it=1,3
                  imap(i,j,it)=0
                  imap(j,i,it)=0
               enddo
            enddo
         enddo

         do ij=1,Lch_str(jk)
c     if residue at position ij is doing nc(ij) with rest of
c     template, then..
            if(nc(ij).gt.0)then
               do kk=1,nc(ij)
c     jj is the position of one of the residues contacting ij
                  jj=icon(kk,ij)
                  it=lor(kk,ij)
                  imap(ij,jj,it)=1
               enddo
            endif
         enddo
         
ccccc calculate aobs(i,j,1-3), weighted number of contacts between (i,j)
         do 333 ij=1,Lch_str(jk)
c     iseq1 is residue-type (0-19) at position "ij" of template "jk"
            iseq1=icode(ij,jk)
c     go to next position in the template if target sequence contains no
c     residue of type iseq1
            if(nex(iseq1).le.0)goto 333
c     go through every contact of residue "ij"
            do 444 kk=1,nc(ij)
c     jj is residue-type (0-19) at position "kk" of template "jk"
               jj=icon(kk,ij)
c     go to next contacting residue if we have a local contact between
c     contacting positions "ij" and "jj"
               if(iabs(jj-ij) .le.iwin)goto 444
c     iseq2 is residue-type (0-19) at position "jj" of template "jk"
               iseq2=icode(jj,jk)
c     go to next contacting residue if target sequence contains no
c     residue of type iseq2
               if(nex(iseq2).le.0)goto 444
c     go through all the sequence positions in the target sequence such
c     that a residue of type iseq1 exists
               do 555 i1=1,nex(iseq1)
c     "i" is the sequence position in the target of one residue that is
c     of type iseq1
                  i=ipos(i1,iseq1)
c     go to next sequence position in the target if the secondary
c     structure at position "i" in the target is different than the
c     secondary structure at position "ij" in the template
                  if(isec1(i).ne.isec2(ij,jk))goto 555
c     go through all the sequence positions such that a residue of type
c     iseq2 exists in the target sequence
                  do 666 j1=1,nex(iseq2) !circle on target
c     "j" is  the sequence position in the target of one residue that is
c     of type iseq2
                     j=ipos(j1,iseq2)
c     go to next sequence position in the target if the secondary
c     structure at position "j" in the target is different than the
c     secondary structure at position "jj" in the template
                     if(isec1(j).ne.isec2(jj,jk))goto 666
c     go to next sequence position in the target if positions "i" and
c     "j" are too close in sequence separation (note that "i" and "j" do
c     not need to make a contact. It is corresponding residues "ij" and
c     "jj" in the template the ones that actually do)
                     if(iabs(i-j).le.iwin)goto 666

ccccc circle for each sequence in MSA
                     do 2222 ihom=1,n_msa
                        if(.not.nogap(i,ihom))goto 2222
                        if(.not.nogap(j,ihom))goto 2222

ccccc calculate similarity between segment of homologue (centered at
c     "i") and segment of template (centered at "ij")
                        isim1=0
                        do ii=-iwin,iwin
                           iv=i+ii
                           jv=ij+ii
                           if(iv .gt.0 .and. iv.le.Lch)then
                              if(jv .gt. 0 .and. jv.le.Lch_str(jk))then
                                 is1=kcode(iv,ihom)
                                 js1=icode(jv,jk)
                                 isim1=isim1+imut(is1,js1)
                              end if
                           end if
                        end do

ccccc calculate similarity between segment of homologue (centered at
c     "j") and segment of template (centered at "jj")
                        isim2=0
                        do ii=-iwin,iwin
                           iv=j+II
                           jv=jj+ii
                           if(iv .gt.0 .and. iv.le.Lch)then
                              if(jv .gt.0 .and. jv.le.Lch_str(jk))then
                                 is1=kcode(iv,ihom)
                                 js1=icode(jv,jk)
                                 isim2=isim2+imut(is1,js1)
                              end if
                           end if
                        end do

c     go to next template if the similarities are negative
                        if(isim1.le.0)goto 2222
                        if(isim2.le.0)goto 2222
                        if(isim1.ge.ileft(i,ihom) .or.
     $                       isim2.ge.ileft(j,ihom)    )then
                           asim1=isim1
                           asim2=isim2
c     sf is the product of the two similarities in "units" of inat, and
c     inat was log probability of the score of the unmutated sequence
c     minus the log probability of the lower cutoff for similarity
c     (ilef)
                           sf=asim1*asim2/(inat(i,ihom)*inat(j,ihom))
                           sf=0.2+sqrt(sf)
c     open a segment of 11 residues in the homologue centered in
c     position "i" and in the template centered at "ij"
                           do 777 iv=-5,5
                              iz=i+iv
                              ijz=ij+iv
                              if(iz.lt.1)goto 777
                              if(ijz.lt.1)goto 777
                              if(iz.gt.Lch)goto 777
                              if(ijz.gt.Lch_str(jk))goto 777
c     iz1 is the residue-type (0-19) at position iz in the homologue
                              iz1=kcode(iz,ihom)
c     ijz1 is the residue-type (0-19) at position ijz in the homologue
                              ijz1=icode(ijz,jk)
c     if the two residue-types are not evolutionary related, skip to
c     next position in the segment
                              if(imut(iz1,ijz1).le.0)goto 777
c     open a segment of 11 residues in the homologue centered in
c     position "j" and in the template centered at "jj"
                              do 888 jv=-5,5
                                 jz=j+jv
                                 jjz=jj+jv
                                 if(jz.lt.1)goto 888
                                 if(jz.gt.Lch)goto 888
                                 if(jjz.lt.1)goto 888
                                 if(jjz.gt.Lch_str(jk))goto 888
c     jz1 is the residue-type (0-19) at position jz in the homologue
                                 jz1=kcode(jz,ihom)
c     jjz1 is the residue-type (0-19) at position jjz in the homologue
                                 jjz1=icode(jjz,jk)
c     if the two residue-types are not evolutionary related, skip to
c     next position in the segment
                                 if(imut(jz1,jjz1).le.0)goto 888
c     if residues at positions iz and jz are in contact in template jk,
c     then add weight sf to the weigthed contact map aobs(.., .., ..)
                                 do it=1,3
                                    aobs(iz,jz,it)=aobs(iz,jz,it) +
     $                                   sf*imap(ijz,jjz,it)
                                    aobs(jz,iz,it)=aobs(jz,iz,it) +
     $                                   sf*imap(ijz,jjz,it)
                                 enddo
c     go to next sequence position in the segment centered at "j"/"jj"
 888                          continue
c     go to next sequence position in the segment centered at "i"/"ii"
 777                       continue
                        endif
c     go to next homologue
 2222                continue
c     go to next sequence position in the target sequence such that the
c     residue type is the same that the type of residue at position "jj"
c     in the template
 666              continue
c     go to next sequence position in the target sequence such that the
c     residue type is the same that the type of residue at position "ij"
c     in the template
 555           continue
c     go to next residue contacting residue at position "ij"
 444        continue
c     go to next position along the template sequence
 333     continue
c     go to next template
 1111 continue

ccccc an(i,j,it) normalization, atot is the total number of weighted
c     contacts for all orientations
      atot=0.
      do i=1,Lch
         do j=1,Lch
            do it=1,3
               atot=atot+aobs(i,j,it)
            enddo
         enddo
      enddo
      anorm=atot/(3*Lch*Lch)
     
      do i=1,Lch
         do j=1,Lch
            do it=1,3
c     
               if(aobs(i,j,it).gt.0)then
                  an(i,j,it)=-alog(aobs(i,j,it)/anorm)
                  an(j,i,it)=-alog(aobs(i,j,it)/anorm)
               else
c     if aobs(.., .., ..) negative, then the contact potential is
c     secondary structure independent.
                  an(i,j,it)=an0(i,j,it)
                  an(j,i,it)=an0(j,i,it)
               endif
            enddo
         enddo
      enddo
ccccc pair.3 done

      return
      end


