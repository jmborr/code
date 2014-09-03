ccccc list of subroutines
c     
c     seq_id(character target(ndim), character(ndim) template, integer
c     gap_open, integer gap_extn, real seqid)
c
ccccc 

ccccc
c     This is a standard Needleman-Wunsch dynamic program (by Y. Zhang 2005)
c     1. Count multiple-gap.
c     2. The gap penality W(k)=Go+Ge*k1+Go+Ge*k2 if gap open on both sequences
c
c     Input: score(i,j), gap_open, gap_extn
c     Output: j2i(j)
c     idir(i,j)=1,2,3, from diagonal, horizontal, vertical
c     val(i,j) is the cumulative score of (i,j)
c
c     suggested: gap_open=-11,gap_extn=-1
c
c     f2py -c -m align_yang_f77 align_yang.f
ccccc 

      function seq_id(target,template,gap_open,gap_extn)
Cf2py integer optional,intent(in) :: gap_open = -11
Cf2py integer optional,intent(in) :: gap_extn = -1
c     maximum number of residues in either sequence
      parameter(ndim=3000)
      common/dpc/score(ndim,ndim),j2i(ndim)
      common/lengths/nseq1,nseq2
      common/seq12/seq1(ndim),seq2(ndim)
      character target*3000,template*3000
      integer gap_open,gap_extn
      real seq_id

ccccc find length of sequences nseq1,nseq2
      call getlengths(target,template)
ccccc initialize seq1 and seq2, which are sequence in numeric code
      call init_seq12(target,template)
ccccc initialize blosum mattrix
      call matrix
ccccc do global alignment
      call dp(gap_open,gap_extn)
ccccc calculate sequence identity
      L_id=0
      L_ali=0
      do j=1,nseq2
         if(j2i(j).gt.0)then
            i=j2i(j)
            L_ali=L_ali+1
            if(seq1(i).eq.seq2(j)) then 
               L_id=L_id+1
            endif
         endif
      enddo
      seq_id=(L_id*1.0)/L_ali
      return
      end
ccccc seq_id finished

ccccc return many things, not only sequence identity
c     nid number of identical residues in the aligment
c     nali number of aligned residues
c     return the tuple (nid,nali)
      subroutine align(target,template,gap_open,gap_extn,nid,nali)      
Cf2py integer optional,intent(in) :: gap_open = -11
Cf2py integer optional,intent(in) :: gap_extn = -1
Cf2py integer intent(out) :: nid
Cf2py integer intent(out) :: nali
c     maximum number of residues in either sequence
      parameter(ndim=3000)
      common/dpc/score(ndim,ndim),j2i(ndim)
      common/lengths/nseq1,nseq2
      common/seq12/seq1(ndim),seq2(ndim)
      character target*3000,template*3000
      integer gap_open,gap_extn

ccccc find length of sequences nseq1,nseq2
      call getlengths(target,template)
ccccc initialize seq1 and seq2, which are sequence in numeric code
      call init_seq12(target,template)
ccccc initialize blosum mattrix
      call matrix
ccccc do global alignment
      call dp(gap_open,gap_extn)
ccccc calculate several things
      L_id=0
      L_ali=0
      do j=1,nseq2
         if(j2i(j).gt.0)then
            i=j2i(j)
            L_ali=L_ali+1
            if(seq1(i).eq.seq2(j)) then 
               L_id=L_id+1
            endif
         endif
      enddo
      nid=L_id
      nali=L_ali
      return
      end
ccccc align finished

ccccc get sequence lengths and nseq1 and nseq2
      subroutine getlengths(target,template)
      common/lengths/nseq1,nseq2
      character target*3000,template*3000
      nseq1=0
      do i=1,3000
         if(target(i:i).eq.' ') goto 18
         nseq1=i
      enddo
      nseq2=0
 18   do i=1,3000
         if(template(i:i).eq.' ') goto 19
         nseq2=i
      enddo
 19   continue
      return
      end
ccccc getlengths finished

ccccc initialize seq1 and seq2
      subroutine init_seq12(target,template)
      parameter(ndim=3000)
      parameter(naa=23)
      common/lengths/nseq1,nseq2
      common/seq12/seq1(ndim),seq2(ndim)
      character target*3000,template*3000
      character seqw(naa)
      data seqw/'A','R','N','D','C','Q','E','G','H','I','L','K',
     $     'M','F','P','S','T','W','Y','V','B','Z','X'/
      do 96 i=1,nseq1
         do j=1,naa
            if(target(i:i).eq.seqw(j)) then
               seq1(i)=j
               goto 96
            endif
         enddo
 96   continue
      do 97 i=1,nseq2
         do j=1,naa
            if(template(i:i).eq.seqw(j)) then
               seq2(i)=j
               goto 97
            endif
         enddo
 97   continue    
      return
      end
ccccc init_seq12 finished


ccccc dynamic programming
      subroutine dp(gap_open,gap_extn)
      parameter(ndim=3000)
      parameter(naa=23)
      common/lengths/nseq1,nseq2
      common/dpc/score(ndim,ndim),j2i(ndim)
      common/matra/imut(naa,naa)
      common/seq12/seq1(ndim),seq2(ndim)

      dimension val(0:ndim,0:ndim),idir(0:ndim,0:ndim)
      dimension jpV(0:ndim,0:ndim),jpH(0:ndim,0:ndim)
      dimension preV(0:ndim,0:ndim),preH(0:ndim,0:ndim)

      real D,V,H
      integer gap_open,gap_extn

ccccc initialize pair-score from mutation matrix
      do i=1,nseq1
         do j=1,nseq2
            score(i,j)=imut(seq1(i),seq2(j))
         enddo
      enddo

ccccc initialize different matrixes
      val(0,0)=0.0
      do i=1,nseq1
         val(i,0)=0.0
         idir(i,0)=0
         jpV(i,0)=1
         jpH(i,0)=i
         preV(i,0)=-1000.0
      enddo
      do j=1,nseq2
         val(0,j)=0.0
         idir(0,j)=0
         jpV(0,j)=j
         jpH(0,j)=1
         preH(0,j)=-1000.0
      enddo

ccccc perform dynamic programming
      do 111 j=1,nseq2
         do 222 i=1,nseq1
c     D=VAL(i-1,j-1)+SCORE(i,j)
c     from diagonal, val(i,j) is val(i-1,j-1)
            D=val(i-1,j-1)+score(i,j) 
c     H=H+gap_open
            jpH(i,j)=1
c     gap_open from both D and V
            val1=val(i-1,j)+gap_open 
c     gap_extn from horizontal
            val2=preH(i-1,j)+gap_extn 
c     last step from D or V
            if(val1.gt.val2) then 
               H=val1
c     last step from H
            else                
               H=val2
c     record long-gap
               if(i.gt.1)jpH(i,j)=jpH(i-1,j)+1 
            endif
c     V=V+gap_open
            jpV(i,j)=1
            val1=val(i,j-1)+gap_open
            val2=preV(i,j-1)+gap_extn
            if(val1.gt.val2) then
               V=val1
            else
               V=val2
               if(j.gt.1)jpV(i,j)=jpV(i,j-1)+1
            endif
c     unaccepted H
            preH(i,j)=H        
c     unaccepted V
            preV(i,j)=V         
            
            if(D.gt.H.and.D.gt.V)then
               idir(i,j)=1
               val(i,j)=D
            elseif(H.gt.V)then
               idir(i,j)=2
               val(i,j)=H
            else
               idir(i,j)=3
               val(i,j)=V
            endif
 222     continue
 111  continue
c     alignment score
      score0=val(nseq1,nseq2)   

c     tracing back the pathway:
      do j=1,nseq2
c     all are not aligned
         j2i(j)=-1              
      enddo
      i=nseq1
      j=nseq2
      do while(i.gt.0.and.j.gt.0)
c     from diagonal
         if(idir(i,j).eq.1)then 
            j2i(j)=i
            i=i-1
            j=j-1
c     from horizonal
         elseif(idir(i,j).eq.2)then 
            do me=1,jpH(i,j)
               if(i.gt.0) then
                  i=i-1
               endif
            enddo
         else
            do me=1,jpV(i,j)
               if(j.gt.0) then
                  j=j-1
               endif
            enddo
         endif
      enddo
      return
      end
ccccc dp finished


ccccc read matrix
      subroutine matrix
c     number of amino acid
      parameter(naa=23) 
c     b,z,x are additional
      common/matra/imut(naa,naa) 

c     folowing from BLOSUM62 used in BLAST:
      imut(1,1)=4
      imut(1,2)=-1
      imut(1,3)=-2
      imut(1,4)=-2
      imut(1,5)=0
      imut(1,6)=-1
      imut(1,7)=-1
      imut(1,8)=0
      imut(1,9)=-2
      imut(1,10)=-1
      imut(1,11)=-1
      imut(1,12)=-1
      imut(1,13)=-1
      imut(1,14)=-2
      imut(1,15)=-1
      imut(1,16)=1
      imut(1,17)=0
      imut(1,18)=-3
      imut(1,19)=-2
      imut(1,20)=0
      imut(1,21)=-2
      imut(1,22)=-1
      imut(1,23)=0
      imut(2,1)=-1
      imut(2,2)=5
      imut(2,3)=0
      imut(2,4)=-2
      imut(2,5)=-3
      imut(2,6)=1
      imut(2,7)=0
      imut(2,8)=-2
      imut(2,9)=0
      imut(2,10)=-3
      imut(2,11)=-2
      imut(2,12)=2
      imut(2,13)=-1
      imut(2,14)=-3
      imut(2,15)=-2
      imut(2,16)=-1
      imut(2,17)=-1
      imut(2,18)=-3
      imut(2,19)=-2
      imut(2,20)=-3
      imut(2,21)=-1
      imut(2,22)=0
      imut(2,23)=-1
      imut(3,1)=-2
      imut(3,2)=0
      imut(3,3)=6
      imut(3,4)=1
      imut(3,5)=-3
      imut(3,6)=0
      imut(3,7)=0
      imut(3,8)=0
      imut(3,9)=1
      imut(3,10)=-3
      imut(3,11)=-3
      imut(3,12)=0
      imut(3,13)=-2
      imut(3,14)=-3
      imut(3,15)=-2
      imut(3,16)=1
      imut(3,17)=0
      imut(3,18)=-4
      imut(3,19)=-2
      imut(3,20)=-3
      imut(3,21)=3
      imut(3,22)=0
      imut(3,23)=-1
      imut(4,1)=-2
      imut(4,2)=-2
      imut(4,3)=1
      imut(4,4)=6
      imut(4,5)=-3
      imut(4,6)=0
      imut(4,7)=2
      imut(4,8)=-1
      imut(4,9)=-1
      imut(4,10)=-3
      imut(4,11)=-4
      imut(4,12)=-1
      imut(4,13)=-3
      imut(4,14)=-3
      imut(4,15)=-1
      imut(4,16)=0
      imut(4,17)=-1
      imut(4,18)=-4
      imut(4,19)=-3
      imut(4,20)=-3
      imut(4,21)=4
      imut(4,22)=1
      imut(4,23)=-1
      imut(5,1)=0
      imut(5,2)=-3
      imut(5,3)=-3
      imut(5,4)=-3
      imut(5,5)=9
      imut(5,6)=-3
      imut(5,7)=-4
      imut(5,8)=-3
      imut(5,9)=-3
      imut(5,10)=-1
      imut(5,11)=-1
      imut(5,12)=-3
      imut(5,13)=-1
      imut(5,14)=-2
      imut(5,15)=-3
      imut(5,16)=-1
      imut(5,17)=-1
      imut(5,18)=-2
      imut(5,19)=-2
      imut(5,20)=-1
      imut(5,21)=-3
      imut(5,22)=-3
      imut(5,23)=-2
      imut(6,1)=-1
      imut(6,2)=1
      imut(6,3)=0
      imut(6,4)=0
      imut(6,5)=-3
      imut(6,6)=5
      imut(6,7)=2
      imut(6,8)=-2
      imut(6,9)=0
      imut(6,10)=-3
      imut(6,11)=-2
      imut(6,12)=1
      imut(6,13)=0
      imut(6,14)=-3
      imut(6,15)=-1
      imut(6,16)=0
      imut(6,17)=-1
      imut(6,18)=-2
      imut(6,19)=-1
      imut(6,20)=-2
      imut(6,21)=0
      imut(6,22)=3
      imut(6,23)=-1
      imut(7,1)=-1
      imut(7,2)=0
      imut(7,3)=0
      imut(7,4)=2
      imut(7,5)=-4
      imut(7,6)=2
      imut(7,7)=5
      imut(7,8)=-2
      imut(7,9)=0
      imut(7,10)=-3
      imut(7,11)=-3
      imut(7,12)=1
      imut(7,13)=-2
      imut(7,14)=-3
      imut(7,15)=-1
      imut(7,16)=0
      imut(7,17)=-1
      imut(7,18)=-3
      imut(7,19)=-2
      imut(7,20)=-2
      imut(7,21)=1
      imut(7,22)=4
      imut(7,23)=-1
      imut(8,1)=0
      imut(8,2)=-2
      imut(8,3)=0
      imut(8,4)=-1
      imut(8,5)=-3
      imut(8,6)=-2
      imut(8,7)=-2
      imut(8,8)=6
      imut(8,9)=-2
      imut(8,10)=-4
      imut(8,11)=-4
      imut(8,12)=-2
      imut(8,13)=-3
      imut(8,14)=-3
      imut(8,15)=-2
      imut(8,16)=0
      imut(8,17)=-2
      imut(8,18)=-2
      imut(8,19)=-3
      imut(8,20)=-3
      imut(8,21)=-1
      imut(8,22)=-2
      imut(8,23)=-1
      imut(9,1)=-2
      imut(9,2)=0
      imut(9,3)=1
      imut(9,4)=-1
      imut(9,5)=-3
      imut(9,6)=0
      imut(9,7)=0
      imut(9,8)=-2
      imut(9,9)=8
      imut(9,10)=-3
      imut(9,11)=-3
      imut(9,12)=-1
      imut(9,13)=-2
      imut(9,14)=-1
      imut(9,15)=-2
      imut(9,16)=-1
      imut(9,17)=-2
      imut(9,18)=-2
      imut(9,19)=2
      imut(9,20)=-3
      imut(9,21)=0
      imut(9,22)=0
      imut(9,23)=-1
      imut(10,1)=-1
      imut(10,2)=-3
      imut(10,3)=-3
      imut(10,4)=-3
      imut(10,5)=-1
      imut(10,6)=-3
      imut(10,7)=-3
      imut(10,8)=-4
      imut(10,9)=-3
      imut(10,10)=4
      imut(10,11)=2
      imut(10,12)=-3
      imut(10,13)=1
      imut(10,14)=0
      imut(10,15)=-3
      imut(10,16)=-2
      imut(10,17)=-1
      imut(10,18)=-3
      imut(10,19)=-1
      imut(10,20)=3
      imut(10,21)=-3
      imut(10,22)=-3
      imut(10,23)=-1
      imut(11,1)=-1
      imut(11,2)=-2
      imut(11,3)=-3
      imut(11,4)=-4
      imut(11,5)=-1
      imut(11,6)=-2
      imut(11,7)=-3
      imut(11,8)=-4
      imut(11,9)=-3
      imut(11,10)=2
      imut(11,11)=4
      imut(11,12)=-2
      imut(11,13)=2
      imut(11,14)=0
      imut(11,15)=-3
      imut(11,16)=-2
      imut(11,17)=-1
      imut(11,18)=-2
      imut(11,19)=-1
      imut(11,20)=1
      imut(11,21)=-4
      imut(11,22)=-3
      imut(11,23)=-1
      imut(12,1)=-1
      imut(12,2)=2
      imut(12,3)=0
      imut(12,4)=-1
      imut(12,5)=-3
      imut(12,6)=1
      imut(12,7)=1
      imut(12,8)=-2
      imut(12,9)=-1
      imut(12,10)=-3
      imut(12,11)=-2
      imut(12,12)=5
      imut(12,13)=-1
      imut(12,14)=-3
      imut(12,15)=-1
      imut(12,16)=0
      imut(12,17)=-1
      imut(12,18)=-3
      imut(12,19)=-2
      imut(12,20)=-2
      imut(12,21)=0
      imut(12,22)=1
      imut(12,23)=-1
      imut(13,1)=-1
      imut(13,2)=-1
      imut(13,3)=-2
      imut(13,4)=-3
      imut(13,5)=-1
      imut(13,6)=0
      imut(13,7)=-2
      imut(13,8)=-3
      imut(13,9)=-2
      imut(13,10)=1
      imut(13,11)=2
      imut(13,12)=-1
      imut(13,13)=5
      imut(13,14)=0
      imut(13,15)=-2
      imut(13,16)=-1
      imut(13,17)=-1
      imut(13,18)=-1
      imut(13,19)=-1
      imut(13,20)=1
      imut(13,21)=-3
      imut(13,22)=-1
      imut(13,23)=-1
      imut(14,1)=-2
      imut(14,2)=-3
      imut(14,3)=-3
      imut(14,4)=-3
      imut(14,5)=-2
      imut(14,6)=-3
      imut(14,7)=-3
      imut(14,8)=-3
      imut(14,9)=-1
      imut(14,10)=0
      imut(14,11)=0
      imut(14,12)=-3
      imut(14,13)=0
      imut(14,14)=6
      imut(14,15)=-4
      imut(14,16)=-2
      imut(14,17)=-2
      imut(14,18)=1
      imut(14,19)=3
      imut(14,20)=-1
      imut(14,21)=-3
      imut(14,22)=-3
      imut(14,23)=-1
      imut(15,1)=-1
      imut(15,2)=-2
      imut(15,3)=-2
      imut(15,4)=-1
      imut(15,5)=-3
      imut(15,6)=-1
      imut(15,7)=-1
      imut(15,8)=-2
      imut(15,9)=-2
      imut(15,10)=-3
      imut(15,11)=-3
      imut(15,12)=-1
      imut(15,13)=-2
      imut(15,14)=-4
      imut(15,15)=7
      imut(15,16)=-1
      imut(15,17)=-1
      imut(15,18)=-4
      imut(15,19)=-3
      imut(15,20)=-2
      imut(15,21)=-2
      imut(15,22)=-1
      imut(15,23)=-2
      imut(16,1)=1
      imut(16,2)=-1
      imut(16,3)=1
      imut(16,4)=0
      imut(16,5)=-1
      imut(16,6)=0
      imut(16,7)=0
      imut(16,8)=0
      imut(16,9)=-1
      imut(16,10)=-2
      imut(16,11)=-2
      imut(16,12)=0
      imut(16,13)=-1
      imut(16,14)=-2
      imut(16,15)=-1
      imut(16,16)=4
      imut(16,17)=1
      imut(16,18)=-3
      imut(16,19)=-2
      imut(16,20)=-2
      imut(16,21)=0
      imut(16,22)=0
      imut(16,23)=0
      imut(17,1)=0
      imut(17,2)=-1
      imut(17,3)=0
      imut(17,4)=-1
      imut(17,5)=-1
      imut(17,6)=-1
      imut(17,7)=-1
      imut(17,8)=-2
      imut(17,9)=-2
      imut(17,10)=-1
      imut(17,11)=-1
      imut(17,12)=-1
      imut(17,13)=-1
      imut(17,14)=-2
      imut(17,15)=-1
      imut(17,16)=1
      imut(17,17)=5
      imut(17,18)=-2
      imut(17,19)=-2
      imut(17,20)=0
      imut(17,21)=-1
      imut(17,22)=-1
      imut(17,23)=0
      imut(18,1)=-3
      imut(18,2)=-3
      imut(18,3)=-4
      imut(18,4)=-4
      imut(18,5)=-2
      imut(18,6)=-2
      imut(18,7)=-3
      imut(18,8)=-2
      imut(18,9)=-2
      imut(18,10)=-3
      imut(18,11)=-2
      imut(18,12)=-3
      imut(18,13)=-1
      imut(18,14)=1
      imut(18,15)=-4
      imut(18,16)=-3
      imut(18,17)=-2
      imut(18,18)=11
      imut(18,19)=2
      imut(18,20)=-3
      imut(18,21)=-4
      imut(18,22)=-3
      imut(18,23)=-2
      imut(19,1)=-2
      imut(19,2)=-2
      imut(19,3)=-2
      imut(19,4)=-3
      imut(19,5)=-2
      imut(19,6)=-1
      imut(19,7)=-2
      imut(19,8)=-3
      imut(19,9)=2
      imut(19,10)=-1
      imut(19,11)=-1
      imut(19,12)=-2
      imut(19,13)=-1
      imut(19,14)=3
      imut(19,15)=-3
      imut(19,16)=-2
      imut(19,17)=-2
      imut(19,18)=2
      imut(19,19)=7
      imut(19,20)=-1
      imut(19,21)=-3
      imut(19,22)=-2
      imut(19,23)=-1
      imut(20,1)=0
      imut(20,2)=-3
      imut(20,3)=-3
      imut(20,4)=-3
      imut(20,5)=-1
      imut(20,6)=-2
      imut(20,7)=-2
      imut(20,8)=-3
      imut(20,9)=-3
      imut(20,10)=3
      imut(20,11)=1
      imut(20,12)=-2
      imut(20,13)=1
      imut(20,14)=-1
      imut(20,15)=-2
      imut(20,16)=-2
      imut(20,17)=0
      imut(20,18)=-3
      imut(20,19)=-1
      imut(20,20)=4
      imut(20,21)=-3
      imut(20,22)=-2
      imut(20,23)=-1
      imut(21,1)=-2
      imut(21,2)=-1
      imut(21,3)=3
      imut(21,4)=4
      imut(21,5)=-3
      imut(21,6)=0
      imut(21,7)=1
      imut(21,8)=-1
      imut(21,9)=0
      imut(21,10)=-3
      imut(21,11)=-4
      imut(21,12)=0
      imut(21,13)=-3
      imut(21,14)=-3
      imut(21,15)=-2
      imut(21,16)=0
      imut(21,17)=-1
      imut(21,18)=-4
      imut(21,19)=-3
      imut(21,20)=-3
      imut(21,21)=4
      imut(21,22)=1
      imut(21,23)=-1
      imut(22,1)=-1
      imut(22,2)=0
      imut(22,3)=0
      imut(22,4)=1
      imut(22,5)=-3
      imut(22,6)=3
      imut(22,7)=4
      imut(22,8)=-2
      imut(22,9)=0
      imut(22,10)=-3
      imut(22,11)=-3
      imut(22,12)=1
      imut(22,13)=-1
      imut(22,14)=-3
      imut(22,15)=-1
      imut(22,16)=0
      imut(22,17)=-1
      imut(22,18)=-3
      imut(22,19)=-2
      imut(22,20)=-2
      imut(22,21)=1
      imut(22,22)=4
      imut(22,23)=-1
      imut(23,1)=0
      imut(23,2)=-1
      imut(23,3)=-1
      imut(23,4)=-1
      imut(23,5)=-2
      imut(23,6)=-1
      imut(23,7)=-1
      imut(23,8)=-1
      imut(23,9)=-1
      imut(23,10)=-1
      imut(23,11)=-1
      imut(23,12)=-1
      imut(23,13)=-1
      imut(23,14)=-1
      imut(23,15)=-2
      imut(23,16)=0
      imut(23,17)=0
      imut(23,18)=-2
      imut(23,19)=-1
      imut(23,20)=-1
      imut(23,21)=-1
      imut(23,22)=-1
      imut(23,23)=-1
        
      return
      end
