ccccc 
c     TM-align version 0.1. A small bug of two-point superposition was 
c     fixed at June 1, 2005, compared with previous version.
c
c     This program identifies the best alignment of two protein
c     structures to give the best TM-score. By default, TM-score is
c     normalized by the second protein. The program can be freely copied
c     or modified or redistributed.
c
c     Reference:
c     Yang Zhang, Jeffrey Skolnick, Nucl. Acid Res. 2005 33: 2303-9
c     (For comments, please email to: zhang6@buffalo.edu)      

      function tmalign_yang(fnam1,fnam2,alignf)
Cf2py character*512 optional,intent(in) :: alignf = ''
      parameter(nmax=3000)
      parameter(nmax2=6000)
      parameter(nmax3=9000)
      
      common/backbone/xa(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/length/nseq1,nseq2,n8_al
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002

      character*512 alignf
      character*200 fnam1,fnam2,pdb(100),outname
      character*3 aa(-1:20),aanam,ss1(nmax),ss2(nmax)
      character*100 s,du
      character seq1(0:nmax),seq2(0:nmax)
      character aseq1(nmax2),aseq2(nmax2),aseq3(nmax2)

      dimension xx(nmax),yy(nmax),zz(nmax)
      dimension m1(nmax),m2(nmax)
      dimension mm1(nmax),mm2(nmax)
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/init/invmap_i(nmax)

      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d8/d8
      
ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms
      data w /nmax*1.0/
ccc   

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     $     'PRO','MET','ASP','ASN','LEU',
     $     'LYS','GLU','GLN','ARG',
     $     'HIS','PHE','TYR','TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C','V','T','I',
     $     'P','M','D','N','L','K','E','Q','R',
     $     'H','F','Y','W','C'/

ccccc options 
c     decided output
      m_out=-1                  
c     fixed length-scale only for output
      m_fix=-1                  
c     using average length
      m_ave=-1                  
c     diminum d0 for search
      m_d0_min=-1          
c     given d0 for both search and output
      m_d0=-1                   

      pdb(1)=fnam1
      pdb(2)=fnam2
ccccc read data from first CA file (template) seq1(i), xa(1,i,j),
c     ss1(i), and nseq1
      open(unit=10,file=pdb(1),status='old')
      i=0
      do while (.true.)
         read(10,9001,end=1010) s
         if(i.gt.0.and.s(1:3).eq.'TER')goto 1010
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or
     &           .s(13:16).eq.'  CA')then
               i=i+1
c     coordinates of template in xa(j,i,0)
               read(s,9000)du,aanam,du,mm1(i),du,
     $              xa(1,i,0),xa(2,i,0),xa(3,i,0)
               do j=-1,20
c     residue types of template in seq1(i), in numeric code
                  if(aanam.eq.aa(j))seq1(i)=slc(j)
               enddo
c     residue types of template in ss1(i), in one-letter code
               ss1(i)=aanam
               if(i.ge.nmax)goto 1010
            endif
         endif
      enddo
 1010 continue
 9000 format(A17,A3,A2,i4,A4,3F8.3)
 9001 format(A100)
      close(10)
      nseq1=i


cccc  read data from the second CA file:
      open(unit=10,file=pdb(2),status='old')
      i=0
      do while (.true.)
         read(10,9001,end=1011) s
         if(i.gt.0.and.s(1:3).eq.'TER')goto 1011
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.
     &           s(13:16).eq.'  CA')then
               i=i+1
               read(s,9000)du,aanam,du,mm2(i),du,
     $              xa(1,i,1),xa(2,i,1),xa(3,i,1)
               do j=-1,20
                  if(aanam.eq.aa(j))seq2(i)=slc(j)
               enddo
               ss2(i)=aanam
               if(i.ge.nmax)goto 1011
            endif
         endif
      enddo
 1011 continue
      close(10)
      nseq2=i


ccccc Scale of TM-score
      d0_min=0.5
c     if we passed a minimum d0 cutoff, then..
      if(m_d0_min.eq.1)then
         d0_min=d0_min_input    
      endif

c     length for defining TMscore in search based on smaller sequence
      anseq_min=min(nseq1,nseq2)
      anseq=anseq_min           

c     remove pairs with dis>d8 during search & final
      d8=1.5*anseq_min**0.3+3.5 

c     use formula for d0 when length of smallest sequence bigger than 15
c     residues
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 
      else
         d0=d0_min
      endif

      if(d0.lt.d0_min)d0=d0_min

c     if we passed a particular d0, use it. (useful when comparing a
c     target against a bunch of templates)
      if(m_d0.eq.1)d0=d0_fix

c     for quickly calculate TM-score in searching
      d00=d0                    
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
      d002=d00**2
      nseq=max(nseq1,nseq2)
      do i=1,nseq
         n1(i)=i
         n2(i)=i
      enddo
      
c     do alignment to find invmap(j)
      call super_align          
      
c     resuperpose to find residues of dis<d8
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
c     for recording residue order
            m1(n_al)=i          
            m2(n_al)=j
         endif
      enddo
      d0_input=d0
c     TM-score with dis<d8 only
      call TMscore8(d0_input,n_al,xtm1,ytm1,ztm1,n1,n_al,
     $     xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) 

ccccc Output TM-score is based on the second protein
c     for output
      d0_min=0.5                
c     length for defining final TMscore
      anseq=nseq2               
c     <L>
      if(m_ave.eq.1)anseq=(nseq1+nseq2)/2.0 
      if(anseq.lt.anseq_min)anseq=anseq_min
c     input length
      if(m_fix.eq.1)anseq=L_fix 
      if(anseq.gt.15)then
c     scale for defining TM-score
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      
ccccc remove dis>d8 in normal TM-score calculation for final report
      j=0
      n_eq=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     $        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)
            m1(j)=m1(i)
            m2(j)=m2(i)
            if(ss1(m1(i)).eq.ss2(m2(i)))then
               n_eq=n_eq+1
            endif
         endif
      enddo
      seq_id=float(n_eq)/(n_al+0.00000001)
      n8_al=j
      d0_input=d0
c     normal TMscore
      call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     $     xtm2,ytm2,ztm2,n2,tm8,Rcomm,Lcomm) 
      rmsd8_al=Rcomm
c     TM-score after cutoff
      tm8=tm8*n8_al/anseq     
      tm8sc=tm8  
      tmalign_yang=tm8

ccccc output aligned portions to file alignf, if required
c     count length of filename alignf
      do i=1,512
         if(alignf(i:i).eq.' ')then
            goto 116
         endif
      enddo
 116  i=i-1
c     if non-zero length, then we want the aligned portions
 1236 format('HEADER tm=',F5.3,' coverage=',F4.2,' rmsd=',F5.2)
 1237 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
 1238 format('TER')
 1239 format('CONECT',I5,I5)
c     output template chain
      if(i.gt.0)then
         open(unit=9,file=alignf,status='unknown')
         write(9,1236)tm8,float(n8_al)/nseq2,rmsd8_al
         do i=1,n8_al
            write(9,1237)m1(i),ss1(m1(i)),mm1(m1(i)),xtm1(i),ytm1(i)
     $           ,ztm1(i)
         enddo
c     write TER line
         write(9,1238)
c     connect atoms for proper rendering
         do i=2,n8_al
            write(9,1239)m1(i-1),m1(i) 
         enddo
c     output target chain, begin writing after residue 5000.
         do i=1,n8_al
            write(9,1237)5000+m2(i),ss2(m2(i)),mm2(m2(i)),
     $           xtm2(i),ytm2(i),ztm2(i)
         enddo
         write(9,1238)
         do i=2,n8_al
            write(9,1239)5000+m2(i-1),5000+m2(i)
         enddo
         close(9)
      endif
      return
      end
ccccc tmalign_Yang finished



ccccc structural superposition from a mixture of secondary structure
c     alignment and TM-score
      SUBROUTINE super_align
      PARAMETER(nmax=3000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      dimension gapp(100)

      TMmax=0
      n_gapp=2
      gapp(1)=-0.6
      gapp(2)=0

c      n_gapp=11
c      do i=1,n_gapp
c         gapp(i)=-(n_gapp-i)
c      enddo

ccccc get initial alignment from gapless threading
c     obtain aligment with best estimation of TM-score, invmap_i(i)
      call get_initial          
c     store alignment in invmap(i)
      do i=1,nseq2
         invmap(i)=invmap_i(i)  
      enddo
c     calculate TM with more precision, starting with the previous best
c     estimation. Calculate also matrix score(i,j)
      call get_score
c     refine the aligment if we obtained a better score than the best
c     estimation
      if(TM.gt.TMmax)then
         TMmax=TM
c     store best alignment from threading method in invmap0(j)
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

ccccc iterative alignment, for different gap_open:
c     different gap penalties (actually, n_gapp=2 only)
      do 111 i_gapp=1,n_gapp	
         gap_open=gapp(i_gapp)  
c     iterate between TM-score of a set of aligned residues, producing a
c     new scoring matrix, and dynamic programing producing a new
c     alignment
         do 222 id=1,30         
c     produce alignment invmap(j). we use as score matrix score(i,j). It
c     is as though one residue of the template makes contact with only
c     one residue of the target after structural superposition. The
c     global alignment results in a new invmap(j)
            call dp(nseq1,nseq2) 
c     calculate again TM-score, score(i,j) with new invmap(j)
            call get_score      
c     record the best alignment in whole search
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
c     stop iteration if the difference in TM is too small
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 33
            endif
            TM_old=TM
 222     continue
 33      continue
 111  continue
c     we finished the threading with TM-score

ccccc get initial alignment from secondary structure alignment
      call get_initial2         
      do i=1,nseq2
c     with highest zcore
         invmap(i)=invmap_i(i)  
      enddo
c     find TM score of the aligned residues with the dynamic programming
c     alignment of assigned secondary structure. Calculate also new
c     score(i,j) matrix
      call get_score   
c     update the best aligment if TM-score turns out to be better
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

ccccc iterative alignment, for different gap_open:
c     different gap panalties
      do 1111 i_gapp=1,n_gapp	
c     gap panalty
         gap_open=gapp(i_gapp)  
c     maximum interation is 200
         do 2222 id=1,30	
c     produce alignment invmap(j)
            call dp(nseq1,nseq2) 
c     Input: score(i,j), and gap_open
c     Output: invmap(j)
c     calculate TM-score and score(i,j)
            call get_score      
c     write(*,21)gap_open,rmsd_al,n_al,TM
c     record the best alignment in whole search
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 333
            endif
            TM_old=TM
 2222    continue
 333     continue
 1111 continue
c     we finished the superposition having begun with an alignment of
c     secondary structures
      
c     get initial alignment from invmap0+SS, the output is stored in
c     invmap_i(i)
      call get_initial3         
c     load invmap(i) with previously computed alignment
      do i=1,nseq2
         invmap(i)=invmap_i(i)  
      enddo
c     refien aligment by calculating TM-score, and produce new
c     score(i,j) and invmap(j)
      call get_score            
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

ccccc iterative alignment, for different gap_open:
c     different gap panalties
      do 1110 i_gapp=1,n_gapp	
c     gap panalty
         gap_open=gapp(i_gapp)  
c     maximum interation is 200
         do 2220 id=1,30	
c     produce alignment invmap(j)
            call dp(nseq1,nseq2) 
c     Input: score(i,j), and gap_open
c     Output: invmap(j)

c     calculate TM-score, score(i,j)
            call get_score      
c     record the best alignment in whole search
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 330
            endif
            TM_old=TM
 2220    continue
 330     continue
 1110 continue
c     we finished looking for the best alignment from a combination of
c     secondary structure and rmsd-derived score

      return
      end
ccccc super_align finished


ccccc get initial alignment invmap0(i) from gapless threading. A quick
c     estimation of TM-score is the score for the threading. Store
c     alignment with best scoring in invmap_i(i)
      subroutine get_initial
      parameter(nmax=3000)
      common/backbone/xa(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)

      aL=min(nseq1,nseq2)
c     minimum size of considered fragment
      idel=aL/2.5        
c     5 residues is absolute minimum size
      if(idel.le.5)idel=5

ccccc we fix the target and slide the template from ishift=n1 postiion
c     to ishift=n2 position
c                            0          nseq2
c     -----------------------xxxxxxxxxxxxxx--------------target fixed 
c              -nseq2+idel     idel
c     --------xxxxxxxxxxxxxxxxxxxxxx------------- template at ishif=n1
c                                nseq1-idel
c    ------------------------------xxxxxxxxxxxxxxxxxxxx---    ishif=n2
c     for each ishift, the aligned fragment is the projection of the
c     template onto the target
      n1=-nseq2+idel
      n2=nseq1-idel
      GL_max=0
      do ishift=n1,n2
         L=0
         do j=1,nseq2
            i=j+ishift
            if(i.ge.1.and.i.le.nseq1)then
               L=L+1
               invmap(j)=i
            else
               invmap(j)=-1
            endif
         enddo
         if(L.ge.idel)then
c     estimate TM-score
            call get_GL(GL)
            if(GL.gt.GL_max)then
               GL_max=GL
               do i=1,nseq2
                  invmap_i(i)=invmap(i)
               enddo
            endif
         endif
      enddo

      return
      end
ccccc get_initial finished



ccccc compute initial alignment invmap0(i) from secondary structure
      subroutine get_initial2
      PARAMETER(nmax=3000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)
      
ccccc assign secondary structures
c     1->coil, 2->helix, 3->turn, 4->strand
      do i=1,nseq1
         isec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq1)then
c     dis13 is the square of the distance between position j1 and j3 in
c     the template
            dis13=diszy(0,j1,j3)
            dis14=diszy(0,j1,j4)
            dis15=diszy(0,j1,j5)
            dis24=diszy(0,j2,j4)
            dis25=diszy(0,j2,j5)
            dis35=diszy(0,j3,j5)
c     we need a set of distances to assess the secondary structure at
c     position "i"
            isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
      do i=1,nseq2
         jsec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq2)then
            dis13=diszy(1,j1,j3)
            dis14=diszy(1,j1,j4)
            dis15=diszy(1,j1,j5)
            dis24=diszy(1,j2,j4)
            dis25=diszy(1,j2,j5)
            dis35=diszy(1,j3,j5)
            jsec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
c     smooth the assignment
      call smooth               
      
ccccc build score matrix, from identity matrix:
c        C H T S
c     C  1 0 0 0
c     H  0 1 0 0
c     T  0 0 1 0
c     S  0 0 0 1
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=1
            else
               score(i,j)=0
            endif
         enddo
      enddo

ccccc find initial alignment: invmap(j)
c     should be -1
      gap_open=-1.0             
c     produce alignment invmap(j)
      call dp(nseq1,nseq2) 
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

      return
      end
ccccc get_initial2 finished



ccccc get initial alignment invmap0(i) from secondary structure 
c     and previous alignments
      subroutine get_initial3
      PARAMETER(nmax=3000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)

c     load the best aligment up to now
      do i=1,nseq2
         invmap(i)=invmap0(i)
      enddo

c     get score(i,j) using RMSD matrix
      call get_score1           

      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
c     reinforce score is secondary structure coincides
               score(i,j)=0.5+score(i,j)
            else
               score(i,j)=score(i,j)
            endif
         enddo
      enddo

ccccc find initial alignment: invmap(j)
c     should be -1
      gap_open=-1.0             
c     produce alignment invmap(j) with score(i,j)
      call dp(nseq1,nseq2)      
c     store alignment in invmap_i(i)
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

      return
      end
ccccc get_initial3 finished



ccccc smooth the secondary structure assignment
      subroutine smooth
      PARAMETER(nmax=3000)
      common/sec/isec(nmax),jsec(nmax)
      common/length/nseq1,nseq2
      
c     smooth single
c     --x-- => -----
c     if position "i" is helix or strand, then set as coil if
c     surrounding four positions are not of the same secondary structure
c     type
      do i=1,nseq1
         if(isec(i).eq.2.or.isec(i).eq.4)then
            j=isec(i)
            if(isec(i-2).ne.j)then
               if(isec(i-1).ne.j)then
                  if(isec(i+1).ne.j)then
                     if(isec(i+1).ne.j)then
                        isec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2.or.jsec(i).eq.4)then
            j=jsec(i)
            if(jsec(i-2).ne.j)then
               if(jsec(i-1).ne.j)then
                  if(jsec(i+1).ne.j)then
                     if(jsec(i+1).ne.j)then
                        jsec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo

c   smooth double -------------->
c     --xx-- => ------
c     if position "i+2" and "i+4" is helix or strand, then set as coil
c     if surrounding four positions are not of the same secondary
c     structure type
      do i=1,nseq1
         if(isec(i).ne.2)then
            if(isec(i+1).ne.2)then
               if(isec(i+2).eq.2)then
                  if(isec(i+3).eq.2)then
                     if(isec(i+4).ne.2)then
                        if(isec(i+5).ne.2)then
                           isec(i+2)=1
                           isec(i+3)=1
                        endif
                     endif
                  endif
               endif
            endif
         endif

         if(isec(i).ne.4)then
            if(isec(i+1).ne.4)then
               if(isec(i+2).eq.4)then
                  if(isec(i+3).eq.4)then
                     if(isec(i+4).ne.4)then
                        if(isec(i+5).ne.4)then
                           isec(i+2)=1
                           isec(i+3)=1
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).ne.2)then
            if(jsec(i+1).ne.2)then
               if(jsec(i+2).eq.2)then
                  if(jsec(i+3).eq.2)then
                     if(jsec(i+4).ne.2)then
                        if(jsec(i+5).ne.2)then
                           jsec(i+2)=1
                           jsec(i+3)=1
                        endif
                     endif
                  endif
               endif
            endif
         endif
         
         if(jsec(i).ne.4)then
            if(jsec(i+1).ne.4)then
               if(jsec(i+2).eq.4)then
                  if(jsec(i+3).eq.4)then
                     if(jsec(i+4).ne.4)then
                        if(jsec(i+5).ne.4)then
                           jsec(i+2)=1
                           jsec(i+3)=1
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo
      
c     connect -------------->
c     x-x => xxx
      do i=1,nseq1
         if(isec(i).eq.2)then
            if(isec(i+1).ne.2)then
               if(isec(i+2).eq.2)then
                  isec(i+1)=2
               endif
            endif
         endif
         
         if(isec(i).eq.4)then
            if(isec(i+1).ne.4)then
               if(isec(i+2).eq.4)then
                  isec(i+1)=4
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2)then
            if(jsec(i+1).ne.2)then
               if(jsec(i+2).eq.2)then
                  jsec(i+1)=2
               endif
            endif
         endif
         
         if(jsec(i).eq.4)then
            if(jsec(i+1).ne.4)then
               if(jsec(i+2).eq.4)then
                  jsec(i+1)=4
               endif
            endif
         endif
      enddo
      
      return
      end
ccccc smooth finished



ccccc square of distance between position i1 and position i2 in sequence
c     "i"
      function diszy(i,i1,i2)
      parameter(nmax=3000)
      common/backbone/xa(3,nmax,0:1)
      diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     $     +(xa(2,i1,i)-xa(2,i2,i))**2
     $     +(xa(3,i1,i)-xa(3,i2,i))**2)
      return
      end
ccccc diszy finished


ccccc assign secondary structure based on a set of distance between
c     surrounding positions along the sequence
c     1->coil, 2->helix, 3->turn, 4->strand
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
c     initialize to coil
      make_sec=1
      delta=2.1

      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
c     it's a helix
                        make_sec=2 
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif

      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
c     strand
                        make_sec=4 
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif

c     turn
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end
ccccc make_sec finished



ccccc quickly calculate TM-score with given invmap(i) in 3 iterations
c     and do not consider subfragments of the aligned residues, but only
c     the whole set of aligned residues
      subroutine get_GL(GL)
      PARAMETER(nmax=3000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d00/d00,d002

      dimension xo1(nmax),yo1(nmax),zo1(nmax)
      dimension xo2(nmax),yo2(nmax),zo2(nmax)
      dimension dis2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms
      data w /nmax*1.0/
ccc   

ccccc calculate RMSD between aligned structures and rotate the
c     structures

c     n_al:number of aligned residues
      n_al=0
      do j=1,nseq2
c     j in target aligned to i in template
         i=invmap(j)            
         if(i.gt.0)then
            n_al=n_al+1
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
            xo1(n_al)=xa(1,i,0)
            yo1(n_al)=xa(2,i,0)
            zo1(n_al)=xa(3,i,0)
            xo2(n_al)=xa(1,j,1)
            yo2(n_al)=xa(2,j,1)
            zo2(n_al)=xa(3,j,1)
         endif
      enddo
c     u rotate r_1 to r_2
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) 
      GL=0
c     transform template coordinates to superimpose to target
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         GL=GL+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration, keep only those aligned residues whose square
ccc   distance below d002t=d0**2
      d002t=d002
 21   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
c     if d002t was very restrictive and few aligned residues were kept,
c     then increase d002t
      if(j.lt.3)then
         d002t=d002t+.5
         goto 21
      endif
      L=j
c     u rotate r_1 to r_2
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) 
      G2=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G2=G2+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration, increase d002t a little bit and repeat
      d002t=d002+1
 22   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3)then
         d002t=d002t+.5
         goto 22
      endif
      L=j
c     u rotate r_1 to r_2
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) 
      G3=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G3=G3+1/(1+dis2(i)/(d0**2))
      enddo

c     select TM-score as biggest of GL,G2,G3
      if(G2.gt.GL)GL=G2
      if(G3.gt.GL)GL=G3
      return
      end
ccccc get_GL finished


ccccc with invmap(i) calculate TM-score and matrix score(i,j) for
cccc  rotation
      subroutine get_score
      PARAMETER(nmax=3000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
c     armsd is real
      double precision u(3,3),t(3),rms,drms 
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the
c     structures
      n_al=0
      do j=1,NSEQ2
c     j position in template aligned to i position in target
         i=invmap(j)            
         if(i.gt.0)then
            n_al=n_al+1
ccc   for TM-score:
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
         endif
      enddo
***   calculate TM-score for the given alignment
      d0_input=d0
c     simplified search engine
      call TMscore8_search(d0_input,n_al,xtm1,ytm1,ztm1,n1,
     $     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) 
c     TM-score
      TM=TM*n_al/anseq          
c     calculate score matrix score(i,j)
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo

ccccc compute score(i,j), which related to the distance between position
c     "i" of template with position "j" of target after best
c     superposition of template to target u rotate r_1 to r_2
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) 
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d0**2)
         enddo
      enddo
      
      return
      end
ccccc get_score finished



ccccc with invmap(i) calculate score(i,j) using RMSD rotation 
      subroutine get_score1
      PARAMETER(nmax=3000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      common/d0min/d0_min
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the
c     structures
      n_al=0
      do j=1,NSEQ2
c     j aligned to i
         i=invmap(j)            
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
c     calculate score matrix score(i,j)
c     u rotate r_1 to r_2
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) 
c     d01 cutoff for contact between one position in the template, and
c     one position in the target
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d02)
         enddo
      enddo

      return
      end
ccccc get_score1 finished



ccccc This is a subroutine to compare two structures and find the 
c     superposition that has the maximum TM-score.
c
c     L1--Length of the first structure
c     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
c     n1(i)--Residue sequence number of i'th residue at the first structure
c     L2--Length of the second structure
c     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
c     n2(i)--Residue sequence number of i'th residue at the second structure
c     TM--TM-score of the comparison
c     Rcomm--RMSD of two structures in the common aligned residues
c     Lcomm--Length of the common aligned regions
c
c     Note: 
c     1, Always put native as the second structure, by which TM-score
c        is normalized.
c     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
c        TM-score superposition.
c dis<8, simplified search engine
      subroutine TMscore8_search(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2, TM
     $     ,Rcomm,Lcomm)
      PARAMETER(nmax=3000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
c     [1,n_ali],align residues for the score
      common/nscore/i_ali(nmax),n_cut 
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)
      dimension iL0(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms
      data w /nmax*1.0/
ccc   

ccccc convert input data
c     because L1=L2 in this special case
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
c     number of aligned residues
      n_ali=L1                  
      Lcomm=L1

c     parameters:
c     d0
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
c     d0_search
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
c     iterative parameters
c     maximum number of iterations
      n_it=20                   
c     for output alignment
      d_output=5                
c     maximum number of L_init
      n_init_max=6              
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
c     the different fragment sizes to compare
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

c     find the maximum score starting from local structures
c     superposition
      score_max=-1              
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        k=0
c     this is the simplification. Instead of sliding the fragment along
c     the aligned residues one by one, we do jumps of 40 residues
        do i=1,iL_max,40        
          k=k+1
          iL0(k)=i
        enddo
        if(iL0(k).lt.iL_max)then
          k=k+1
          iL0(k)=iL_max
        endif
        n_shift=k
        do 300 i_shift=1,n_shift
           iL=iL0(i_shift)
           LL=0
           ka=0
           do i=1,L_init
c     [1,n_ali] common aligned
              k=iL+i-1          
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
c     u rotate r_1 to r_2
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
c     global superposition
           if(i_init.eq.1)then  
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
c     init, get scores, n_cut+i_ali(i) for iteration
           call score_fun8       
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
ccccc iteration for extending
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
c     [1,n_ali]
                 m=i_ali(i)     
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
c     u rotate r_1 to r_2
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
c     get scores, n_cut+i_ali(i) for iteration
              call score_fun8    
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
c     for iteration
 301       continue             
 302       continue
c     for shift
 300    continue                
c     for initial length, L_ali/M
 333  continue                  

ccccc return the final rotation
      LL=0
      do i=1,ka0
c     record of the best alignment
         m=k_ali0(i)            
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
c     u rotate r_1 to r_2
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

      return
      END
ccccc TMscore8_search finished



c     This is a subroutine to compare two structures and find the 
c     superposition that has the maximum TM-score.
c
c     L1--Length of the first structure
c     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
c     n1(i)--Residue sequence number of i'th residue at the first structure
c     L2--Length of the second structure
c     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
c     n2(i)--Residue sequence number of i'th residue at the second structure
c     TM--TM-score of the comparison
c     Rcomm--RMSD of two structures in the common aligned residues
c     Lcomm--Length of the common aligned regions
c
c     Note: 
c     1, Always put native as the second structure, by which TM-score
c        is normalized.
c     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
c        TM-score superposition.
c     dis<8, but same search engine
      subroutine TMscore8(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     $     TM,Rcomm,Lcomm)
      PARAMETER(nmax=3000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
c     [1,n_ali],align residues for the score
      common/nscore/i_ali(nmax),n_cut 
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
c     armsd is real
      double precision u(3,3),t(3),rms,drms 
      data w /nmax*1.0/
ccc   
      
ccccc convert input data
c     because L1=L2 in this special case
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
c     number of aligned residues
      n_ali=L1                  
      Lcomm=L1

c     parameters:
c     d0
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
c     d0_search
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
c     iterative parameters
c     maximum number of iterations
      n_it=20                   
c     for output alignment
      d_output=5          
c      maximum number of L_init
      n_init_max=6              
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

c     find the maximum score starting from local structures
c     superposition
c     TM-score
      score_max=-1              
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
c     on aligned residues, [1,nseqA]
        do 300 iL=1,iL_max    
           LL=0
           ka=0
           do i=1,L_init
c     [1,n_ali] common aligned
              k=iL+i-1          
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
c     u rotate r_1 to r_2
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
c     global superposition
           if(i_init.eq.1)then  
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
c     init, get scores, n_cut+i_ali(i) for iteration
           call score_fun8       
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
c     iteration for extending
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
c     [1,n_ali]
                 m=i_ali(i)     
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
c     u rotate r_1 to r_2
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
c     get scores, n_cut+i_ali(i) for iteration
              call score_fun8    
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
c     for iteration
 301       continue             
 302       continue
c     for shift
 300    continue                
c     for initial length, L_ali/M
 333  continue                  
      
ccccc return the final rotation
      LL=0
      do i=1,ka0
c     record of the best alignment
         m=k_ali0(i)            
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
c     u rotate r_1 to r_2
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

      return
      END

ccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccc
      subroutine score_fun8
      PARAMETER(nmax=3000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
c     [1,n_ali],align residues for the score
      common/nscore/i_ali(nmax),n_cut 
      common/scores/score
      double precision score,score_max
      common/d8/d8

c     number of residue-pairs dis<d, for iteration
      n_cut=0                   
c     TMscore
      score_sum=0               
      do k=1,n_ali
c     [1,nseqA] reoder number of structureA
         i=iA(k)                
c     [1,nseqB]
         j=iB(k)                
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d)then
            n_cut=n_cut+1
c     [1,n_ali], mark the residue-pairs in dis<d
            i_ali(n_cut)=k      
         endif
         if(dis.le.d8)then
            score_sum=score_sum+1/(1+(dis/d0)**2)
         endif
      enddo
c     TM-score
      score=score_sum/float(nseqB) 

      return
      end

c     This is a subroutine to compare two structures and find the 
c     superposition that has the maximum TM-score.
c
c     L1--Length of the first structure
c     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
c     n1(i)--Residue sequence number of i'th residue at the first structure
c     L2--Length of the second structure
c     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
c     n2(i)--Residue sequence number of i'th residue at the second structure
c     TM--TM-score of the comparison
c     Rcomm--RMSD of two structures in the common aligned residues
c     Lcomm--Length of the common aligned regions
c
c     Note: 
c     1, Always put native as the second structure, by which TM-score
c        is normalized.
c     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
c        TM-score superposition.
c     normal TM-score:
      subroutine TMscore(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     $     TM,Rcomm,Lcomm)
      PARAMETER(nmax=3000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
c     [1,n_ali],align residues for the score
      common/nscore/i_ali(nmax),n_cut 
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
c     armsd is real
      double precision u(3,3),t(3),rms,drms 
      data w /nmax*1.0/
ccc   
      
ccccc convert input data
c     because L1=L2 in this special case
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
c     number of aligned residues
      n_ali=L1                  
      Lcomm=L1

c     parameters:
c     d0
c     d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
c     d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
c     iterative parameters ----->
c     maximum number of iterations
      n_it=20                   
c     for output alignment
      d_output=5          
c     maximum number of L_init
      n_init_max=6              
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

c     find the maximum score starting from local structures
c     superposition
c     TM-score
      score_max=-1              
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
c     on aligned residues, [1,nseqA]
        do 300 iL=1,iL_max      
           LL=0
           ka=0
           do i=1,L_init
c     [1,n_ali] common aligned
              k=iL+i-1          
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
c     u rotate r_1 to r_2
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
c     global superposition
           if(i_init.eq.1)then  
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
c     init, get scores, n_cut+i_ali(i) for iteration
           call score_fun       
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
c     iteration for extending
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
c     [1,n_ali]
                 m=i_ali(i)     
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
c     u rotate r_1 to r_2
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
c     get scores, n_cut+i_ali(i) for iteration
              call score_fun    
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
c     for iteration
 301       continue             
 302       continue
c     for shift
 300    continue                
c     for initial length, L_ali/M
 333  continue                  

ccccc return the final rotation
      LL=0
      do i=1,ka0
c     record of the best alignment
         m=k_ali0(i)            
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
c     u rotate r_1 to r_2
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) 
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

      return
      END

ccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccc
      subroutine score_fun
      PARAMETER(nmax=3000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
c     [1,n_ali],align residues for the score
      common/nscore/i_ali(nmax),n_cut 
      common/scores/score
      double precision score
c     number of residue-pairs dis<d, for iteration
      n_cut=0                   
c     TMscore
      score_sum=0               
      do k=1,n_ali
c     [1,nseqA] reoder number of structureA
         i=iA(k)                
c     [1,nseqB]
         j=iB(k)                
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d)then
            n_cut=n_cut+1
c     [1,n_ali], mark the residue-pairs in dis<d
            i_ali(n_cut)=k      
         endif
         score_sum=score_sum+1/(1+(dis/d0)**2)
      enddo
c     TM-score
      score=score_sum/float(nseqB) 

      return
      end



ccccc Dynamic programming for alignment.
c     Input: score(i,j), and gap_open
c     output: invmap(j)
      subroutine dp(nseq1,nseq2)
      parameter(nmax=3000)
      logical*1 dir
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      dimension dir(0:nmax,0:nmax),val(0:nmax,0:nmax)
      real h,v
      
c     initialize the matrix:
      val(0,0)=0
      do i=1,nseq1
        dir(i,0)=.false.
        val(i,0)=0
      enddo
      do j=1,nseq2
        dir(0,j)=.false.
        val(0,j)=0
        invmap(j)=-1
      enddo

c     decide matrix and path:
      do j=1,nseq2
         do i=1,nseq1
            d=val(i-1,j-1)+score(i,j)
            h=val(i-1,j)
            if(dir(i-1,j))h=h+gap_open
            v=val(i,j-1)
            if(dir(i,j-1))v=v+gap_open
            
          if((d.ge.h).and.(d.ge.v)) then
             dir(i,j)=.true.
             val(i,j)=d
          else
             dir(i,j)=.false.
             if(v.ge.h)then
                val(i,j)=v
             else
                val(i,j)=h
             end if
          endif
       enddo
      enddo
      
c     extract the alignment:
      i=nseq1
      j=nseq2
      do while((i.gt.0).and.(j.gt.0))
         if(dir(i,j))then
            invmap(j)=i
            i=i-1
            j=j-1
         else
            h=val(i-1,j)
            if(dir(i-1,j))h=h+gap_open
            v=val(i,j-1)
            if(dir(i,j-1))v=v+gap_open
            if(v.ge.h) then
               j=j-1
            else
               i=i-1
            endif
         endif
      enddo
      
      return
      end
ccccc dp finished



ccccc Calculate sum of (r_d-r_m)^2
c     w - w(m) is weight for atom pair c m (given) x - x(i,m) are
c     coordinates of atom c m in set x (given) y - y(i,m) are
c     coordinates of atom c m in set y (given) n - n is number of atom
c     pairs (given) mode - 0:calculate rms only (given) 1:calculate
c     rms,u,t (takes longer) rms - sum of w*(ux+t-y)**2 over all atom
c     pairs (result) u - u(i,j) is rotation matrix for best
c     superposition (result) t - t(i) is translation vector for best
c     superposition (result) ier - 0: a unique optimal superposition has
c     been determined(result) -1: superposition is not unique but
c     optimal -2: no result obtained because of negative weights w or
c     all weights equal to zero.
ccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(1), x(3, 1), y(3, 1), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     $e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     $, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     $ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     $, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     $ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     $, e2), (e(3), e3)
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
c     DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
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
c     DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
         e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i))
     $        **2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     $     ) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3)
     $     * ((r(2,1) * r(3,2)) - (r(2,2) * r(3,1))))
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
     $)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     $) + (rr1 * rr3)) - (rr2 * rr2)) / three
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
c****************** ROTATION MATRIX
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     $))
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
     $))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     $ * xc(3))
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

