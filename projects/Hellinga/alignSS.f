c     **************************************************************
c
c     pgf77 -Mextend -O -s -fast -Wl,-static -o alignSS.x alignSS.f
c
c     This program is to do alignment between two fasta-like files with
c     secondary structure assignment (characters are H,E,-)
c
c     **************************************************************
      program alignSS
      PARAMETER(ndim=3000)
c     number of secondary structure types
      parameter(nss=4) 
      common/dpc/score(ndim,ndim),gap_open,gap_extn,j2i(ndim)
     $     ,nseq1,nseq2
c     b,z,x are additional
      common/matra/imut(nss,nss)     

      integer seq1(ndim),seq2(ndim)
      character*10000 fnam1,fnam2
      character*10000 s
      character ss(nss)
      character*100 du,ad
      character sequenceA(ndim),sequenceB(ndim),sequenceM(ndim)

      data ss/'C','H','S','E'/
      call getarg(1,fnam1)
      call getarg(2,fnam2)

 1    format(A10000)
c     read sequence1:
      open(unit=10,file=fnam1,status='old')
c     read past the header
      do while(.true.)
         read(10,1,end=33)s
         if(s(1:1).eq.'>')then
            goto 329
         endif
      enddo
 329  continue
c     read the sequence
      i=0
      do while(.true.)
         read(10,1,end=33)s
         do k=1,10000
            do j=1,nss
               if(s(k:k).eq.ss(j)) then
                  i=i+1
                  seq1(i)=j
                  goto 331
               endif
            enddo
 331        if(i.ge.ndim) goto 33
         enddo
      enddo
 33   continue
      close(10)
      nseq1=i
c     read sequence2:
      open(unit=10,file=fnam2,status='old')
      i=0
      do while(.true.)
         read(10,1,end=44)s
         do k=1,10000
            do j=1,nss
               if(s(k:k).eq.ss(j))then
                  i=i+1
                  seq2(i)=j
                  goto 441
               endif
            enddo
 441        if(i.ge.ndim)goto 44
         enddo
      enddo
 44   continue
      close(10)
      nseq2=i

ccccc read mutation  matrix
      call matrix               

ccccc score
      do i=1,nseq1
         do j=1,nseq2
            score(i,j)=imut(seq1(i),seq2(j))
         enddo
      enddo

ccccc dynamic program:
      gap_open=-1
      gap_extn=-0.333
      call DP2(score0)

ccccc calculate sequence identity
      L_id=0
      L_ali=0
      do j=1,nseq2
         if(j2i(j).gt.0)then
            i=j2i(j)
            L_ali=L_ali+1
            if(seq1(i).eq.seq2(j))L_id=L_id+1
         endif
      enddo

      write(*,*)
      write(*,101)nseq1,fnam1
 101  format('Length of sequence 1: ',I4,' ->',A10)
      write(*,102)nseq2,fnam2
 102  format('Length of sequence 2: ',I4,' ->',A10)
      write(*,103)L_ali
 103  format('Aligned length: ',I4)
      write(*,104)L_id
 104  format('Identical length: ',I4)
      write(*,105)float(L_id)/(L_ali+0.00000001),L_id,L_ali
 105  format('Sequence identity: ',F8.3,' (=',I4,'/',I4,')')
      write(*,*)

ccccc output aligned sequences
c     final aligned order
      k=0                       
c     on sequence 1
      i=1                       
c     on sequence 2
      j=1                       
 800  continue
      if(i.gt.nseq1.and.j.gt.nseq2)goto 802
c     unaligned C on 1
      if(i.gt.nseq1.and.j.le.nseq2)then 
         k=k+1
         sequenceA(k)='-'
         sequenceB(k)=ss(seq2(j))
         sequenceM(k)=' '
         j=j+1
         goto 800
      endif
c     unaligned C on 2
      if(i.le.nseq1.and.j.gt.nseq2)then 
         k=k+1
         sequenceA(k)=ss(seq1(i))
         sequenceB(k)='-'
         sequenceM(k)=' '
         i=i+1
         goto 800
      endif
c     if aligned
      if(i.eq.j2i(j))then    
         k=k+1
         sequenceA(k)=ss(seq1(i))
         sequenceB(k)=ss(seq2(j))
c     identical
         if(seq1(i).eq.seq2(j))then 
            sequenceM(k)=':'
         else
            sequenceM(k)=' '
         endif
         i=i+1
         j=j+1
         goto 800
c     if gap on 1
      elseif(j2i(j).lt.0)then   
         k=k+1
         sequenceA(k)='-'
         sequenceB(k)=ss(seq2(j))
         sequenceM(k)=' '
         j=j+1
         goto 800
c     if gap on 2
      elseif(j2i(j).gt.0)then 
         k=k+1
         sequenceA(k)=ss(seq1(i))
         sequenceB(k)='-'
         sequenceM(k)=' '
         i=i+1
         goto 800
      endif
 802  continue

      write(*,601)(sequenceA(i),i=1,k)
      write(*,601)(sequenceM(i),i=1,k)
      write(*,601)(sequenceB(i),i=1,k)
      write(*,602)(mod(i,10),i=1,k)
 601  format(2000A1)
 602  format(2000I1)
      write(*,*)

 999  END


ccccc
c     This is a standard Needleman-Wunsch dynamic program (by Y. Zhang 2005)
c     1. Count multiple-gap.
c     2. The gap penality W(k)=Go+Ge*k1+Go+Ge*k2 if gap open on both sequences
c
c     Input: score(i,j), gap_open, gap_extn
c     Output: j2i(j)
c     idir(i,j)=1,2,3, from diagonal, horizontal, vertical
c     val(i,j) is the cumulative score of (i,j)
ccccc 
      subroutine DP2(score0)
      PARAMETER(ndim=3000)
      common/dpc/score(ndim,ndim),gap_open,gap_extn,j2i(ndim)
     $     ,nseq1,nseq2
      dimension val(0:ndim,0:ndim),idir(0:ndim,0:ndim)
      dimension jpV(0:ndim,0:ndim),jpH(0:ndim,0:ndim)
      dimension preV(0:ndim,0:ndim),preH(0:ndim,0:ndim)
      real D,V,H

ccccc initializations
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

ccccc DP
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
ccccc DP finished

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccc read matrix
      subroutine matrix
c     number of amino acid
      parameter(nss=4) 
c     b,z,x are additional
      common/matra/imut(nss,nss) 
c     'C'->1 'H'->2, 'S'->3 , 'E'->4

      imut(1,1)=0.333
      imut(2,2)=1
      imut(3,3)=1
      imut(4,4)=1

      imut(1,2)=0
      imut(2,2)=0

      imut(1,3)=0
      imut(3,1)=0

      imut(1,4)=0
      imut(4,1)=0

      imut(2,3)=0
      imut(3,2)=0

      imut(2,4)=-1 !incompatible helix/strand
      imut(4,2)=-1

      return
      end
c     end of subroutine matrix

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccc smooth the secondary structure assignment
      subroutine smooth(seq1,seq2)
      PARAMETER(ndim=3000)
      integer seq1(ndim),seq2(ndim)
      common/dpc/score(ndim,ndim),gap_open,gap_extn,j2i(ndim)
     $     ,nseq1,nseq2
      
c     smooth single
c     --x-- => -----
c     if position "i" is helix or strand, then set as coil if
c     surrounding four positions are not of the same secondary structure
c     type
      do i=1,nseq1
         if(seq1(i).eq.2.or.seq1(i).eq.4)then
            j=seq1(i)
            if(seq1(i-2).ne.j)then
               if(seq1(i-1).ne.j)then
                  if(seq1(i+1).ne.j)then
                     if(seq1(i+1).ne.j)then
                        seq1(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(seq2(i).eq.2.or.seq2(i).eq.4)then
            j=seq2(i)
            if(seq2(i-2).ne.j)then
               if(seq2(i-1).ne.j)then
                  if(seq2(i+1).ne.j)then
                     if(seq2(i+1).ne.j)then
                        seq2(i)=1
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
         if(seq1(i).ne.2)then
            if(seq1(i+1).ne.2)then
               if(seq1(i+2).eq.2)then
                  if(seq1(i+3).eq.2)then
                     if(seq1(i+4).ne.2)then
                        if(seq1(i+5).ne.2)then
                           seq1(i+2)=1
                           seq1(i+3)=1
                        endif
                     endif
                  endif
               endif
            endif
         endif

         if(seq1(i).ne.4)then
            if(seq1(i+1).ne.4)then
               if(seq1(i+2).eq.4)then
                  if(seq1(i+3).eq.4)then
                     if(seq1(i+4).ne.4)then
                        if(seq1(i+5).ne.4)then
                           seq1(i+2)=1
                           seq1(i+3)=1
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(seq2(i).ne.2)then
            if(seq2(i+1).ne.2)then
               if(seq2(i+2).eq.2)then
                  if(seq2(i+3).eq.2)then
                     if(seq2(i+4).ne.2)then
                        if(seq2(i+5).ne.2)then
                           seq2(i+2)=1
                           seq2(i+3)=1
                        endif
                     endif
                  endif
               endif
            endif
         endif
         
         if(seq2(i).ne.4)then
            if(seq2(i+1).ne.4)then
               if(seq2(i+2).eq.4)then
                  if(seq2(i+3).eq.4)then
                     if(seq2(i+4).ne.4)then
                        if(seq2(i+5).ne.4)then
                           seq2(i+2)=1
                           seq2(i+3)=1
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
         if(seq1(i).eq.2)then
            if(seq1(i+1).ne.2)then
               if(seq1(i+2).eq.2)then
                  seq1(i+1)=2
               endif
            endif
         endif
         
         if(seq1(i).eq.4)then
            if(seq1(i+1).ne.4)then
               if(seq1(i+2).eq.4)then
                  seq1(i+1)=4
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(seq2(i).eq.2)then
            if(seq2(i+1).ne.2)then
               if(seq2(i+2).eq.2)then
                  seq2(i+1)=2
               endif
            endif
         endif
         
         if(seq2(i).eq.4)then
            if(seq2(i+1).ne.4)then
               if(seq2(i+2).eq.4)then
                  seq2(i+1)=4
               endif
            endif
         endif
      enddo
      
      return
      end
ccccc smooth finished

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
