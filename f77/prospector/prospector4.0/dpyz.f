*     idir(i,j)=1,2,3, from diagonal, horizontal, vertical
*     val(i,j) is the cumulative score of (i,j)
********************************************************************
      subroutine DP3(score0)
      PARAMETER(ndim=3000)
      common/dpc/score(ndim,ndim),gap_open,gap_extn,j2i(ndim),nseq1,nseq2
      
      dimension val(0:ndim,0:ndim),idir(0:ndim,0:ndim)
      dimension preV(0:ndim,0:ndim),preH(0:ndim,0:ndim),preD(0:ndim,0:ndim)
      dimension idirH(0:ndim,0:ndim),idirV(0:ndim,0:ndim)
      
      val(0,0)=0.0
      do i=1,nseq1
         val(i,0)=0
         idir(i,0)=0
         preD(i,0)=0.0
         preH(i,0)=-1000.0
         preV(i,0)=-1000.0
      enddo
      do j=1,nseq2
         val(0,j)=0
         idir(0,j)=0
         preD(0,j)=0.0
         preH(0,j)=-1000.0
         preV(0,j)=-1000.0
      enddo
      
      do 111 i=1,nseq1
         do 222 j=1,nseq2
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
 222     continue
 111  continue
      score0=val(nseq1,nseq2)
      
      do j=1,nseq2
         j2i(j)=-1
      enddo
      i=nseq1
      j=nseq2
      do while(i.gt.0.and.j.gt.0)
         if(idir(i,j).eq.1)then
            j2i(j)=i
            i=i-1
            j=j-1
         elseif(idir(i,j).eq.2)then
            i=i-1
            idir(i,j)=idirH(i,j)
         else
            j=j-1
            idir(i,j)=idirV(i,j)
         endif
      enddo
      
      return
      end
      
