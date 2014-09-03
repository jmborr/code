	SUBROUTINE DP(NSEQ1,NSEQ2,RESULT,ils,ilen)
	PARAMETER(nmax=1000)
	PARAMETER(ncases=5700)
        parameter (maxres=500)       
	LOGICAL*1 DIR
       COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *    GAP_OPEN, GAP_EXTN
	common/dir2/DIR(0:MAXRES,0:NMAX)         
	common/mappings/jcode(-1:maxres)
	common/mappings2/jk,lst
	common/alignment/invmap(0:nmax,0:ncases,4)		
        dimension idir(0:maxres,0:nmax)
        dimension preV(0:maxres,0:nmax),preH(0:maxres,0:nmax),preD(0:maxres,0:nmax)
        dimension idirH(0:maxres,0:nmax),idirV(0:maxres,0:nmax)

c	DO DYNAMIC PROGRAMMING FIRST
c	===========================================================
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

		IF(val(i,j).lt.0)then
		val(i,j)=0.
		end if

 222     continue
 111  continue


c	===========================================================

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
		if(idir(i,j).eq.1)then
  		if(val(i,j).gt.0)then
		invmap(j,jk,lst)=i
        	     i=i-1
	             j=j-1
                 else
			go to 10	
		end if
           ELSEif(idir(i,j).eq.2)then
		if(val(i-1,j).gt.0)then
                i=i-1
	       idir(i,j)=idirH(i,j)
		else
		go to 10
		end if
	else 	
		if(val(i,j-1).gt.0)then
	        j=j-1
               idir(i,j)=idirV(i,j)
		else
		go to 10
		end if
             ENDIF
         ENDDO
10	result=best
	return
       END

c	update the score matrix for residues above
c	residue j is updated
	return
	end 

c
