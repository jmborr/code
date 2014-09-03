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
	common/path/ix(0:maxres,0:nmax),iy(0:maxres,0:nmax)
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
	iix=ix(i,j)
	iiy=iy(i,j)
		if(val(i+iix,j+iiy).gt.0)then
                i=i+iix
		j=j+iiy
		else
		go to 10
		end if
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
	common/path/ix(0:maxres,0:nmax),iy(0:maxres,0:nmax)
	common/cont2/ior(15,0:nmax,0:ncases)
	common/gv/go,gv
       REAL H, V
C Init matrices

        do j=1,nseq2
         do i=1,nseq1	 
	    dir(i,j)=.false.
	    val(i,j)=0.
	    ix(i,j)=0
	    iy(i,j)=0
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
		ix(i,j)=-1
		iy(i,j)=0
               ELSE
                 VAL(i,j)=V
                 DIR(i,j)=.true.
		ix(i,j)=0
		iy(i,j)=-1
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
