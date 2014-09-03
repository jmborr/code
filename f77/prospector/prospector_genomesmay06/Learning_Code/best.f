c	based on c10b6c7b7.f
c	4/7/01  Jeffrey Skolnick
	parameter(nmax=1600)
	parameter(ncases=4200)
	character*5 name(0:ncases),name2
	character*3 nametarg,name3
	logical*1 dir
	dimension icode(nmax,0:ncases),nres(0:ncases)
	common/nam/name3

	open(unit=9,file='stat_best')
	call r14stat(name,nres,icode)
	open(unit=1,file='finish_best')
	rewind(1)
	write(1,*)'1'
	close(1)
	stop
	end


	
	SUBROUTINE r14STAT(name,NRES,ICODE)
	PARAMETER(NMAX=1600)
	PARAMETER(NCASES=4200)
	PARAMETER(MAXRES=1800)
	PARAMETER(REJECT=25.)	
	CHARACTER*5 NAME(0:NCASES),NAMESELECT2(10000),NAMESELECT
	CHARACTER*5 NAMETARG
	CHARACTER*3 NAME3,TER
	character*20 a33
	character*1 a3,qdir
	CHARACTER*6 a6
	CHARACTER*60 head
	CHARACTER*255 NAMEP,namev(ncases)
	CHARACTER*1 AA(0:19),ad
	CHARACTER*21 TEXT
	CHARACTER*5 name5
	character*4 genome
	CHARACTER*6 prodin
	LOGICAL*1 DIR,kselect,hommap(NCASES)
	double precision u,t,umin(3,3),tm(3)
	common/matrix/u(3,3),t(3)
c	DIMENSION AMUT(-1:19,0:NMAX)	
	DIMENSION XSC(3,NMAX),XSCT(3,NMAX)
	DIMENSION XSC2(3,NMAX),XSCT2(3,NMAX)	
	COMMON/C/NDATA,NDATA2
	DIMENSION XA(3,NMAX,0:NCASES)
        COMMON /ALIGNPARAM/ SCORE(0:MAXRES,0:NMAX), 
     *   VAL(0:MAXRES,0:NMAX),
     *   LIGNMENT(0:2*MAXRES), GAP_OPEN, GAP_EXTN
	COMMON/DIR2/DIR(0:NMAX,0:NMAX)    
	COMMON/SEQALL/IMUT(-1:19,-1:19)	
	COMMON/ALIGNMENT/INVMAP(0:NMAX,0:NCASES,4)
	COMMON/MAPPINGS2/JK,LST
	COMMON/nam/name3
	COMMON/MAPPINGS/JCODE(-1:MAXRES),IMAP(0:MAXRES)     	
	dimension isz(3000)	
	DIMENSION ICODE(NMAX,0:NCASES)
	DIMENSION NRES(0:NCASES)
	DIMENSION IMAPOLD(-5:NMAX,NCASES)
	DIMENSION AMUT(-1:19,0:NMAX)
	dimension ist(0:2),ist2(0:2)
	common/rotate/mres,xs(3,NMAX)
        dimension energ(0:ncases),nm(0:ncases),item(0:ncases)
	dimension zminrms(0:ncases),sa(3),sb(3)
	dimension zminrmso(0:ncases)
	dimension ires1(0:nmax),jres1(0:nmax),rt(0:ncases)
        parameter  (maxseq =1600)	
        common/seq/nresseq,nseq_ori,kcode(maxres,maxseq)	
	
	DATA AA/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/
	open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)

	open(unit=55,file='nametarg1')
        read(55,13)prodin
13      format(a6)
	OPEN(unit=4,file='LIST.targ'//prodin)
	read(4,*)nd,genome
	open(unit=63,file='listbest')
	OPEN(unit=33, file='STRUCTURE_rap3orienrev')
	read(33,*)ndata2
	itot=0
	do 1001 IIk=1,NDATA2
        read(33,9432)iselected,namep(1:nsize),nameselect,zselect,head
9432    format(i6,3x,a,3x,a5,1f7.3,,3x,a60)
        OPEN(unit=3,file=namep(1:nsize)//'.threadrap3orienrev')
	read(3,*)nth
	if(nth.eq.0)go to 1001
	read(3,25)name5,zth
	if(zth .gt. 25)go to 22
	if(nth .gt.1) go to 1001
25      format(a5,1x,1f7.3)
        if(zth .lt.11 )then
        go to 1001
	end if
22	write(9,*)namep(1:nsize),zth,name5
        itot=itot+1
        namev(itot)(1:nsize)=namep(1:nsize)

1001	continue
        write(63,*)itot
        do ik=1,itot
        write(63,63)namev(ik)(1:nsize)
63      format(a)
        end do

        write(9,*)'total number of proteins examined',itot
        write(9,*)'average fraction of coverage=',atc/itot


	RETURN
	END 	

	subroutine rmscalc(lenf2,xsc,xsct,rrms,drms)	
	PARAMETER(NMAX=1600)
	PARAMETER(NCASES=4200)	
	dimension xsc(3,nmax),xsct(3,nmax),x3(3)
c this program will compare (on the rms level) the PDB file and
c output from lattice simulations 
c

c            version 1.0 (19 Novemebr 1990)
c
c   written in Novebmer 1990, based of programs FINDSITE
c   and various generations of earlier programs (SEARCH, KONTACT etc)
c
c  REMBER TO UPGRADE VERSION NUMBER AFTER INTRODUCING CHANGES !!!
c      1.0 - just started
c---------------------------------------------------
c
c                      INPUT SECTION
c
c---------------------------------------------------
c      program rms
	parameter(ndim=1600)
c
c all variables with names ending in q are for query protein
c                                    t are for target protein
c maximum number of atoms and residues are defined as parameters
	real rrms, drms
     	double precision w(ndim),u, t, rms, sigma
c	double precision xxta(3,ndim)
	double precision xxq(3,ndim),xxt(3,ndim)
c	double precision xxt1(3,ndim),xxq1(3,0:ndim)
c	double precision xxt2(3,ndim),xx1(3,ndim)
	common/rotate/mres,xs(3,NMAX)
	common/matrix/u(3,3),t(3)
    	 version = 1.0

      do i = 1, 1000
        w(i) = 1.0
      end do

	ind=0
	do i=1,LENF2
		do j=1,3
		xxt(j,i)=xsc(j,i)
		xxq(j,i)=xsct(j,I)
		end do
	end do
	ind=lenf2
	dij=0.
	ic=0
	do i=1,ind-1
	do j=i+1,ind
	d1=0.
	d2=0.
	ic=ic+1
		do l=1,3
		D1=D1+(XXT(L,I)-XXT(L,J))**2
		D2=D2+(XXq(L,I)-XXq(L,J))**2
		end do
	d1=sqrt(d1)
	d2=sqrt(d2)
	dij=dij+(d1-d2)**2
	end do
	end do
	dij=sqrt(dij/ic)
	drms=dij
	ntot=ind
      call u3b(w, xxq, xxt, ntot, 1, rms, u, t, ier)
      do j=1,3
      	do l=1,3
	if(abs(u(j,l)) .lt. .000001)u(j,l)=0.
	end do
      end do
	do i=1,lenf2
		do j=1,3
		xs(j,i)=t(j)
		do l=1,3
		xs(j,i)=xs(j,i)+u(j,l)*xsct(l,i)
		end do
		end do
	end do
	rms=dsqrt(rms/ntot)
	rrms=rms
	return
      end

c.... this version copied July 1986. DO NOT REDISTRIBUTE.
c.... If you want this routine, ask Wolfgang Kabsch !!
c**** CALCULATES A BEST ROTATION & TRANSLATION BETWEEN TWO VECTOR SETS
c**** SUCH THAT U*X+T IS THE CLOSEST APPROXIMATION TO Y.
c**** THE CALCULATED BEST SUPERPOSITION MAY NOT BE UNIQUE AS INDICATED
c**** BY A RESULT VALUE IER=-1. HOWEVER IT IS GARANTIED THAT WITHIN
c**** NUMERICAL TOLERANCES NO OTHER SUPERPOSITION EXISTS GIVING A
c**** SMALLER VALUE FOR RMS.
c**** THIS VERSION OF THE ALGORITHM IS OPTIMIZED FOR THREE-DIMENSIONAL
c**** REAL VECTOR SPACE.
c**** USE OF THIS ROUTINE IS RESTRICTED TO NON-PROFIT ACADEMIC
c**** APPLICATIONS.
c**** PLEASE REPORT ERRORS TO
c**** PROGRAMMER:  W.KABSCH   MAX-PLANCK-INSTITUTE FOR MEDICAL RESEARCH
c        JAHNSTRASSE 29, 6900 HEIDELBERG, FRG.
c**** REFERENCES:  W.KABSCH   ACTA CRYST.(1978).A34,827-828
c           W.KABSCH ACTA CRYST.(1976).A32,922-923
c
c  W    - W(M) IS WEIGHT FOR ATOM PAIR  c M           (GIVEN)
c  X    - X(I,M) ARE COORDINATES OF ATOM c M IN SET X       (GIVEN)
c  Y    - Y(I,M) ARE COORDINATES OF ATOM c M IN SET Y       (GIVEN)
c  N    - N IS NUMBER OF ATOM PAIRS             (GIVEN)
c  MODE  - 0:CALCULATE RMS ONLY              (GIVEN)
c      1:CALCULATE RMS,U,T   (TAKES LONGER)
c  RMS   - SUM OF W*(UX+T-Y)**2 OVER ALL ATOM PAIRS        (RESULT)
c  U    - U(I,J) IS   ROTATION  MATRIX FOR BEST SUPERPOSITION  (RESULT)
c  T    - T(I)   IS TRANSLATION VECTOR FOR BEST SUPERPOSITION  (RESULT)
c  IER   - 0: A UNIQUE OPTIMAL SUPERPOSITION HAS BEEN DETERMINED(RESULT)
c     -1: SUPERPOSITION IS NOT UNIQUE BUT OPTIMAL
c     -2: NO RESULT OBTAINED BECAUSE OF NEGATIVE WEIGHTS W
c      OR ALL WEIGHTS EQUAL TO ZERO.
c
c-----------------------------------------------------------------------
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
c**** REDUCE CUBIC TO STANDARD FORM Y**3HY+2G=0 BY PUTTING X=Y+SPUR
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
	SUBROUTINE DP(NSEQ1,NSEQ2,RESULT,ils,ilen)
C
C  Dynamic Programming for Sequence Alignment, Fortran version
C
	PARAMETER(NMAX=1600)
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
	PARAMETER(NMAX=1600)
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
	score(nseq1,i)=-20
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
	PARAMETER(NMAX=1600)
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
      PARAMETER     (MAXSEQ =1000)
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
