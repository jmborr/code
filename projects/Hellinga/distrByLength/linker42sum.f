c
c       pgf77 -Mextend -O -s -fast -Wl,-static -o linker42sum.x linker42sum.f
c	secassign. Assigns the dominant secondary secondary elements 
c	based on the r14 criterion:
C	*********************************************************
C	BASEDON A SMOOTHED DISTRIBTION USES THE NEW 94 database*
C	*               ENER14G.F 	3/2/93			*
C	*	CALCULATES THE R14 DISTRIBUTION FUNCTION IN THE *
C	* 	DATABASE					*
C	*	FOR THE 311 LATTICE BUT WITH A SMOOTHED 	*
C	*	DISTRIBUTION					*
C	*							*
C	*********************************************************
	PARAMETER(nmax=2000)
	PARAMETER(ibins=6)
	DIMENSION XA(3,NMAX),ICODE(NMAX)	
	CHARACTER*3 AA(20)
	CHARACTER*1 A1
	CHARACTER*255 NAME
	character*5 name5
	common/is/ibb(-87:91)	
	data aa/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &		     'MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP'/


        do i=-86,91
        if(i .lt. -56)			ibb(i)=4
        if(i .ge. -56 .and. i .lt. -25)	ibb(i)=3
        if(i .ge. -25 .and. i .lt. 0)	ibb(i)=3
        if(i .ge. 0 .and. i .lt.  26)	ibb(i)=2
        if(i .ge. 26  .and. i .lt.  56)	ibb(i)=3
        if(i .ge. 56)			ibb(i)=4
        end do
					ibb(0)=0

	OPEN(unit=20,file='/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/filtered/LIST.small')
	open(unit=6,file='out_linker42sum')
	READ(20,*)NDATA
	NTOT=0
	open(unit=23,file='secmaster')
	DO IK=1,NDATA
        READ(20,1)NAME
1       format(a255)
        do i=255,1,-1
        if(name(i:i) .eq. ' ')nsize=i
        end do
        nsize=nsize-1
	a1='_'
	open(unit=10,file='/gpfs1/scratch/jose/Hellinga/distrByLengthDecreasedLoopLoopInteraction/linker42sum/filtered/'//name(1:nsize))
	CALL  PDB(NAME,XA,nres,icode,a1)
	write(6,*)'NO of residues =',nres

	call r14s(nres,icode,xa,name)

	close(6)
	close(16)
	END DO	
	CLOSE(20)
	close(23)
	stop
	end



	SUBROUTINE R14S(NRES,ICODE,XA,name)
C	ALL  DIMENSIONS ARE IN LATTICE UNITS IN THIS PROGRAM
	PARAMETER(nmax=2000)

	LOGICAL QSIDE,itest,iaddl,iaddr
	integer sed(-2:nmax)
	DIMENSION QSIDE(NMAX)	
        character*5 struct(4)	
	character struct1(4)
	character*255 name
	common/is/ibb(-87:91)
	dimension icode(nmax),mbegin(0:nmax)
	dimension mend(0:nmax),lbegin(40),lend(0:40)
	dimension istate(40),lstate(40)
	dimension iconf(0:nmax)
	dimension y(3,-7:nmax),il(nmax),ir(nmax)
	dimension rc(3),xa(3,nmax),xb(3,-10:nmax)
 	dimension rad(nmax),ased(nmax)
	dimension x(3),z(3)
	dimension awt(0:5),jsec(nmax)
	double precision xxt(3,60),xxq(3,60),rms
	dimension kbegin(40),kend(40),kstate(40)
	data struct /' coil','helix',' turn',' beta'/
	data struct1 /'-','H','-','E'/

	OPEN(Unit=11,file='helix')
c	write(6,*)'helix template:'
	do i=1,5
	read(11,*)id,(xxt(j,i),j=1,3)
	do j=1,3
	xxt(j,i)=xxt(j,i)/1.22
	end do
c	write(6,*)i,(xxt(j,i),j=1,3)
	end do
	close(11)
	awt(5)=1./243.
	awt(4)=5./243.
	awt(3)=15./243.
	awt(2)=30./243.
	awt(1)=45./243.
	awt(0)=51./243
	DO I=nres+1,nmax
	sed(i)=0
	qside(i)=.false.
	end do
	DO I=1,NRES
	sed(i)=0
	qside(i)=.true.
	end do
	CA=3.785/1.22
	CA2=CA*CA
C
c
c	cos(180-valence_angle) is in the range [-.30 0.901] including.
c	which corresponds to the range of val_ang 72.5 to 154
c


	
	DO I=2,NRES-1
	R21=0	
	R22=0.
	DP=0.
		DO J=1,3
		R22=R22+(XA(J,I+1)-XA(J,I))**2 ! (v(i+1)-v(i))^2
		R21=R21+(XA(J,I)-XA(J,I-1))**2 ! (v(i)-v(i-1))^2
		DP=DP+(XA(J,I+1)-XA(J,I))*(XA(J,I)-XA(J,I-1))  ! (v(i+1)-v(i))*(v(i)-v(i-1))
		END DO
		IF( ABS(R21/CA2)-1. .GT. 0.2)THEN
			QSIDE(I)=.FALSE.
			GO TO 101
		END IF

C	AAA IS (180-VALANCE ANGLE)
		AAA=DP/SQRT(R22*R21)
		IF(AAA.LT.-0.30 .OR. AAA.GT.0.901) THEN
		QSIDE(I)=.FALSE.
		ELSE
		QSIDE(I)=.TRUE.	
		END IF
101 	CONTINUE
	END DO

c 	


	do 140 i=2,nres-2
	ip1=i+1
	if(qside(i).and. qside(ip1))then
	im1=i-1
	ip2=i+2
	r14=(xa(1,im1)-xa(1,ip2))**2 +(xa(2,im1)-xa(2,ip2))**2 +
     &		(xa(3,im1)-xa(3,ip2))**2 !(v(i+2)-v(i-1))^2
C
C	CHECK THAT R14 IS NOT UNREASONABLE

	IR14=INT(R14+0.5)
	IF(IR14 .GT. 90)IR14=90
	WX1=XA(1,I)-XA(1,IM1) !w1=v(i)-v(i-1)
	WY1=XA(2,I)-XA(2,IM1)
	WZ1=XA(3,I)-XA(3,IM1)
	WX2=XA(1,IP1)-XA(1,I) !w2=v(i+1)-v(i)
	WY2=XA(2,IP1)-XA(2,I)
	WZ2=XA(3,IP1)-XA(3,I)
	PX=WY1*WZ2-WY2*WZ1    !p=w1xw2
	PY=WZ1*WX2-WZ2*WX1
	PZ=WX1*WY2-WX2*WY1

	WX3=XA(1,IP2)-XA(1,IP1) !w3=v(i+2)-v(i+1)
	WY3=XA(2,IP2)-XA(2,IP1)
	WZ3=XA(3,IP2)-XA(3,IP1)
	IHAND=PX*WX3+PY*WY3+PZ*WZ3 !p*w3
	IF(IHAND .LT. 0) IR14=-IR14
	IF(Ir14 .lt. -86)IR14=-86
	I14=ibb(IR14)
	iconf(i)=i14
	if(i .lt. nres-2 .and. i.gt.2)then
	iv=0
	do id=i-2,i+2
	iv=iv+1
	do j=1,3
	xxq(j,iv)=xa(j,id)
	end do
	end do
	ind=iv
	call rmscalc(ind,xxq,xxt,rms)
	if(rms .lt.0.25)then
	iconf(i)=2
	else
	if(iconf(i).eq.2)iconf(i)=3
	end if
	end if

c	write(6,*)i,ir14,rms
	SED(i)=I14
12 	CONTINUE	
	end if
140 	continue
	iconf(1)=iconf(2)
	iconf(nres-1)=iconf(nres-2)
	iconf(nres)=iconf(nres)
c	set up the location of the secondary elements:
c	renumber if there are less than 3 residues per beta state

c	create average tube:

	do i=1,nres
	   do j=1,3
	      xb(j,i)=xa(j,i)
	   end do
	end do

	do i=-4,0
	   do j=1,3
	      xb(j,i)=xb(j,1)
	   enddo
	end do

	do i=1,5
	   do j=1,3
	      xb(j,i+nres)=xb(j,nres)
	   end do
	end do

	do i=1,nres
	   do j=1,3
	      y(j,i)=0.
	      do iw=-5,5
		 iwa=abs(iw)
		 y(j,i)=y(j,i)+awt(iwa)*xb(j,i+iw)
	      end do
	   end do
	end do

	do i=1,nres
	   do j=1,3
	      xa(j,i)=y(j,i)
	   end do
	end do

	do i=3,nres-3	
	   rad(i)=0.
	   sum=0.
	   anorm=0.
	   do j=1,3
	      anorm=anorm+(y(j,i+1)-y(j,i-1))**2
	      rc(j)=y(j,i-2)+y(j,i+2)-2.*y(j,i)
	      sum=sum+rc(j)**2/16.
	   end do
	   anorm=sqrt(anorm)
	   sum=sum/(anorm+0.0001)
	   rad(i)=sum
c	write(6,*)'rc**_1',i,rad(i)
	   if(sum .lt. .09)then
	      sed(i)=0
	   else
	      sed(i)=1
	   end if
 22	   format(i4,3f6.4)
	end do


	icount=0
100 	continue
	SED(0)=0
	SED(1)=0
	SED(2)=0
	sed(3)=0
	sed(nres-2)=0
	sed(nres-1)=0
	sed(nres)=0
	icount=icount+1
	if(icount .gt. 2) goto 200 !maximum number of iterations
	nloop=0
	do i=1,nres
c	if(sed(i).eq.1 .and. sed(i-1) .ne. 1 .and. sed(i-2) .ne.1)then
	   if(sed(i).eq.1 .and. sed(i-1) .eq. 0)then
	      nloop=nloop+1
	      istate(nloop)=1
	      mbegin(nloop)=i
	      mend(nloop)=i
c 	elseif(sed(i).eq.1.and.sed(i+1).ne.4.and.sed(i+2).ne.4)then
	   elseif(sed(i).eq.1 .and.sed(i+1).eq.0)then
	      mend(nloop)=i
	   end if
	end do
c	write(6,*)'no. of loops =', nloop
c	do k=1,nloop
c	write(6,*)'begin=',mbegin(k),'  end=',mend(k)
c	end do

	nelmt=0 !number of secondary structure elements
	mend(0)=0
	mbegin(nloop+1)=nres
	do k=1,nloop
	   if(mbegin(k)-mend(k-1).gt. 1)then
	      nelmt=nelmt+1
	      lbegin(nelmt)=mend(k-1)+1
	      lend(nelmt)=mbegin(k)-1
	      lstate(nelmt)=istate(k)
	      if(lend(nelmt) .eq. nres-1)lend(nelmt)=nres-2
	   else
	      lend(nelmt)=mend(k)
	   end if
	end do
c	write(6,*)'the no. of secondary structural elements=',nelmt

	lend(0)=0
	lbegin(nelmt+1)=mend(nloop)+1
	lend(nelmt+1)=nres
	do k=1,nelmt+1
	   ib=lbegin(k) !beginning of element
	   ie=lend(k)   !end of element
	   a=0.
	   is=0
	   do ii=ib+1,ie
	      is=is+1
	      do j=1,3
		 a=a+(xa(j,ii)-xa(j,ii-1))**2
	      end do
	   end do
	   a=sqrt(a)/is !average CA_{i}-CA_{i+1} distance
	   if(a .lt. 0.6)then
c	write(6,*)'element=',k,a,ib,ie,' helix'
	      lstate(k)=2
	   else if(a .lt. 1.1)then
c	write(6,*)'element=',k,a,ib,ie,' beta'
	      lstate(k)=4
	   else
c	write(6,*)'element=',k,a,ib,ie,'extended  beta/turn'
	      lstate(k)=4
	   end if
	end do

	itest=.false.	
	do k=1,nelmt+1
c	   write(*,*)k
	   ib=lbegin(k)
	   ie=lend(k)
	   n=ie-ib+1
c	write(6,*)k,'length=',ncheck
	   if(lstate(k) .eq.2)then
	      ncheck=25 !maximum length of helix
	   elseif(lstate(k).eq.4)then
	      ncheck=11 !maximum length of strand
	   end if
	   if(n .gt. ncheck)then
	      do i=ib,ie	
		 if(i .gt.4 .and. i .lt.nres-3)then
c	write(6,*)i,rad(i)
		    if(rad(i) .gt. 0.07)then
		       sed(i)=1 !break the long helix or strand
		       itest=.true.
		    end if
		 end if
	      end do
	   end if
	end do
	if(itest)goto 100 !begin again the apportioning of the loops
200 	continue

	do i=1,nres
	   jsec(i)=3
	end do

	lbegin(nelmt+2)=nres
	do k=1,nelmt+1
	   ib=lbegin(k)
	   ie=lend(k)
c	write(6,*)k,'initialset ',ib,ie
	   is=lstate(k)
	   iaddr=.true.
	   iaddl=.true.
c	end of test for addition
c	now we begin to delete states
	   do ii=ib,ie-1
	      if(iconf(ii).ne.is)then
		 lbegin(k)=ii+1
		 iaddl=.false.
	      else
		 go to 105
	      end if
	   end do
 105	   continue
	   do ii=ie,ib+1,-1
	      if(iconf(ii) .ne.is)then
		 lend(k)=ii-1
		 iaddr=.false.
	      else
		 go to 106
	      end if
	   end do
 106	   continue
c	addition of helical states
	   if(iaddl)then
	      if(ib -lend(k-1).gt.1)then
		 do iv=ib-1,lend(k-1)+1,-1
		    if(iconf(iv).eq.is)then
		       lbegin(k)=iv
		    else
		       go to 103
		    end if
		 end do
	      end if
	   end if
 103	   continue	
	   if(iaddr)then
	      if(lbegin(k+1) -ie.gt.1)then
		 do iv=ie+1,lbegin(k+1)-1
		    if(iconf(iv).eq.is)then
		       lend(k)=iv
		    else
		       go to 104
		    end if
		 end do
	      end if
	   end if
 104	   continue
c	write(6,*)k,': final set ',lbegin(k),lend(k)
	end do
	nbk=0
	do k=1,nelmt+1
	   if(lstate(k) .eq.2 .or. lstate(k).eq.4)then
	      idiff=lend(k)-lbegin(k)
	      if(idiff .gt.1)then
		 nbk=nbk+1
		 kend(nbk)=lend(k)
		 kbegin(nbk)=lbegin(k)
		 kstate(nbk)=lstate(k)
	      else
c	write(6,*)'eliminating block ',k, 'too short',
c     &	lbegin(k),lend(k)
	      end if
	   end if
	end do
c	write(16,*)nbk,'  blocks'
c	write(6,*)nbk,'  blocks'
	do i=1,nres
	   jsec(k)=3
	end do
	do k=1,nbk
c	write(6,*)k,' ',kbegin(k),' ',kend(k),kstate(k)
c	write(16,*)k,' ',kbegin(k),' ',kend(k),kstate(k)
	   do ik=kbegin(k),kend(k)
	      jsec(ik)=kstate(k)
	   end do
	end do
c	write(23,233)nres,name
	write(23,'a,5a')'>',name(1:5)
 233	format(i4,1x,a9)
c	write(23,53)(jsec(i),i=1,nres)
	write(23,'60a') (struct1(jsec(i)),i=1,nres)
	write(23,'a') ''
 53	format(50(i1,1x))

c	all elements are assigned
	do i=1,nres
	   ased(i)=0.
	end do
	return
	end 	


	subroutine dot(x,z,dp)
	dimension x(3),z(3)
	dp=0.
	do j=1,3
	dp=dp+x(j)*z(j)
	end do
	return
	end

	SUBROUTINE PDB(NAME,XA,nres,icode,a1)
c	==================================================
c       this program reads in a pdb file extracts the ca =
c	coordinates					 =
c	all units are in lattice units(1lu=1.22 Angstroms=
c	==================================================
	PARAMETER(CFS=0.8197,nmax=2000)
	CHARACTER*4 NAME,ITYPE
	CHARACTER*1 A1,chain
	character*3 ARES,AA(20)
	character*6 atom
	DIMENSION X(3),XA(3,NMAX),ICODE(NMAX)
	data aa/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &		     'MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP'/
	DO I=1,NMAX
	ICODE(I)=0
		DO J=1,3
		XA(J,I)=0
		END DO
	END DO
	write(6,*)NAME//'.pdb is being analyzed'
	ICOUNT=0
c	OPEN(UNIT=10,FILE='/library/pdb/pdb'//NAME//'.ent')
	NTEST=0
	NRES=0
	INUMOLD=-200
	ntest2=0
	ICOUNT=0	
1	CONTINUE
	ntest=ntest+1
c	if(ntest .gt.nmax) go to 20
100	READ(10,10,ERR=1,END=20)ATOM,ID,ITYPE,
     &	ARES,CHAIN,INUM,X(1),X(2),X(3)
10	FORMAT(A6,I5,1X,A4,1x,a3,1X,a1,I4,1x,3X,3F8.3)
c 	IF(ATOM.EQ. 'END' .OR. ATOM .EQ. 'TER')THEN
	IF(ATOM.EQ. 'END' )THEN
		WRITE(6,*)'END OF FILE'
		GO TO 20
	END IF
	if(a1 .ne.'_')then
	if(CHAIN .ne. a1)go to 1
	if(CHAIN .eq.a1)then
	if(atom .eq. 'ENDMDL')go to 20
	if(atom .eq.'TER')go to 20	
	end if
	else
	if(atom .eq.'TER')go to 20
	if(atom .eq. 'ENDMDL')go to 20
	end if
	IF(ATOM .EQ. 'HETATM')GO TO 1
	if(ares .EQ. 'ACE')GO TO 1
	IF(ATOM .NE. 'ATOM  ')GO TO 1
C	BEGIN READ IN OF NEW RESIDUE
	IF(ITYPE .EQ. ' CA')THEN	
	ICOUNT=ICOUNT+1
	XA(1,ICOUNT)=X(1)*CFS
	XA(2,ICOUNT)=X(2)*CFS
	XA(3,ICOUNT)=X(3)*CFS
	DO IRES=1,20
	   IF(ARES .EQ. AA(IRES))THEN
	      ICODE(ICOUNT)=IRES
	      GO TO 100
	   END IF
	END DO
	END IF
	go to 100
20	CONTINUE
200	continue
	NRES=ICOUNT
	write(6,*)name,' WHICH HAS',NRES,' RESIDUES is finished'
	close(10)
	RETURN
	END
	
	subroutine rmscalc(ind,xxq,xxt,rms)
c this program will compare (on the rms level) the PDB file and
c output from lattice simulations 
c
c               program under development
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
c
c all variables with names ending in q are for query protein
c                                    t are for target protein
c maximum number of atoms and residues are defined as parameters
      	double precision w(2000),u(3, 3), t(3), rms, sigma
	double precision xxq(3,60),xxt(3,60)
	dimension rms2(4000)
    	 version = 1.0

      do i = 1, 2000
        w(i) = 1.0
      end do

	ntot=ind
      call u3b(w, xxq, xxt, ntot, 0, rms, u, t, ier)
	rms=dsqrt(rms/ntot)
	return
      
      end
c this subroutine will read the job informations for the FindSite progra
cm
c all parameters are passed in commons
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

