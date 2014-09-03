c	secassign. assigns the dominant secondary secondary elements 
c	based on the r14 criterion:
c	*********************************************************
c	basedon a smoothed distribtion uses the new 94 database*
c	*               ener14g.f 	3/2/93			*
c	*	calculates the r14 distribution function in the *
c	* 	database					*
c	*	for the 311 lattice but with a smoothed 	*
c	*	distribution					*
c	*							*
c	*********************************************************
c       the program is in `~skolnick/linker/linker42sum.f
c       a. kolinski, j. skolnick, a. godzik and w-p hu.  a method
c       for the prediction of surface u-turns and transglobular
c       connections in small proteins.  proteins 1997:27: 290-308.
c       w-p. hu, a. kolinski and j. skolnick.  improved method for
c       the prediction of the protein backbone u-turn positions
c       and the major secondary structures between the u-turns.
c       proteins 1997:29: 443-460

	parameter(nmax=2000)
	parameter(ibins=6)
	dimension xa(3,nmax),icode(nmax)	
	character*3 aa(20)
	character*1 a1
	character*6 name
	common/is/ibb(-87:91)	
	data aa/ 'gly','ala','ser','cys','val','thr','ile','pro',
     &		     'met','asp','asn','leu',
     &		     'lys','glu','gln','arg',
     &		     'his','phe','tyr','trp'/


        do i=-86,91
	   if(i .lt. -56)			ibb(i)=4
	   if(i .ge. -56 .and. i .lt. -25)	ibb(i)=3
	   if(i .ge. -25 .and. i .lt. 0)	ibb(i)=3
	   if(i .ge. 0 .and. i .lt.  26)	ibb(i)=2
	   if(i .ge. 26  .and. i .lt.  56)	ibb(i)=3
	   if(i .ge. 56)			ibb(i)=4
        end do
					ibb(0)=0
	open(unit=20,file='list.longer5')
	read(20,*)ndata
	ntot=0
	do ik=1,ndata
	   read(20,1)name
c	open(unit=6,file=name//'.link')
	   open(unit=16,file=name//'.blocksb')	
 1	   format(a6)
	   call pdb(name,xa,nres,icode,a1) !reads a pdb file and extract ca coordinates
c	write(6,*)'no of residues =',nres
	   if(nres .gt. nmax) then
c	write(6,121)name,nres
 121	      format(1x,a6, 3x,i4,'nside exceeds nmax abort program')
	      stop
	   end if		
	   open(unit=23,file=name//'.sec')
	   call r14s(nres,icode,xa)
	   close(23)
	   close(6)
	   close(16)
	end do	
	close(20)
	stop
	end
	subroutine r14s(nres,icode,xa)
c	all  dimensions are in lattice units in this program
	parameter(nmax=2000)

	logical qside,itest,iaddl,iaddr
c       sed stores the assignment of a residue as turn
	integer sed(-2:nmax) 
	dimension qside(nmax)	
        character*5 struct(4)	
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
	open(unit=11,file='helix')
c	write(6,*)'helix template:'
	do i=1,5
	read(11,*)id,(xxt(j,i),j=1,3)
	do j=1,3
	xxt(j,i)=xxt(j,i)/1.22
	end do
c	write(6,*)i,(xxt(j,i),j=1,3)
	end do
	close(11)

c       weights for smoothing CA coordinates
	awt(5)=1./243.
	awt(4)=5./243.
	awt(3)=15./243.
	awt(2)=30./243.
	awt(1)=45./243.
	awt(0)=51./243

	do i=nres+1,nmax
	   sed(i)=0
	   qside(i)=.false.
	end do

	do i=1,nres
	   sed(i)=0
	qside(i)=.true.
	end do

c       1.22angstroms is the unit distance
	ca=3.785/1.22
	ca2=ca*ca
c
c
c	cos(180-valence_angle) is in the range [-.30 0.901] including.
c	which corresponds to the range of val_ang 72.5 to 154
c


	
	do i=2,nres-1
	   r21=0	
	   r22=0.
	   dp=0.
	   do j=1,3
	      r22=r22+(xa(j,i+1)-xa(j,i))**2
	      r21=r21+(xa(j,i)-xa(j,i-1))**2
	      dp=dp+(xa(j,i+1)-xa(j,i))*(xa(j,i)-xa(j,i-1))
	   end do
	   if( abs(r21/ca2)-1. .gt. 0.2)then
	      qside(i)=.false.
	      go to 101
	   end if

c	aaa is (180-valance angle)
	   aaa=dp/sqrt(r22*r21)
	   if(aaa.lt.-0.30 .or. aaa.gt.0.901) then
	      qside(i)=.false.
	   else
	      qside(i)=.true.	
	   end if
101	   continue
	end do

c 	


	do 140 i=2,nres-2
	   ip1=i+1
	   if(qside(i).and. qside(ip1))then
	      im1=i-1
	      ip2=i+2
	      r14=(xa(1,im1)-xa(1,ip2))**2 +(xa(2,im1)-xa(2,ip2))**2 +
     &		(xa(3,im1)-xa(3,ip2))**2 
c
c	check that r14 is not unreasonable

	      ir14=int(r14+0.5)
	      if(ir14 .gt. 90)ir14=90
	      wx1=xa(1,i)-xa(1,im1)
	      wy1=xa(2,i)-xa(2,im1)
	      wz1=xa(3,i)-xa(3,im1)
	      wx2=xa(1,ip1)-xa(1,i)
	      wy2=xa(2,ip1)-xa(2,i)
	      wz2=xa(3,ip1)-xa(3,i)
	      px=wy1*wz2-wy2*wz1
	      py=wz1*wx2-wz2*wx1
	      pz=wx1*wy2-wx2*wy1

	      wx3=xa(1,ip2)-xa(1,ip1)
	      wy3=xa(2,ip2)-xa(2,ip1)
	      wz3=xa(3,ip2)-xa(3,ip1)
	      ihand=px*wx3+py*wy3+pz*wz3
	      if(ihand .lt. 0) ir14=-ir14
	      if(ir14 .lt. -86)ir14=-86
	      i14=ibb(ir14)
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
	      sed(i)=i14
12	      continue	
	   end if
140	continue
	iconf(1)=iconf(2)
	iconf(nres-1)=iconf(nres-2)
	iconf(nres)=iconf(nres)
c	set up the location of the secondary elements:
c	renumber if there are less than 3 residues per beta state

c ***create average tube:

c store coordinates in xb, which is an extended version of xa,
c since it allows insertion of extra residues at the beginning
c and at the end
	do i=1,nres
	   do j=1,3
	      xb(j,i)=xa(j,i)
	   end do
	end do

c add virtual residues at beginning
	do i=-4,0
	   do j=1,3
	      xb(j,i)=xb(j,1)
	   enddo
	end do

c       add virtual residues at end
	do i=1,5
	   do j=1,3
	      xb(j,i+nres)=xb(j,nres)
	   end do
	end do

c       create smoothed backbone, by averaging over coordinates
c       of preceding and succeeding five residues
	do i=1,nres
	   do j=1,3
	      y(j,i)=0.
	      do iw=-5,5
c                iwa is associated weigth
		 iwa=abs(iw) 
		 y(j,i)=y(j,i)+awt(iwa)*xb(j,i+iw)
	      end do
	   end do
	end do

c       substitute backbone with smoothed chain
	do i=1,nres
	   do j=1,3
	      xa(j,i)=y(j,i)
	   end do
	end do

c       calculate radius of curvature passing through
c       each residue, store in rad()
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
c	   write(6,*)'rc**_1',i,rad(i)
c          sed() is signal for loop presence
	   if(sum .lt. .09)then
	      sed(i)=0
	   else
	      sed(i)=1
	   end if
 22	   format(i4,3f6.4)
	end do

	icount=0
100 	continue

c ***calculate number of loops (nloop), and residues
c    where each loop begins mbegin() and ends mend()

c       No loop can be defined in the first three and the
c       last three residues in the chain
	sed(0)=0
	sed(1)=0
	sed(2)=0
	sed(3)=0
	sed(nres-2)=0
	sed(nres-1)=0
	sed(nres)=0
	icount=icount+1
	if(icount .gt. 2) goto 200
	nloop=0
	do i=1,nres
c	if(sed(i).eq.1 .and. sed(i-1) .ne. 1 .and. sed(i-2) .ne.1)then
c       the if signals the possible beginning of a loop
	   if(sed(i).eq.1 .and. sed(i-1) .eq. 0)then
	      nloop=nloop+1
	      istate(nloop)=1
	      mbegin(nloop)=i
	      mend(nloop)=i
c 	elseif(sed(i).eq.1.and.sed(i+1).ne.4.and.sed(i+2).ne.4)then
c       the if signals the possible end of a loop
	   elseif(sed(i).eq.1 .and.sed(i+1).eq.0)then
	      mend(nloop)=i
	   end if
	end do

c	write(6,*)'no. of loops =', nloop
c	do k=1,nloop
c	write(6,*)'begin=',mbegin(k),'  end=',mend(k)
c	end do


c ***find number of secondary structure elements, nelmt,
c    and where each element begins, lbegin(), and ends, lend()
	nelmt=0
	mend(0)=0
	mbegin(nloop+1)=nres
	do k=1,nloop
c          the if signals a loop bigger than one residue
	   if(mbegin(k)-mend(k-1).gt. 1)then
	      nelmt=nelmt+1
	      lbegin(nelmt)=mend(k-1)+1
	      lend(nelmt)=mbegin(k)-1
	      lstate(nelmt)=istate(k)
	      if(lend(nelmt) .eq. nres-1)lend(nelmt)=nres-2
	   else
c             a loop of one residue is no loop. Thus, extend
c             the end of the previous element to include the
c             one residue loop
	      lend(nelmt)=mend(k)
	   end if
	end do
c	write(6,*)'the no. of secondary structural elements=',nelmt

	lend(0)=0
	lbegin(nelmt+1)=mend(nloop)+1
	lend(nelmt+1)=nres
	do k=1,nelmt+1
	   ib=lbegin(k)
	   ie=lend(k)
	   a=0.
	   is=0
	   do ii=ib+1,ie
	      is=is+1
	      do j=1,3
		 a=a+(xa(j,ii)-xa(j,ii-1))**2
	      end do
	   end do
	   a=sqrt(a)/is
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
	ib=lbegin(k)
	ie=lend(k)
	n=ie-ib+1
c	write(6,*)k,'length=',ncheck
	if(lstate(k) .eq.2)then
	ncheck=15
	elseif(lstate(k).eq.4)then
	ncheck=11
	end if
	if(n .gt. ncheck)then
	do i=ib,ie	
	if(i .gt.4 .and. i .lt.nres-3)then
c	write(6,*)i,rad(i)
	if(rad(i) .gt. 0.07)then
	sed(i)=1
	itest=.true.
	end if
	end if
	end do
	end if
	end do
	if(itest)goto 100
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
	if(iconf(ii) .ne.is)then
	lbegin(k)=ii+1
	iaddl=.false.
	else
	go to 105
	end if
	end do
105	continue
	do ii=ie,ib+1,-1
	if(iconf(ii) .ne.is)then
		lend(k)=ii-1
	iaddr=.false.
		else
		go to 106
		end if
		end do
106 	continue
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
103	continue	
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
104	continue
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
	write(16,*)nbk,'  blocks'
c	write(6,*)nbk,'  blocks'
	do i=1,nres
	jsec(k)=3
	end do
	do k=1,nbk
c	write(6,*)k,' ',kbegin(k),' ',kend(k),kstate(k)
	write(16,*)k,' ',kbegin(k),' ',kend(k),kstate(k)
	do ik=kbegin(k),kend(k)
	jsec(ik)=kstate(k)
	end do
	end do
	write(23,*)nres
	write(23,*)(jsec(i),i=1,nres)
	write(23,*)(icode(i),i=1,nres)	
	
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

	subroutine pdb(name,xa,nres,icode,a1)
c	==================================================
c       this program reads in a pdb file extracts the ca =
c	coordinates					 =
c	all units are in lattice units(1lu=1.22 angstroms=
c	==================================================
	parameter(cfs=0.8197,nmax=2000)
	character*6 name
	character*4 itype
	character*1 a1,chain
	character*3 ares,aa(20)
	character*6 atom
	dimension x(3),xa(3,nmax),icode(nmax)
	data aa/ 'gly','ala','ser','cys','val','thr','ile','pro',
     &		     'met','asp','asn','leu',
     &		     'lys','glu','gln','arg',
     &		     'his','phe','tyr','trp'/
	do i=1,nmax
	icode(i)=0
		do j=1,3
		xa(j,i)=0
		end do
	end do
c	write(6,*)name//'.pdb is being analyzed'
	icount=0
	open(unit=10,file='/nfs/users/yzhang/for1/template_lib/'//
     &	name//'.pdb')
	ntest=0
	nres=0
	inumold=-200
	ntest2=0
	icount=0	
1	continue
	ntest=ntest+1
c	if(ntest .gt.nmax) go to 20
100	read(10,10,err=1,end=20)atom,id,itype,
     &	ares,chain,inum,x(1),x(2),x(3)
10	format(a6,i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3)
c 	if(atom.eq. 'END' .or. atom .eq. 'TER')then
	if(atom.eq. 'END' )then
c		write(6,*)'end of file'
		go to 20
	end if
	if(atom .eq. 'ENDMDL')go to 20
	if(atom .eq.'TER')go to 20	
	if(atom .eq. 'HETATM')go to 1
	if(ares .eq. 'ACE')go to 1
	if(atom .ne. 'ATOM  ')go to 1
c	begin read in of new residue
	if(itype .eq. ' CA')then	
	icount=icount+1
	xa(1,icount)=x(1)*cfs
	xa(2,icount)=x(2)*cfs
	xa(3,icount)=x(3)*cfs
		do ires=1,20
			if(ares .eq. aa(ires))then
			icode(icount)=ires
			go to 100
			end if
		end do
	end if
	go to 100
20	continue
200	continue
	nres=icount
c	write(6,*)name,' which has',nres,' residues is finished'
	close(10)
	return
	end
	

	subroutine rmscalc(ind,xxq,xxt,rms)
c this program will compare (on the rms level) the pdb file and
c output from lattice simulations 
c
c               program under development
c
c            version 1.0 (19 novemebr 1990)
c
c   written in novebmer 1990, based of programs findsite
c   and various generations of earlier programs (search, kontact etc)
c
c  rember to upgrade version number after introducing changes !!!
c      1.0 - just started
c---------------------------------------------------
c
c                      input section
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
c this subroutine will read the job informations for the findsite progra
cm
c all parameters are passed in commons
c.... this version copied july 1986. do not redistribute.
c.... if you want this routine, ask wolfgang kabsch !!
c**** calculates a best rotation & translation between two vector sets
c**** such that u*x+t is the closest approximation to y.
c**** the calculated best superposition may not be unique as indicated
c**** by a result value ier=-1. however it is garantied that within
c**** numerical tolerances no other superposition exists giving a
c**** smaller value for rms.
c**** this version of the algorithm is optimized for three-dimensional
c**** real vector space.
c**** use of this routine is restricted to non-profit academic
c**** applications.
c**** please report errors to
c**** programmer:  w.kabsch   max-planck-institute for medical research
c        jahnstrasse 29, 6900 heidelberg, frg.
c**** references:  w.kabsch   acta cryst.(1978).a34,827-828
c           w.kabsch acta cryst.(1976).a32,922-923
c
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs             (given)
c  mode  - 0:calculate rms only              (given)
c      1:calculate rms,u,t   (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs        (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier   - 0: a unique optimal superposition has been determined(result)
c     -1: superposition is not unique but optimal
c     -2: no result obtained because of negative weights w
c      or all weights equal to zero.
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
c**** determine centroids of both vector sets x and y
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
c**** determine correlation matrix r between vector sets y and x
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
c**** calculate determinant of r(i,j)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** form upper triangle of transposed(r)*r
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** eigenvalues *****************************************
c**** form characteristic cubic  x**3-3*spur*x**2+3*cof*x-det=0
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
c**** reduce cubic to standard form y**3-3hy+2g=0 by putting x=y+spur
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** solve cubic. roots are e1,e2,e3 in decreasing order
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
c.....handle special case of 3 identical roots
c 224 "rms.for"
      if (mode) 10, 50, 10
c**************** eigenvectors *****************************************
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
c****************** rotation matrix ************************************
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
c****************** translation vector *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** rms error **************************************
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
c.....end u3b...........................................................
c----------------------------------------------------------
c                       the end
c----------------------------------------------------------
c 338 "rms.for"
      end

