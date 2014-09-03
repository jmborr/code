c	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c       pgf77 -Mextend -O -s -fast -Wl,-static -o contactsev3.x contactsev3.f
c
c	seqset.f calculates
c	1. sequence files
c	2. contact map files for the fisher database
C	*********************************************************
C	*	ALL coordinates are in ANGSTROMS!!		*
C	*********************************************************
	PARAMETER(nmax=3000)
	CHARACTER*3 AA(0:21)
	CHARACTER*255 NAME
	CHARACTER*10 NAME2
	LOGICAL QSIDE	
	COMMON qside(-3:nmax)
	dimension iside(-1:21)
	COMMON /STAT/nc(nmax),icon(50,nmax),ior(50,nmax)
	DIMENSION XA(3,4,NMAX),ICODE(-NMAX:NMAX) ,XS(3,30,NMAX),
     &	SCM(3,NMAX)
	character*1 aseq1,aa1(0:19),a1
	character*3 a3
	parameter(jvec=30)	
	common/vectors/ivec,kind(jvec,jvec),ivcover(-1:19,-1:19)
	common/vectors2/avx(jvec),avy(jvec),avz(jvec)

	dimension touch(0:nmax,0:nmax),nsec(nmax),sb(3,nmax)
	common/surface/ireshist(0:21),icover(0:21,0:jvec)
	common/contacts/sidedist(0:19,0:19),backside(0:19,0:19)
	dimension idrl(0:19,0:19),idrbs(0:19,0:19)
	dimension cover(0:21,0:jvec),std(0:21,0:jvec)
	dimension ibin(0:jvec),xm(0:19),a(3),b(3)
C	========================================
C	FOR HYDROPHOBIC MOMENT UPDATE
C	========================================
	DATA ISIDE/4,1,1, 2, 2, 3, 3, 4, 3,
     &                  4,4,4,4,
     &                  5,5,5,7,
     &			6,7,8,10,2,6/

	DATA AA/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &		     'MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP','CYX','HIX'/
C     &		     'H','F','Y','W','X'/
	DATA AA1/ 'G','A','S','C','V','T','I','P',
     &		     'M','D','N','L',
     &		     'K','E','Q','R',
     &		     'H','F','Y','W'/

	
	OPEN(UNIT=6,FILE='out_contacts')
	OPEN(UNIT=20,FILE='LIST.structures')
	OPEN(UNIT=19,FILE='DISTFIT',STATUS='OLD')
	open(UNIT=9,FILE='stat_contactsev3')
c	uses optimized set of fits
	read(19,*)
	do i=0,19
	read(19,104)a3,(idrl(i,j),j=0,19)
104	format(a3,2x,20(i3,3x))
	do j=0,19
	sidedist(i,j)=idrl(i,j)**2
	sidedist(i,j)=1.6*sidedist(i,j)/(10000)
	end do
	end do
	close(19)

	READ(20,*)NDATA
	NTOT=0
	ncases=0
	DO IK=1,NDATA
	READ(20,1)NAME
1	format(a255)
        do i=255,1,-1
        if(name(i:i) .eq. ' ')nsize=i
        end do
	nsize=nsize-1

c	\\\\\\\\\\\\\\\\\\\/////////////////////////////////
	a1='_'
	
c	open(unit=10,file='/gpfs1/scratch/jose/Hellinga/filtered/'//name(1:nsize))
	open(unit=10,file='/gpfs1/scratch/jose/Hellinga/filtered/'//name(1:nsize))
	CALL PDBSIDE(XA,XS,SCM,NRES,ICODE,ibegin,a1,sb)
	close(10)
	WRITE(9,*)NAME(1:nsize)
c	\\\\\\\\\\\\\\\\\\\/////////////////////////////////
	do i=1,nres
	if(icode(i).gt.0)icode(i)=4
	end do
22	format(a5,1x,i5,1x,i5)
	close(31)
	ncases=ncases+1
	IF(NRES .GT. NMAX) THEN
		WRITE(6,121)NAME,NRES
121		FORMAT(1X,A4, 3X,I4,'NSIDE EXCEEDS NMAX ABORT PROGRAM')
		STOP
		END IF
	CALL  STATCALC(NRES,ICODE,XA,XS,SCM,sb)
	OPEN(unit=71,file='150_3040fit3/'//NAME(1:nsize)//'.FIT')
	do ii=1,nres
	if(nc(ii).gt.15)nc(ii)=15
	write(71,543)ii,nc(ii),(icon(j,ii),j=1,nc(ii))
	write(71,543)ii,nc(ii),(ior(j,ii),j=1,nc(ii))
	end do
	
543	format(26(i4,1x))
	WRITE(6,*)NAME,'.pdb is now analyzed with ',nres,'RESIDUES'	
	NTOT=NTOT+NRES
	close(31)
1000	continue
	END DO		
	CLOSE(20)
	WRITE(6,*)'       =============     '
	WRITE(6,*)' THE NUMBER OF CASES ANALYZED=',NCASES
	write(6,*)'total number of residues',ntot
	WRITE(6,*)'       =============     '
	STOP
	END
	SUBROUTINE STATCALC(NRES,ICODE,XA,XS,SCM,sb)
C	ALL  DIMENSIONS ARE IN ANGSTROMS IN THIS PROGRAM
	PARAMETER(nmax=3000,AMIN=20.25)
	LOGICAL QSIDE,TOUCH
	COMMON qside(-3:nmax),itest
	COMMON /STAT/nc(nmax),icon(50,nmax),ior(50,nmax)
	DIMENSION xa(3,4,NMAX),ICODE(-NMAX:NMAX) ,XS(3,30,NMAX)
	dimension sb(3,nmax)
	parameter(jvec=30)
	common/vectors/ivec,kind(jvec,jvec),ivcover(-1:19,-1:19)
	common/vectors2/avx(jvec),avy(jvec),avz(jvec)


	PARAMETER(ar=4.5)
c	ar is the radius of each heavy atom

	common/surface/ireshist(0:21),icover(0:21,0:jvec)
	common/contacts/sidedist(0:19,0:19),backside(0:19,0:19)
	DIMENSION ISIDE(0:21)
	DIMENSION SCM(3,NMAX),x(3)
	character*3 aa(0:21)
	DATA ISIDE/1,1, 2, 2, 3, 3, 4, 3,
     &                  4,4,4,4,
     &                  5,5,5,7,
     &			6,7,8,10,2,6/

	DATA AA/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &		     'MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP','CYX','HIX'/

	CA=3.785
	CA2=CA*CA
C
C
C	COS(180-VALENCE_ANGLE) IS IN THE RANGE [-.30 0.901] INCLUDING.
C	WHICH CORRESPONDS TO THE RANGE OF VAL_ANG 72.5 TO 154
C

	do i=-3,0
	qside(i)=.false.
	end do
	DO I=NRES+1,NMAX
	QSIDE(I)=.FALSE.
	END DO
	do i=1,nres
	qside(i)=.true.
	end do
	DO I=-NRES,NRES
C	SET UP OF THE BACKBONE INTERACTION SET:
	IF(I.LT.0)ICODE(I)=-1
C	SET UP OF THE SIDECHAIN INTERACTION SET
	END DO
	DO I=2,nres-1
	IF(qside(i))then
	R21=0	
	R22=0.
	DP=0.
		DO J=1,3
		R22=R22+(XA(J,3,I+1)-XA(J,3,I))**2
		R21=R21+(XA(J,3,I)-XA(J,3,I-1))**2
		DP=DP+(XA(J,3,I+1)-XA(J,3,I))*(XA(J,3,I)-XA(J,3,I-1))
		END DO
		IF( ABS(R21/CA2)-1. .GT. 0.2)THEN
			QSIDE(I)=.FALSE.
			GO TO 101
		END IF

C	AAA IS (180-VALANCE ANGLE)
		AAA=DP/SQRT(R22*R21)
		IF(AAA.LT.-0.30 .OR. AAA.GT.0.901) THEN
		QSIDE(I)=.FALSE.
		END IF
	end if
101 	CONTINUE
	END DO
c		
c	residues 1 and nres are undefined by the cb procedure

c	QSIDE(NRES)=.FALSE.
c	QSIDE(1)=.FALSE.

C
C	==========================================================
c	LONG  RANGE CONTACTS
	DO 140 I=1,NRES
	nc(i)=0
	IF(QSIDE(I)) THEN
	ISEQ1=ICODE(I)
	DO L=1,NRES
 	IF(QSIDE(L))THEN
	IF(IABS(L-I) .gt.2)then
		ISEQ2=ICODE(L)
C	+++++++++++++++++++++++++++++++++++++++++++++++
c	side chain - side chain contacts
	touch=.false.
	DP=0.0
		DO J=1,3
		x(j)=SCM(J,L)-SCM(J,I)
		DP=DP+x(j)**2
		END DO
		IF(DP .LT. SIDEDIST(ISEQ1,ISEQ2))THEN
		DP1=0.
		s1=0.
		s2=0.
		do j=1,3
		s1=s1+sb(j,l)**2
		s2=s2+sb(j,i)**2
		DP1=DP1+SB(J,L)*SB(J,I)
		end do	
		dp1=dp1/(sqrt(s1)*sqrt(s2))
		nc(i)=nc(i)+1
			if(dp1 .le. -0.5)then	
			icon(nc(i),i)=L
			ior(nc(i),i)=1
			elseif(dp1 .gt. -0.5 .and. dp1 .lt.0.5)then
			icon(nc(i),i)=L
			ior(nc(i),i)=2
			else
			icon(nc(i),i)=L
			ior(nc(i),i)=3
			end if
	
	END IF
	END IF
	END IF
	END DO
	END IF
140 	CONTINUE
	RETURN
	END 
	
C	===================================================================
	SUBROUTINE PDBSIDE(XA,XS,SCM,NRES,ICODE,ibegin,a1,sb)
	PARAMETER(nmax=3000)
	LOGICAL QSIDE
	DIMENSION X(3),xa(3,4,NMAX),ICODE(-NMAX:NMAX),
     &	XS(3,30,NMAX),SCM(3,NMAX)
	character*1 c1,c1a,c1b,c2,chain,a1
	CHARACTER*4 NAMES,ITYPE,name1
	CHARACTER*6 ATOM
	CHARACTER*3 ARES,AA(0:19),bden(4)
	character*6 ad6
	character*7 ad7
	character*3 ad3
	dimension ic1(40),ic2(40)
	dimension ic1s(40),ic2s(40),sb(3,nmax)

	DIMENSION ISIDE(0:21),kside(NMAX),scma(3,nmax)
	CHARACTER*4 JDEN(0:19,10)
	COMMON qside(-3:nmax)
	DATA AA/ 'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &		     'MET','ASP','ASN','LEU',
     &		     'LYS','GLU','GLN','ARG',
     &		     'HIS','PHE','TYR','TRP'/
	DATA ISIDE/0,1, 2, 2, 3, 3, 4, 3,
     &                  4,4,4,4,
     &                  5,5,5,7,
     &			6,7,8,10,2,6/
	DATA (JDEN(0,J),J=1,10)/10*' '/
	DATA (JDEN(1,J),J=1,10)/'CB',9*' '/
	DATA (JDEN(2,J),J=1,10)/'CB','OG',8*' '/
	DATA (JDEN(3,J),J=1,10)/'CB','SG',8*' '/
	DATA (JDEN(4,J),J=1,10)/'CB','CG1','CG2',7*' '/
	DATA (JDEN(5,J),J=1,10)/'CB','OG1','CG2',7*' '/
	DATA (JDEN(6,J),J=1,10)/'CB','CG1','CG2','CD1',6*' '/
	DATA (JDEN(7,J),J=1,10)/'CB','CG','CD',7*' '/
	DATA (JDEN(8,J),J=1,10) /'CB','CG','SD','CE',6*' '/
	DATA (JDEN(9,J),J=1,10)/'CB','CG','OD1','OD2',6*' '/
	DATA (JDEN(10,J),J=1,10)/'CB','CG','OD1','ND2',6*' '/
	DATA (JDEN(11,J),J=1,10)/'CB','CG','CD1','CD2',6*' '/
	DATA (JDEN(12,J),J=1,10)/'CB','CG','CD','CE','NZ',5*' '/
	DATA (JDEN(13,J),J=1,10)/'CB','CG','CD','OE1','OE2',
     &	5*' '/
	DATA (JDEN(14,J),J=1,10)/'CB','CG','CD','OE1',
     &	'NE2',5*' '/
	DATA (JDEN(15,J),J=1,10)/'CB','CG','CD','NE',
     &	'CZ','NH1','NH2',3*' '/
	DATA (JDEN(16,J),J=1,10)/'CB','CG','ND1','CD2',
     &	'CE1','NE2',4*' '/
	DATA (JDEN(17,J),J=1,10)/'CB','CG','CD1','CD2',
     &	'CE1','CE2','CZ',3*' '/
	DATA (JDEN(18,J),J=1,10)/'CB','CG','CD1','CD2',
     &	'CE1','CE2','CZ','OH',2*' '/
	DATA(JDEN(19,J),J=1,10)/'CB','CG','CD1','CD2','NE1','CE2'
     &	,'CE3','CZ2','CZ3',
     &				'CH2'/
	data bden/' N', ' C',' CA',' O'/
	

	NBACK=0		
	DO I=1,NMAX
	QSIDE(I)=.TRUE.
	ICODE(I)=0
		do jk=1,4
			DO J=1,3
			XA(J,jk,I)=0.
			END DO
		END DO
	END DO
	ICOUNT=0
C	OPEN(UNIT=8,FILE=NAME//'side.pdb')
	NTEST=0
	NRES=0
	INUMOLD=-200
	ntest2=0
1	CONTINUE
	ntest=ntest+1
c	if(ntest .gt.nmax) go to 20
9	READ(10,10,ERR=1,END=20)ATOM,ID,ITYPE,
     &	ARES,CHAIN,INUM,X(1),X(2),X(3)
10	FORMAT(A6,I5,1X,A4,1x,a3,1X,a1,I4,1x,3X,3F8.3)
c 	IF(ATOM.EQ. 'END' .OR. ATOM .EQ. 'TER')THEN
	IF(ATOM.EQ. 'END' )THEN
		WRITE(6,*)'END OF FILE'
		GO TO 20
	END IF
	if(a1 .ne.'_')then
	if(CHAIN .ne. a1)then
	go to 1
	else
	if(atom .eq.'TER')go to 20
	if(atom .eq.'ENDMDL')go to 20
	end if
	else
	if(atom .eq.'TER')go to 20
	if(atom .eq.'ENDMDL')go to 20
	end if
	IF(ATOM .EQ. 'HETATM')GO TO 1
	if(ares .EQ. 'ACE')GO TO 1
	IF(ATOM .NE. 'ATOM  ')GO TO 1
	if(ntest2 .eq. 0)then
	ibegin=inum-1
	write(6,*)'ibegin=',inum,ibegin
	end if
	ntest2=ntest2+1
	IF(ITYPE .EQ. BDEN(3) .AND. INUM .NE. INUMOLD)THEN
C	BEGIN READ IN OF NEW RESIDUE
C		CHECK IF THE OLD SIDECHAIN IS COMPLETE
			IF(NRES .GT.0)THEN
			kside(nres)=icount
			IF(ICOUNT .NE. ISIDE(JRES))THEN
			QSIDE(NRES)=.FALSE.
			end if		
  		END IF
	nres=nres+1			
225	JCOUNT=0
	ICOUNT=0
	INUMOLD=INUM
			DO IRES=0,19
				IF(ARES .EQ. AA(IRES))THEN
				JRES=IRES
				ICODE(NRES)=IRES
				GO TO 15
			END IF
			END DO
15 	CONTINUE

		js=3
			IF(JRES .EQ. 0)THEN
				DO J=1,3
				XS(J,1,NRES)=X(J)
				END DO
			END IF
				DO J=1,3
				XA(J,js,NRES)=X(J)
				END DO

		NBACK=NBACK+1
	END IF
	IF(jres.ne. 0  .and. INUM .EQ. INUMOLD)then
	DO JS=1,ISIDE(JRES)
	IF(ITYPE .EQ. JDEN(JRES,JS) .
     &	OR. ITYPE .EQ. ' '//JDEN(JRES,JS))THEN
	ICOUNT=ICOUNT+1
		DO J=1,3
c	!!!  the xs coordinates are not relative but absolute
		XS(J,JS,NRES)=X(J)
			END DO
C	write(8,10)ATOM,ID,JDEN(jres,js),AA(jres),NRES,(Xs(J,js,nres),J=1,3)
		GO TO 22
	END IF
		END DO
	end if
22	CONTINUE
	GO TO 9
20	continue
c		IF(ntest .gt. nmax)then
c	write(6,*)ntest,'exceeds nmax ABORT!'
c	stop
c	end if

	if(icode(nres).ne.0)then
	kside(nres)=icount
	else
	kside(nres)=1
	end if
	do I=1,NRES
	IRES=ICODE(I)	
	NSIDE=ISIDE(IRES)
	if(nside .ne. kside(i))then

	qside(i)=.false.
	end if
	IF(qside(i))then
		IF(IRES .NE. 0)THEN
						do j=1,3
				SCM(j,I)=0.
				SCMA(j,I)=0.				
					do js=1,nside
					SCM(J,I)=SCM(J,I)+xs(j,js,i)
					end do
				scma(j,i)=scm(j,I)+xa(j,3,i)
				scma(j,i)=scma(j,i)/(kside(I)+1)
				scm(j,i)=scm(j,i)/kside(i)
				sb(j,i)=scm(j,i)-xa(j,3,i)
				end do
		ELSE
		kside(i)=1
			DO J=1,3
			SCM(J,I)=XA(J,3,I)
			SCMA(J,I)=XA(J,3,I)			
			END DO
	if(i .gt. 1 .and. i .lt. nres)Then
			do j=1,3
			SB(j,i)=2*xa(j,3,i)-xa(j,3,i-1)-xa(j,3,i+1)
			end do
	else if(i .eq.1)then
			do j=1,3
		SB(j,1)=-2*xa(j,3,2)+xa(j,3,1)+xa(j,3,3)
			end do
	else if(i .eq.nres)then
			Dov j=1,3
		SB(j,nres)=-2*xa(j,3,nres-1)+xa(j,3,nres-2)+xa(j,3,nres)
			end do
	end if		
		END IF
	END IF
	END DO
	close(10)
	RETURN
	END
C
c	##########################################################
	FUNCTION IR(X)
	DIF=X-INT(X)
	IF(DIF .LE. 0.5) THEN
	IDIF=0
	ELSE
	IDIF=1
	END IF
	IR=INT(X) +IDIF
	RETURN
	END
	function dot(a,b)
	dimension a(3), b(3)
	dot=0.
	do j=1,3
	dot=dot+a(j)*b(j)
	end do
	return
	end 
