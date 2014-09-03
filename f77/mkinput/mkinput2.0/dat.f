c       compiling line
c       pgf77 -Mextend -O -s -fast -Wl,-static -o dat.x dat.f f2kcli.f
c       calculate par.dat (pair pot) comb.dat (contact restr ) for SG
c       dist.dat distL.dat (dist restr)  combCA.dat (contact restr ) for CA
c       ./dat.x ntemp,irest xxxxxrap3orienrev5s.pdb
c       where ntemp is number of templates in xxxxxrap3orienrev5s.pdb
c       and irest is target sequence length. xxxxxrap3orienrev5s.pdb is
c       output from prospector (list of templates):
c
c       1ak6_  138  28.111   4
c       ATOM     1   CA  MET     1       3.382  19.011  11.345    7 SER
c       ...
c       ATOM     1   CA  ASN   139      41.211  21.608   2.514  162 GLY
c       TER   0.348
c       1ahq_  132  25.765   4  
c       ...
c
c       where 20 means number of templates, and 150 is sequence length.
c       Next there are 20 lines, each providing info for a template.
c       
	program dat

	implicit real*8 (a-h,o-z)
        parameter (maxres=600,maxmd=200,iwin=9,maxp=20000,
     &  maxd=maxres*maxmd,maxdp=20000000,maxq=600)
        common/input/ntemp,irest,tempin,ltempin !command line arguments
        character whole*80,ctmp6*6,frgfil*255,exclfil*255,
     &  resn(20)*3,mypdb*10,rname*3,line*80,mytemp(maxp)*6
     &  ,exclist(maxp)*6,datbase*255,atmn*3,resnyz(20)*3
     &  ,tempin*255,temppdb(maxp)*6,whole1(maxres)*80
        integer ifrg(maxres,2),itemp(maxres,maxmd,2),myseq(maxres),
     &  tempid(maxp),tempres(maxp),tempnum(maxres),templst(maxres,maxmd)
     &  ,tempused(maxp),temppos(maxp),tempmap(maxd,2),excind(maxp),
     &  excnum,pdbselected(maxp),isort(maxp),ibase(maxdp),frgpos(maxdp),
     &  invmap(maxres),resmap(20),ndist(maxq,maxq),rnum,ltempin,
     &  ndista(maxq,maxq),tempseq(maxres),exclid(maxp)
        real*8  tempscore(maxres,maxmd),tempxyz(maxres,maxmd,iwin,3),
     &  cacnt(maxres,maxres),zscore(maxp),cacnt0(maxres,maxres)
     &  ,sgcut(20,20),cacut(20,20),sgcut_dev(20,20),cacut_dev(20,20)
     &  ,frgsc(maxdp),scmx(maxres,maxmd),scmy(maxres,maxmd),
     &   scmz(maxres,maxmd),cax(maxres,maxmd),cay(maxres,maxmd),
     &   caz(maxres,maxmd),cx,cy,cz,freq_comb(maxq,maxq),
     &   freq0_par(maxres,maxres),dist(maxq,maxq,maxmd),
     &   dista(maxq,maxq,maxmd),freq_combCA(maxq,maxq),ee(maxq,maxq)
     &   ,disave(maxp),disdelta(maxp),cax1(0:maxres),cay1(0:maxres),
     &    caz1(0:maxres),scmx1(maxres),scmy1(maxres),scmz1(maxres)
     &   ,pij(maxres,maxres)
         integer irestyp(maxres),jtemp(maxres), 
     &   mee(maxq,maxq),iid(maxp),jjd(maxp)
        data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/
        data resnyz/'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     &  'MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE',
     &  'TYR','TRP'/

	call process_command_line() !read input arguments
c       unpack common files
	call system('tar jxf /gpfs1/active/jose/code/f77/tasser/tasser.1.0/common.tbz2 -C .')
c	contact.comm contains average distances between side-groups
	open(unit=10,file='contact.comm',status='old',err=1000)
	read(10,*) 
	do i=1,20
	   read(10,*) rname,(sgcut(i,j),j=1,20)                     
	enddo
c       contact.comm contains standard deviation around the averages
	read(10,*) 
	do i=1,20
	   read(10,*) rname,(sgcut_dev(i,j),j=1,20)                     
	enddo
c       generous cut-off
	do i=1,20
	   do j=1,20
	      sgcut(i,j)=sgcut(i,j)+sgcut_dev(i,j)*2.5
	   enddo
	enddo
c       contact.comm contains average distances between CA's
	read(10,*) 
	do i=1,20
	   read(10,*) rname,(cacut(i,j),j=1,20)                     
	enddo
c       contact.comm contains standard deviation around the averages
	read(10,*) 
	do i=1,20
	   read(10,*) rname,(cacut_dev(i,j),j=1,20)                     
	enddo
c       generous cut-off
	do i=1,20
	   do j=1,20
	      cacut(i,j)=cacut(i,j)+cacut_dev(i,j)*2.5
	   enddo
	enddo
	close(10)
c       initialize several variables
	do i=1,maxq
	   do j=1,maxq
	      ndist(i,j)=0
	      ndista(i,j)=0
	      freq_comb(i,j)=0.0
	      freq0_par(i,j)=0.0
	   enddo
	enddo
	pcut=0.3 !a cut-off probability
c       read number of templates and target sequence length
c	read in template  coordinates:
	do i=1,maxres
	   jtemp(i)=0
	enddo
c       open file containing the templates
	open(unit=10,file=tempin,status='old',err=1010)
	ibs=index(tempin,' ')-1 !number of characters in template file name
c       process each template
	do im=1,ntemp
	  do i=1,irest
	     frgpos(i)=-1 !initialize
	  enddo
	  id=1
	  ires0=-100 !initialize ires0
	  rnum=0
c       read up to the first ATOM CA entry. whole1(id) contains all the
c       ATOM CA entries. One entry is of the form
c
c       ATOM     1   CA  THR     1      21.901  21.313  12.081   40 ALA
c
c       40 ALA indicates the coordinates correspond to ALA residue at
c       sequence position 40 in the template sequence
15	  read(10,'(a80)',end=10) whole1(id)
	  if(whole1(id)(1:3).eq.'TER') goto 10
	  if(whole1(id)(1:3).eq.'END') goto 10
	  if(whole1(id)(1:4).ne.'ATOM') goto 15
	  if(whole1(id)(14:15).ne.'CA') goto 15	! CA only
	  
	  read(whole1(id)(23:26),'(i4)') ires1 !residue number
	  if(ires1.ne.ires0) then
             if(ires0.lt.0) ist=ires1 !starting position read
             ires0=ires1 !ires0 is residue number of residue just read
             rname=whole1(id)(18:20) !3-letter code for residue
             do ll=1,20
c       tempseq(ires1) is numeric code for the amino acid at position ires1
		if(resnyz(ll).eq.rname) tempseq(ires1)=ll
c       irestyp(ires1) is another numeric code
		if(resn(ll).eq.rname) irestyp(ires1)=ll
             enddo
c       frgpos(ires1) also indicates that the template is covering this
c       position
	     frgpos(ires1)=ires1
c       jtemp(ires1) is number templates that cover position ires1
	     jtemp(ires1)=jtemp(ires1)+1
	  endif             

	  read(whole1(id)(31:54),'(3F8.3)') xx,yy,zz
	  cax1(ires1)=xx
	  cay1(ires1)=yy
	  caz1(ires1)=zz
	  id=id+1
	  goto 15 !read next ATOM CA entry
 10	  iend=ires0 !finished reading template, mark ending position read

c       estimate position of side-groups based on CA positions. Recall that
c       irest is target sequence length, ist is the starting
c       amino acid position covered by the template, iend is the ending
c       amino acid position covered by the template, irestyp is numeric
c       code for amino acid type.
	  call getsg(irest,ist,iend,irestyp,frgpos,cax1,cay1,caz1,scmx1,
     &    scmy1,scmz1)          
	 
c       freq_comb(l,j) number of times l and j side-groups below
c       distance cut-off. In order to normalize and give a probability,
c       we divide by freq0_par(l,j)
	  do l=1,irest-1
	     l1=frgpos(l)
	     do j=l+1,irest
		j1=frgpos(j)
		if(l1.gt.0.and.j1.gt.0) then
		   dis=sqrt((scmx1(l)-scmx1(j))**2
     &                      +(scmy1(l)-scmy1(j))**2
     &                      +(scmz1(l)-scmz1(j))**2)
		   freq0_par(l,j)=freq0_par(l,j)+1
		   if(dis.lt.sgcut(tempseq(l1),tempseq(j1))) then
		      freq_comb(l,j)=freq_comb(l,j)+1
		   endif		   
		endif
	     enddo
	  enddo

	  do l=1,irest-1
	     l1=frgpos(l)
	     do j=l+1,irest
		j1=frgpos(j)
		if(l1.gt.0.and.j1.gt.0) then
		   dis=sqrt((cax1(l)-cax1(j))**2
     &                      +(cay1(l)-cay1(j))**2
     &                      +(caz1(l)-caz1(j))**2)
c       dista(l,j, ndista(l,j)) lists all distances between l and j. The
c       number of such distances is ndista(l,j)
		   ndista(l,j)=ndista(l,j)+1
		   dista(l,j, ndista(l,j))=dis
c       short range distance (sequence separation below 6):
		   if((j-l).le.6) then
		      do iu=l,j
c       there's an internal gap, so we don't know what is the real
c       sequence separation
			 if(frgpos(iu).lt.0) goto 301 
		      enddo
c       we require consecutive CA atoms below 4.1Angstrom distance. The
c       reason is that we may not have gaps in the sequence, but there
c       may be real gaps because the PDB file for the template contained
c       mistakes as hidden gaps.
		      do iu=l+1,j
			 dis2=sqrt((cax1(iu-1)-cax1(iu))**2
     &                      +(cay1(iu-1)-cay1(iu))**2
     &                      +(caz1(iu-1)-caz1(iu))**2)
			 if(dis2.gt.4.1) goto 301
		      enddo
c       dist(l,j, ndist(l,j)) lists all distances between l ad j such
c       that j-l<6, ie, local contacts
		      ndist(l,j)=ndist(l,j)+1
		      dist(l,j, ndist(l,j))=dis
		   endif
 301		   continue
		endif
	     enddo!matches do j=l+1,irest
	  enddo !matches do l=1,irest-1
	enddo !matches 	do im=1,ntemp
 	close(10) !we finished reading the file containing the templates

	do i=1,irest
	   jtemp(i)=ntemp !assume residue i covered by all ntemp templates
	enddo 

c       output comb.dat for SG
	nres=0			!number of restraints
	do i=1,irest-5		!these are long-range restraints (j-i>6)
	   do j=i+5,irest
	      pij(i,j)=freq_comb(i,j)/sqrt(float(jtemp(i)*jtemp(j)))
c       recall pcut=0.3
              if(pij(i,j).gt.pcut) nres=nres+1
	   enddo
	enddo
	open(unit=20,file='comb.dat',status='replace')
	write(20,*) nres
	do i=1,irest-5
	   do j=i+5,irest
              if(pij(i,j).gt.pcut) write(20,130) i,j,pij(i,j)
	   enddo
	enddo
 130	format(5x,i5,3x,i5,3x,f5.3)
	close(20)

c       output par.dat for SG:
	ee_a=0 !summ of all ee(i,j)
	n_ee=0 !number of (i,j) pairs
        do i=1,irest-1
           do j=i+1,irest
	      mee(i,j)=-1
	      if(freq0_par(i,j).gt.0) then
		 frequence=freq_comb(i,j)/(freq0_par(i,j)
     &           /(0.8*sqrt(float(jtemp(i)*jtemp(j)))))
		 if(frequence.gt.0.0001) then
c       logarithm of the expected frequencies for the (i,j) pair
		    ee(i,j)=-log(frequence)
		 else
		    ee(i,j)=0.0
		 endif
		 n_ee=n_ee+1
		 mee(i,j)=1 !Marks there is a (i,j) pair
		 ee_a=ee_a+ee(i,j)
	      endif
           enddo !matches do j=i+1,irest
        enddo !matches do i=1,irest-1

c       shift ee with respect to ee_a
	ee_a=ee_a/n_ee
	do i=1,irest-1
	   do j=i+1,irest
              if(mee(i,j).gt.0) then
		 ee(i,j)=ee(i,j)-ee_a
              else
		 ee(i,j)=0.0
              endif
              ee(j,i)=ee(i,j)
	   enddo
	enddo

	open(unit=20,file='par.dat',status='replace')
	do i=1,irest
           write(20,'(5x,i5,1x,15a)') i,"==============="
           write(20,'(1x, 10f8.3)')(ee(i,j),j=1,irest)
	enddo
	close(20)
           
c       output combCA.dat for CA
	nres=0 !number of restraints
	do i=1,irest-5
	   do j=i+5,irest
	      nn=ndista(i,j)
	      freq_combCA(i,j)=0
	      do k=1,nn !summ all contacts with distances less than 6Angstroms
                 if(dista(i,j,k).lt.6.0) then
		    freq_combCA(i,j)=freq_combCA(i,j)+1
                 endif
	      enddo
	      pij(i,j)=freq_combCA(i,j)/sqrt(float(jtemp(i)*jtemp(j)))
	      if(pij(i,j).gt.pcut) nres=nres+1
	   enddo
	enddo
	  
	open(unit=20,file='combCA.dat',status='replace')
	write(20,*) nres
	do i=1,irest-5
	   do j=i+5,irest
              if(pij(i,j).gt.pcut) write(20,130) i,j,pij(i,j)
	   enddo
	enddo   
	close(20)

c       output dist.dat  for CA
	nres=0
	xdist0=2
	if(ntemp.lt.2) xdist0=1
	do i=1,irest-1
	   jm=i+6
	   if(jm.gt.irest) jm=irest       
	   do j=i+2,jm	      
              kk0=ndist(i,j) !number contacts between i and j for all templates
              if(kk0.ge.xdist0) then
		 dist_a=0
		 dist2_a=0
                 do k=1,kk0
		    dist_a=dist_a+dist(i,j,k) !average distance
		    dist2_a=dist2_a+dist(i,j,k)**2 !average of the square
                 enddo
		 dist_a=dist_a/kk0
		 dist2_a=dist2_a/kk0
		 delta2=dist2_a-dist_a**2
		 if(delta2.lt.0.00001) then
		    delta=0.0
		 else          
		    delta=sqrt(delta2) !standard dev around average distance
		 endif
		 nres=nres+1
		 iid(nres)=i                  
		 jjd(nres)=j                  
		 disave(nres)=dist_a !store average distance
		 disdelta(nres)=delta!store standar deviation
              endif
	   enddo
	enddo

	open(unit=20,file='dist.dat',status='replace')
	write(20,'(5x,i5)') nres
	do i=1,nres
           write(20,'(1x,i5,1x,i5,1x,i5,1x,f8.3,1x,f8.3)') 
     &     iid(i),jjd(i),ndist(iid(i),jjd(i)),disave(i),disdelta(i)
	enddo
	close(20)           

c       output distL.dat  for CA
	nres=0			!number of distance restraints
	do i=1,irest-10		!output only if sequence separation above 10 res.
	   do j=i+10,irest,10
	      nn=ndista(i,j)
	      do k=1,nn
		 nres=nres+1
	      enddo
	   enddo
	enddo

	open(unit=20,file='distL.dat',status='replace')
	write(20,*) nres
	do i=1,irest-10
	   do j=i+10,irest,10
	      nn=ndista(i,j)
	      do k=1,nn
		 write(20,140) i,j,dista(i,j,k)
	      enddo
	   enddo
	enddo   
 140	format(5x,i5,3x,i5,3x,f7.3)
	close(20)
	call cleanup() !remove created temporary files
	call exit(0) !exit with success
c       error messages and exit with error exit code
 1000	write(*,*) 'ERROR in dat.x: file contact.comm not found'
	goto 2000
 1010	write(*,*) 'ERROR in dat.x: file ',tempin(1:ltempin),' not found'
 2000	call cleanup()
	call exit(1)
	end !end of program dat.f
          
c       ###########################################
          
	subroutine cleanup()
	call system('/bin/rm *.comm')
	end

c       ###########################################

	subroutine process_command_line()
c       variables related to command line arguments
        common/input/ntemp,irest,tempin,ltempin
	integer ntemp,irest,ltempin
	character tempin*255
c       variables intrinsic to the subroutine
	integer ireq,required,narg,iarg,length,istatus
	character onearg*256,option

	include 'f2kcli.inc' !it is neccessary to include this header file

	required=3 !number of required arguments
	narg=command_argument_count()
	ireq=0
	iarg=0
 5	continue
	iarg=iarg+1
	call get_command_argument(iarg,onearg,length,istatus)
	option=onearg(2:2)
c       number of templates
	if(option.eq.'a')then
	   ireq=ireq+1
	   iarg=iarg+1
	   call get_command_argument(iarg,onearg,length,istatus)
	   read(onearg,*) ntemp
c       number of residues
	else if(option.eq.'b')then
	   ireq=ireq+1
	   iarg=iarg+1
	   call get_command_argument(iarg,onearg,length,istatus)
	   read(onearg,*) irest
c       file containing the templates
	else if(option.eq.'c')then
	   ireq=ireq+1
	   iarg=iarg+1
	   call get_command_argument(iarg,tempin,length,istatus)
	   ltempin=length
c	asking for help
	else if(option.eq.'h')then
	   call message()
c       wrong option
	else
	   write(*,*)'Wrong option passed'
	   call message()
	endif
c       check if we finished processing the arguments
	if(iarg.lt.narg) then
	   goto 5
	endif
	if(ireq.lt.required)then !not enough number of required arguments
	   call message()
	endif
	end

c	###########################################

	subroutine message()
	write(*,*)'Usage: ./dat.x -a ntemp -b irest -c tempin'
	write(*,*)'required arguments:'
	write(*,*)' ntemp: number of templates in tempin to be read'
	write(*,*)' irest: number of residues for query sequence'
	write(*,*)' tempin: file containing templates. Format is same as prospector output'
	call exit(1) !exit with error
	end

c	###########################################

	subroutine setcacnt(mypdb,cacnt,zs_i)
        parameter(maxres=600)
	real*8 cacnt(maxres,maxres),zs_i
        character mypdb*6
           
	ipos=index(mypdb,' ')-1
	if(ipos.le.0) ipos=6
	open(unit=20,file=mypdb(1:ipos)//'_sp3temp.cnt',
     &  status='old')

 10	read(20,*,end=50) i1,i2,ca
	cacnt(i1,i2)=zs_i
	cacnt(i2,i1)=zs_i
	goto 10
 50	close(20)
	return
	end
c     ########################################
	subroutine readexc(excfil,exclist,excnum,excind)
	character excfil*255,exclist(*)*6
	integer excnum,excind(*)

	open(unit=70,file=excfil,status='old')
	excnum=1
 40	read(70,'(a6)',end=50) exclist(excnum)
	idd=index(exclist(excnum),' ')-1
	if(idd.le.0) idd=6
	excind(excnum)=idd
	excnum=excnum+1
	goto 40
 50	close(70)
	excnum=excnum-1
	
	return
        end

***********************************************************************
        subroutine sortzscore(npro,zscore,zscore_sort)
        real*8 zscore(*)
        integer zscore_sort(*)
c       use unix command sort

	open(unit=50,file='_sort.in',status='replace')
	do i=1,npro
	   write(50,*) i,zscore(i)
	enddo
	close(50)
	call system('sort -nr +1 _sort.in > _sort.out')
	open(unit=50,file='_sort.out',status='old')
	do i=1,npro
	   read(50,*) zscore_sort(i)
	enddo
	close(50)
	
	return
        end
c       #########################
	subroutine topnzscore(npro,zscore,isort,nn)
        parameter (maxres=600,maxp=20000)
	real*8 zscore(*),xx,yy
	integer isort(*),idel(maxp)
	do i=1,npro
	   idel(i)=i
	enddo
	
	do i=1,nn
	   isort(i)=i
	   xx=zscore(i)
	   do j=1,npro
	      if(idel(j).gt.0) then
		 yy=zscore(j)
                 if(yy.gt.xx) then 
                    j0=j
                    xx=yy
                 endif   
	      endif
	   enddo  
	   isort(i)=j0
	   idel(j0)=-1
	enddo
	return
	end

c###########################################

c       estimate position of side-groups based on CA positions. Recall that
c       rnum is sequence lenght, ist is the starting
c       amino acid position covered by the template, iend is the ending
c       amino acid position covered by the template, irestyp is numeric
c       code for amino acid type.
        subroutine getsg(rnum,ist,iend,irestyp,frgpos,cax,cay,caz,scmx,
     &  scmy,scmz)
	parameter(maxa=100000,ntyp=20,maxres=600,nvec=312,
     &	maxdp=20000000)
	integer i,j,k,jm,jp,rnum, seq(maxres),frgpos(maxdp),irestyp(*)
	real*8 cax(0:maxres),cay(0:maxres),caz(0:maxres),
     &  scmx(*),scmy(*),scmz(*)
	character base*250,aa(-1:20)*3,resn(20)*3
	data init/0/
	data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
     &          'HIS','ARG','LYS','PRO'/

	data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
c                -1    0     1     2     3     4     5     6
     &     'PRO','MET','ASP','ASN','LEU',
c            7     8     9    10    11
     &     'LYS','GLU','GLN','ARG',
c           12    13    14    15
     &     'HIS','PHE','TYR','TRP','CYX'/
c           16    17    18    19    20

	common/beta3/seq 	
	common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
	common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)

	if(init.eq.0) then
	   init=1
	   open(unit=4,file='sidechain.comm',status='old') !for Sc position.
c       sidechain.comm contains estimated relative positions for side-groups
c       sidechain.comm format goes like this:
c
c       GLY
c       0.000  0.000  0.000    0.000  0.000  0.000
c       ALA
c       0.249 -1.118  0.975    0.113 -0.736  1.294
c       ...    
c       TRP
c       0.476 -1.156  2.541   -0.058 -0.427  2.894
c
c       first three numbers are components of relative coordinates
c       assuming an alpha-helix configuration in the vector base spanned
c       by a orthonormal base. Other three
c       numbers are for other configurations.
	   do k=0,19
	      read(4,*)
	      read(4,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
	   enddo
	   close(4) !sidecent.comm
        endif	    ! init

	do i=ist,iend
	   do l=0,19
	      if(resn(irestyp(i)).eq.aa(l)) seq(i)=l
	   enddo 
	enddo
		
	do i=ist,iend
	   if(frgpos(i).lt.0) goto 10
           j=i
c       jm:first residue to the left of j that is covered by the
c       template. jp is similar but to the right of j
	   call findneighbors(j,jm,jp,irest,frgpos)
c	   write(*,'3(i3)')j,jm,jp
           amx=cax(j)-cax(jm)
           amy=cay(j)-cay(jm)
           amz=caz(j)-caz(jm)
           aaa=sqrt(amx**2+amy**2+amz**2)
c       director cosines from j-1 to j
           amx=amx/aaa      
           amy=amy/aaa
           amz=amz/aaa
           apx=cax(jp)-cax(j)
           apy=cay(jp)-cay(j)
           apz=caz(jp)-caz(j)
           aaa=sqrt(apx**2+apy**2+apz**2)
c       director cosines from j to j+1
           apx=apx/aaa      
           apy=apy/aaa
           apz=apz/aaa
c       (i-1)--(i)--(i+1) angle
           ang=acos(-(amx*apx+amy*apy+amz*apz))*180/3.1415926 

	   aaax=amx+apx
	   aaay=amy+apy
	   aaaz=amz+apz
	   aaa=sqrt(aaax**2+aaay**2+aaaz**2)
c       director cosines for the vector sum of (amx,amy,amz) and (apx,apy,apz)
	   ax=aaax/aaa
	   ay=aaay/aaa
	   az=aaaz/aaa

	   ccx=amx-apx
	   ccy=amy-apy
	   ccz=amz-apz
	   aaa=sqrt(ccx**2+ccy**2+ccz**2)
c       director cosines for the bisector of (amx,amy,amz) and (apx,apy,apz)
	   cx=ccx/aaa
	   cy=ccy/aaa
	   cz=ccz/aaa

	   bbx=amy*apz-amz*apy
	   bby=amz*apx-amx*apz
	   bbz=amx*apy-amy*apx
	   aaa=sqrt(bbx**2+bby**2+bbz**2)
c       director cosines for vector perpend to (amx,amy,amz) and (apx,apy,apz)
	   bx=bbx/aaa
	   by=bby/aaa
	   bz=bbz/aaa
c       (ax,ay,az), (bx,by,bz), (cx,cy,cz) form an orthonormal base
	   k=seq(i)
	   if(ang.lt.105) then	! alpha-helix or turn like
              dx=(axalf(k)*ax+ayalf(k)*bx+azalf(k)*cx)
              dy=(axalf(k)*ay+ayalf(k)*by+azalf(k)*cy)
              dz=(axalf(k)*az+ayalf(k)*bz+azalf(k)*cz)
	   else			! beta-sheet
              dx=(axbet(k)*ax+aybet(k)*bx+azbet(k)*cx)
              dy=(axbet(k)*ay+aybet(k)*by+azbet(k)*cy)
              dz=(axbet(k)*az+aybet(k)*bz+azbet(k)*cz)
	   endif

           scmx(i)=cax(j)+dx 
           scmy(i)=cay(j)+dy
           scmz(i)=caz(j)+dz 
	   
 10	enddo !go to following residue
	end !end of subroutine getsg

c       ###############################################

	subroutine findneighbors(j,jm,jp,irest,frgpos)
        parameter (maxres=600,maxmd=200,iwin=9,maxp=20000,
     &  maxd=maxres*maxmd,maxdp=20000000,maxq=600)
	integer j,jm,jp,irest,frgpos(maxdp)
c	find first residue to the left that is covered
	jm=j-1
	do 10 while(jm.gt.0 .and. frgpos(jm).lt.0)
	   jm=jm-1
 10	continue
c	find first residue to the right that is covered
	jp=j+1
	do 20 while(jp.le.irest .and. frgpos(jp).lt.0)
	   jp=jp+1
 20	continue
c       handling special case jm.eq.0
	if(jm.eq.0) then
	   jm=jp
	   jp=jp+1
	   do 30 while(jp.le.irest .and. frgpos(jp).lt.0)
	      jp=jp+1
 30	   continue   
	endif
c       handling special case jp.eq.irest+1
	if(jp.eq.irest+1) then
	   jp=jm
	   jm=jm-1
	   do 40 while(jm.gt.0 .and. frgpos(jm).lt.0)
	      jm=jm-1
 40	   continue
	endif

	end

c       ###############################################
