	PARAMETER(ncases=3000)
	character*255 namep,namep3(ncases),namep2,namev(ncases)
	character*5 name5
	logical*1 touch(ncases)
	character*4 genome
	character*6 prodin
	character*20 a33
	character*255 home
	dimension x3(3)
c	&&&&&&& home:    &&&&&&&&&&
	open(unit=1,file='homedir')
	read(unit=1,fmt='A255')home
	nhome=0
	do i=255,1,-1
	if(home(i:i) .eq. ' ')nhome=i
	end do
	nhome=nhome-1
c	&&&&&&& home:    &&&&&&&&&&


	open(unit=55,file='nametarg1')
        read(55,13)prodin
13      format(a6)
	open(unit=9,file='stat_bad')
	open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)
	OPEN(unit=4,file='LIST.targ'//prodin)
	read(4,*)nd,genome

	open(unit=22,file='listbest')
	read(22,*)nfz
	do ik=1,nfz
	read(22,11)namep3(ik)(1:nsize)
11	format(a,1x,a15)
	end do
	nn=nfz
	open(unit=22,file='list_verygood')
	read(22,*)nfz
	do ik=1,nfz
	nn=nn+1
	read(22,11)namep3(nn)(1:nsize)
	end do
	open(unit=22,file='list_fair')
	read(22,*)nfz
	do ik=1,nfz
	nn=nn+1
	read(22,11)namep3(nn)(1:nsize)
	end do

	write(9,*)'best, verygood & fair cases ',nn
	open(unit=23,file='listbad')
          open(unit=2,file='STRUCTURE_rap3orienrev')

	read(2,*)ndata2
	ntarg=0
	do 1001 ik=1,ndata2
	read(2,1)namep(1:nsize)
1	format(9x,a)
        OPEN(unit=3,file=namep(1:nsize)//'.threadrap3orienrev')
	read(3,*)nth
        read(3,18)name5,zth
		do jk=1,nn
		if(namep3(jk)(1:nsize).eq.namep(1:nsize))go to 1001
		end do

	nn=nn+1
	namep3(nn)=namep(1:nsize)
	OPEN(unit=59,file=home(1:nhome)//genome//'zpot4a/'//namep(1:nsize)//'.predictedrap3orienrev')
        open(unit=69,file=home(1:nhome)//genome//'qzmediumcontacts/'//namep(1:nsize)//'.contacts')
        read(59,*)nt	
	write(69,*)nt
	do jd=1,nt
	read(59,*)i,j
	write(69,*)i,j
        end do
	close(59)
	close(69)
	 OPEN(unit=34,file=home(1:nhome)//genome//'pdborienrev/'//namep(1:nsize)//'rap3orienrev1.pdb')

        OPEN(unit=35,file=home(1:nhome)//genome//'qzmediumtemplates/'//namep(1:nsize)//'.pdb')
	open(unit=31,file=home(1:nhome)//genome//'qzmediumsummary/'//namep(1:nsize)//'.thread')

	rewind(3)
	write(31,*)nth
	do mk=1,nth
        read(3,18)name5,zth
	write(31,18)name5,zth
	read(34,15)name5,ic
	write(35,15)name5,ic
	do id=1,ic
        read(34,1033)a33,i,(x3(l),l=1,3),j
        write(35,1033)a33,i,(x3(l),l=1,3),j	
	end do
	read(34,*)
        write(35,906)'TER'
c	==============================================================
88	continue
18      format(a5,1x,1f7.3)
1033    FORMAT(a20,2X,I4,1x,3X,3F8.3,1x,i4,1x,a3)
906     format(a3)
15	format(a5,1x,i4,1x,1f7.3,1x,i3)
c

	END DO
c	==============================================================
	ntarg=ntarg+1
	namev(ntarg)(1:nsize)=namep(1:nsize)
1001	continue
	write(23,*)ntarg
	do ik=1,ntarg
	write(23,63)namev(ik)(1:nsize)
63	format(a)
	end do

	close(2)
	open(unit=1,file='finish_bad')
	rewind(1)
	write(1,*)'1'
	close(1)

	stop
	end
