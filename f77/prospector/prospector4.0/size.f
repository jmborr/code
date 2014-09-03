	PARAMETER(NCASES=6000)
	character*255 f3590,f3590aln,namep,namep2(ncases)
	dimension nres(ncases)

c	readin  the directories"
c	&&&&&&& ENTER 3590 set:    &&&&&&&&&&
	open(unit=1,file='3590')
	read(unit=1,fmt='A255')f3590
	n3590=0
	do i=255,1,-1
	if(f3590(i:i) .eq. ' ')n3590=i
	end do
	n3590=n3590-1
	open(unit=1,file='3590aln')
	read(unit=1,fmt='A255')f3590aln
	naln90=0
	do i=255,1,-1
	if(f3590aln(i:i) .eq. ' ')naln90=i
	end do
	naln90=naln90-1
	open(unit=1,file='LISTall')
	read(1,*)ndata
	open(unit=11,file='namesize')
	read(11,*)nsize
	close(11)
	nn=0
	do ik=1,ndata
	read(1,1)namep(1:nsize)
1	format(a)
	OPEN(UNIT=31,FILE=f3590(1:n3590)//namep(1:nsize)//f3590aln(1:naln90))
	read(31,*)nd,mres
	if(mres .lt.501)then
	nn=nn+1
	nres(nn)=mres
	namep2(nn)(1:nsize)=namep(1:nsize)
	end if
	end do
	close(1)
	open(unit=1,file='LIST')
	write(1,*)nn
	do ik=1,nn
	write(1,2)namep2(ik)(1:nsize),nres(ik)
2	format(a,1x,i5)
	end do
	close(1)
	stop
	end 
