	program check
c this program is to check the stability of the Bessel
c function algorithm
c kwh feb 20 2010
c ok this bessel function routine changed to check for the case where
c you might just bounce largely positive because of loss of accuracy
c as the Bessel function absolute value dips below 10^-8

	real*8 bj(0:25),qa,dlim
	real*8 qamin,dqa
	integer i,j,iflag

100	continue
	write(6,*)' enter dlim: '
	read(5,*)dlim
	iflag=0
	qamin=0.1
	dqa=0.0000001
	write(6,*)'qa max:',dqa*200000001
	do 200 i=1,200000001
	  qa=qamin+float(i-1)*dqa
	  call sphbes(bj,qa,dlim)
	  do 190 j=0,25
		if (bj(j) .gt. 1.) then 
	           iflag=1
		end if
190	  continue
	  if (iflag .eq. 1) then
	  write(6,*)' qa:',qa,' dlim:',dlim
	  do 191 j=0,25
	     write(6,*) j,bj(j)
191	  continue
	  iflag=0
	  endif
200	continue

	go to 100
	stop
	end
	

      subroutine sphbes(bj,qa,dlim)
      implicit real*8(a-h,o-z)
	integer flag
      real*8 bj(0:25)
      real*8 qa
      bj(0) = sin(qa)/qa
      bj(1) = sin(qa)/(qa*qa) - cos(qa)/qa
	flag=1
      do 10 i=2,25
	if (dabs(bj(i-1)) .le. dlim) then 
		flag=0
	endif
	if (flag .eq. 1) then
        bj(i) = float(2*i-1) * bj(i-1) / qa - bj(i-2)
	else
		bj(i)=dlim
	endif
10    continue
        
      end
