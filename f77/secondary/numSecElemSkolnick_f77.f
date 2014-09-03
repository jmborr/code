c      program junk
c      parameter(nmax=2000)
c      dimension xa(3,nmax)
c      character*256 name
c      name='256bA.pdb'
c      write(*,*)'nsec=',nsec2(name)
c      end

ccccc ===============================================
c     same as nsec, but with a different argument
c     ===============================================
      function nsec2(name)
      parameter(nmax=2000)
      dimension xa(3,nmax)
      character*256 name
      nres=0
      call import_ca_coords(name,xa,nres)
      nsec2=nsec(nres,xa)
      return 
      end

ccccc ===============================================
c     compute number of secondary structure elements
c     reference(s):
c     a. kolinski, j. skolnick, a. godzik and w-p hu.  a method
c     for the prediction of surface u-turns and transglobular
c     connections in small proteins.  proteins 1997:27: 290-308.
c     w-p. hu, a. kolinski and j. skolnick.  improved method for
c     the prediction of the protein backbone u-turn positions
c     and the major secondary structures between the u-turns.
c     proteins 1997:29: 443-460
c     ===============================================
      function nsec(nres,xa)
      parameter(nmax=2000)
      integer sed(-2:nmax) 
      dimension awt(0:5)
      dimension xa(3,nmax),xb(3,-10:nmax),rc(3)
      dimension y(3,-7:nmax),rad(nmax)
      dimension mbegin(0:nmax)
      dimension mend(0:nmax)
      dimension msbegin(0:nmax)
      dimension msend(0:nmax)
      dimension nsectype(0:nmax)
      dimension msec(0:100)
      dimension minS(1:2)
      dimension ntype(1:2)
c     weights for smoothing CA coordinates
      awt(5)=1./243.
      awt(4)=5./243.
      awt(3)=15./243.
      awt(2)=30./243.
      awt(1)=45./243.
      awt(0)=51./243

ccccc minimal lengths and code types
      minS(1)=6 !helix
      minS(2)=4 !strand
      minl=4    !set to minimum of the two previous
      ntype(1)=2
      ntype(2)=4

ccccc create average tube:

c     store coordinates in xb, which is an extended version of xa,
c     since it allows insertion of extra residues at the beginning
c     and at the end
      do i=1,nres
         do j=1,3
            xb(j,i)=xa(j,i)
         end do
      end do
      
c     add virtual residues at beginning
      do i=-4,0
         do j=1,3
            xb(j,i)=xb(j,1)
         enddo
      end do

c     add virtual residues at end
      do i=1,5
         do j=1,3
            xb(j,i+nres)=xb(j,nres)
         end do
      end do

c     create smoothed backbone, by averaging over coordinates
c     of preceding and succeeding five residues
      do i=1,nres
         do j=1,3
            y(j,i)=0.
            do iw=-5,5
c     iwa is associated weigth
               iwa=abs(iw) 
               y(j,i)=y(j,i)+awt(iwa)*xb(j,i+iw)
            end do
         end do
      end do

c     substitute backbone with smoothed chain
      do i=1,nres
         do j=1,3
            xa(j,i)=y(j,i)
         end do
      end do

ccccc calculate radius of curvature passing through
c     each residue, store in rad()
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
c     sed() is signal for loop presence (1 for loop presence)
c     1/0.11=9angstroms
         if(sum .lt. 0.11)then
            sed(i)=0
         else
            sed(i)=1
         end if
 22      format(i4,3f6.4)
         write(*,'a20,i3,x,f5.3')'radius-of-curvature=',i,rad(i)
      end do
      icount=0
 100  continue

ccccc calculate number of loops (nloop), and residues
c     where each loop begins mbegin() and ends mend()
c     No loop can be defined in the first three and the
c     last three residues in the chain
      sed(0)=0
      sed(1)=0
      sed(2)=0
      sed(3)=0
      sed(nres-2)=0
      sed(nres-1)=0
      sed(nres)=0
      nloop=0
      do i=1,nres
         if(sed(i).eq.1 .and. sed(i-1) .eq. 0)then
            nloop=nloop+1
            mbegin(nloop)=i
c            write(*,'a7,2(i3)')'mbegin=',nloop,i
         endif
         if(sed(i).eq.1 .and.sed(i+1).eq.0 .and. nloop.gt.0)then
            mend(nloop)=i
c            write(*,'a5,2(i3)')'mend=',nloop,i
         end if
      end do

c      do i=1,nloop
c         write(*,'a6,3(x,i3)')'loops=',i,mbegin(i),mend(i)
c      enddo
ccccc calculate where each putative element of secondary structure
ccccc begins and ends
      npsec=0
c     consider third residue as beginning of the first element of
c     secondary structure      
      nprev=3
      if(mbegin(1).lt.nprev)then
         nprev=1
      endif
      do k=1,nloop
c     consider only elements bigger than minl residues
         if(mbegin(k)-nprev.gt.minl)then
            npsec=npsec+1
            msbegin(npsec)=nprev
            msend(npsec)=mbegin(k)-1
         endif
c     update beginning of next secondary structure element
         nprev=mend(k)+1 
      end do
c     consider C-terminal fragment after last loop as another putative
c     secondary element fragment. Do not take into account the last
c     three residues
      if(nres-3+1-nprev.gt.minl)then
         npsec=npsec+1
         msbegin(npsec)=nprev
         msend(npsec)=nres
      endif
      do i=1,npsec
         write(*,'i3,a10,i3,a7,i3,a5,i3')i,'SS  begin=',msbegin(i),
     &         ' end=',msend(i)
      enddo
c      write(*,*)'npsec=',npsec
ccccc assign secondary structure type for each residue
c     1:loop   2:helix   4:strand
      do i=1,nres
         nsectype(i)=getsectype(i,xb)
c         write(*,'A10,2i3')'sec-types=',i,nsectype(i)
      enddo


ccccc assign secondary structure type for each putative secondary
c     structure segment. There are several conditions to be met for
c     special cases, namely:
c     (1) the leading secondary structure assignment for the three
c     central residues is taken as the secondary structure assignment of
c     the whole block
c     (2) the length of the block should be bigger than minimal length
c     required for helix or beta
      do i=1,npsec
c     central residue, minus one
         m=(msbegin(i)+msend(i))/2-1
         write(*,'a12,4(x,i3)')'central res=',m,msbegin(i),msend(i),i
         nh=0
         nc=0
         nb=0
         do j=m,m+2
            if(nsectype(j).eq.1)then
               nc=nc+1
            elseif(nsectype(j).eq.2)then
               nh=nh+1
            else
               nb=nb+1
            endif
         enddo
         write(*,'2x,3(a3,i2,x)')'nh=',nh,'nb=',nb,'nc=',nc
c     store secondary structure type
         if(nh.gt.nc .and. nh.gt.nb)then
            if(1+msend(i)-msbegin(i).ge.minlh)then !require minimal length
               msec(i)=2
            endif
         elseif(nb.gt.nc .and. nb.gt.nh)then
            if(1+msend(i)-msbegin(i).ge.minlb)then !require minimal length
               msec(i)=4
            endif
         else
            msec(i)=1
         endif
      enddo

c     possible splitting of long alpha helixes and long strands
      npsec2=npsec
      do i=1,npsec
         do k=1,2
            minx=minS(k)
            nt=ntype(k)
c            write(*,'a5,i3')'minx=',minx
c     Do we find a long fragment ?
            if(msec(i).eq.nt.and. 1+msend(i)-msbegin(i).gt.1+2*minx)then
               nprev=msbegin(i)
               m1=msbegin(i)+minx
               m2=msend(i)-minx
c               write(*,'a8,i3,a7,i3')'msbeg=',msbegin(i),' msend=',msend(i)
c               write(*,'2(a3,i3,x),x,a6,i3')'m1=',m1,'m2=',m2,
c     &               'nprev',nprev
c     then find kinks in the fragment, ie, residues with rad(i)>0.06. kink
c     has to be at least minx residues away from the fragment termini
               do j=m1,m2
                  r=rad(j)
c                  write(*,'a2,i3,x,a5,f5.3,x,a7,f5.3')'j=',j,'radj=',
c     &                  rad(j),'radj-1=',rad(j-1)
c     signal for the beginning of the kink
                  if(r.ge.0.06 .and. rad(j-1).lt.0.06)then
                     if(nprev.eq.msbegin(i))then
                        nlast=msend(i) !remember end of primary fragment
c                        write(*,'2(a6,x,i3)')'nprev=',nprev,'nlast=',nlast
                        msend(i)=j-1 !new N-terminal fragment
                        nprev=0
                     else
                        npsec2=npsec2+1 !new interior fragment
c                        write(*,'a7,i4')'npsec2=',npsec2
                        msbegin(npsec2)=nprev
                        msend(npsec2)=j-1
                        msec(npsec2)=nt
                        nprev=0
                     endif                  
                  endif
c     signal for the end of the kink
                  if(j.gt.m1.and.r.lt.0.06.and.rad(j-1).ge.0.06)then
                     nprev=j    !mark beginning of a putative new fragment
c                     write(*,'a15,i3')'modified nprev=',nprev
                  endif
               enddo
c     trailing residues may be a new C-terminal fragment
               if(rad(m2+1).lt.0.06 .and. rad(m2).ge.0.06)then               
                  nprev=m2
c                  write(*,'a6,i3')'nprev=',nprev
               endif
               if(nprev.ne.msbegin(i))then
                  npsec2=npsec2+1 !new C-terminal fragment
                  msbegin(npsec2)=nprev
                  msend(npsec2)=nlast
                  msec(npsec2)=nt
c                  write(*,'a20')'new C-terminal fragment'
c                  write(*,'a5,i3,x,a8,i3')'npre=',nprev,'msbeg=',msbegin(i)
c                  write(*,'a6,i3,x,a6,i3')'nprev=',nprev,'nlast=',nlast
               endif
            endif
         enddo
      enddo

      npsec=npsec2 !update npsec
c      write(*,'a6,i3')'npsec=',npsec
      nhelix=0
      ncoil=0
      nbeta=0
      do i=1,npsec
c         write(*,'a5,i1')'msec=',msec(i)
         if(msec(i).eq.2)then
            nhelix=nhelix+1
            write(*,'i3,a,i3,x,a5')msbegin(i),'-',msend(i),'helix'
         elseif(msec(i).eq.4)then
            nbeta=nbeta+1
            write(*,'i3,a,i3,x,a4')msbegin(i),'-',msend(i),'beta'
         else
            ncoil=ncoil+1       !assign as coil all other possible cases
         endif
      enddo
c     number of secondary structures is alpha plus beta
      nsec=nhelix+nbeta
c      write(*,'a5,i3')'nsec=',nsec
      end
ccccc
ccccc end of function nsec(nres,xa)
ccccc =================================================



c     ==================================================
c     this function calculates secondary structure type
c     1:loop   2:helix   4:strand
c     ==================================================
      function getsectype(nresid,xb)
      parameter(nmax=2000)
      dimension xb(3,-10:nmax),r1(3),r2(3),r3(3),r4(3),r5(3)
c     obtain coordinates of several vectors
      do i=1,3
         r1(i)=xb(i,nresid)  - xb(i,nresid-1) !v_{i-1,i}=x_i - x_{i-1}
         r2(i)=xb(i,nresid+1)- xb(i,nresid)   !v_{i,i+1}=x_{i+1} - x_i
         r3(i)=xb(i,nresid+2)- xb(i,nresid+1) !v_{i+1,i+2}=x_{i+2} - x_{i+1}
         r4(i)=xb(i,nresid+2)- xb(i,nresid-1) !v_{i-1,i+2}=v_{i+2} - x_{i-1}
      enddo
c     cross-product of v_{i-1,i} x v_{i,i+1} = r1 x r2
      r5(1)=r1(2)*r2(3)-r1(3)*r2(2)
      r5(2)=r1(3)*r2(1)-r1(1)*r2(3)
      r5(3)=r1(1)*r2(2)-r1(2)*r2(1)
c     distance from i-1 to i+2
      rr=r4(1)*r4(1)+r4(2)*r4(2)+r4(3)*r4(3)
c     signed distance from i-1 to i+2 given by dot product between r5
c     and r3
      if(r5(1)*r3(1)+r5(2)*r3(2)+r5(3)*r3(3).ge.0)then
        rrs=rr
      else
         rrs=-rr
      endif
      write(*,'a7,i3,x,a3,f5.1,x,a4,f6.1')'nresid=',nresid,'rr=',
     &      rr,'rrs=',rrs
c     assign secondary structure based on rr and rrs
      if(rrs.gt.0 .and. rrs.lt.37)then
         getsectype=2 !helix
      elseif(rr.gt.74)then
         getsectype=4 !beta
      else
         getsectype=1 !coil
      endif
      end
ccccc
ccccc end of function getsectype(resid,y)
ccccc =================================================


     
c     ==================================================
c     this subroutine reads in a pdb file extracts the ca
c     coordinates					 
c     ==================================================
      subroutine import_ca_coords(name,xa,nres)
      parameter(nmax=2000)
      character*256 name
      character*4 itype
      character*1 chain
      character*3 ares,aa(20)
      character*6 atom
      dimension x(3),xa(3,nmax)
      data aa/ 'gly','ala','ser','cys','val','thr','ile','pro',
     &     'met','asp','asn','leu',
     &     'lys','glu','gln','arg',
     &     'his','phe','tyr','trp'/
     
      nres=0
      open(unit=10,file=name)

 1    continue
 100  read(10,10,err=1,end=20)atom,id,itype,
     &     ares,chain,inum,x(1),x(2),x(3)
 10   format(a6,i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3)
      if(atom.eq. 'END' )then
         go to 20
      end if
c     stop reading if 'endmdl' or 'ter' is hit
      if(atom .eq. 'ENDMDL')go to 20
      if(atom .eq. 'TER'   )go to 20	
c     read another line for the following cases
      if(atom .eq. 'HETATM')go to 100
      if(ares .eq. 'ACE'   )go to 100
      if(atom .ne. 'ATOM  ')go to 100
c     begin read in of new residue
      if(itype .eq. ' CA')then	
         nres=nres+1
         xa(1,nres)=x(1)
         xa(2,nres)=x(2)
         xa(3,nres)=x(3)
c         write(*,'i3,3f8.3')nres,xa(1,nres),xa(2,nres),xa(3,nres)
      end if
c     read another line of the file
      go to 100
c     we reached end of pdb file
 20   continue
 200  continue
      close(10)
      return
      end
