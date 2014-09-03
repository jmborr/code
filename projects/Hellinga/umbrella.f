ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set of routines to implement umbrella sampling in native rmsd
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     initialize arrays
      subroutine umbrella_init_params(n_repf,natf)
      include 'umbrella.inc.f'
      integer _ER_,i,report_exception,rmsdinf
      character*256 line,natf
      real x,y,z,rmsdmax,rmsdmin
      rmsdmin=2.0
      rmsdmax=8.0
      rmsdinf=1000.0
      !init number of replicas that we store in zrep, "local variable" to the umbrella "namespace"
      zrep=n_repf
      !write(6,'a,i2')'zrep=',zrep
      !init rmsd targets for each replica
      do i=1,numrepmax
c         rmsdtarget(i)=(rmsdmax/rmsdmin)**( (1.0*(i-1))/(zrep-1) ) *
c     &        rmsdmin
         rmsdtarget(i)=rmsdmin+(rmsdmax-rmsdmin)*((1.0*(i-1))/(zrep-1))
         rmsdprev(i)=rmsdinf
         !write(6,'a,i3,a,f9.2')'rmsdtarget(',i,')=',rmsdtarget(i)
      enddo
      !init native coordinates
      nres=0
      _ER_=_ER_NAT_ !set global variable error code
      _fu_nat_=1
      open(unit=_fu_nat_,file=natf,status='old',err=666)
      do while(1)
         read(_fu_nat_,end=10,'a') line
         !write(6,'a,a')'line=',line
         if(line(1:5).eq.'ATOM'.and.line(13:17).eq.' CA ')then
            nres=nres+1
            !write(6,'a,i3')'nres=',nres
            if(nres.gt.numresmax)then !sequence length exceeds allocated length
               _ER_=_ER_TOOBIG_
               goto 666
            endif
            !write(6,'a1,a,a1')'|',line(31:54),'|'
            !call exit(1)
            read(line(31:39),'f8.3') x
            read(line(39:47),'f8.3') y
            read(line(47:54),'f8.3') z
            natR(1,nres)=x
            natR(2,nres)=y
            natR(3,nres)=z
            natxyz(1,nres)=x
            natxyz(2,nres)=y
            natxyz(3,nres)=z
         endif
      enddo
 10   close(_fu_nat_)
      !do i=1,nres
      !   write(6,*) natxyz(1,i),natxyz(2,i),natxyz(3,i)
      !enddo
      !call exit(1)
      return
 666  i=report_exception(_ER_)
      !call exit(i)
      end

c     pass replica number to variable irep local to umbrella "namespace"
      subroutine umbrella_init_rep(itemp)
      include 'umbrella.inc.f'
      irep=itemp
      end

c     accept or reject a new conformation based on new rmsd value
      subroutine umbrella_reject_accept(id,x,y,z,ex,ey,ez,mv)
      include 'umbrella.inc.f'
      parameter(ndim=1000)
      integer ires,id,ier,x(ndim),y(ndim),z(ndim),mv(ndim)
      real ex(ndim),ey(ndim),ez(ndim)
      real rmsd,rmsdtarget_irep,rmsdtarget_irep2,prob,latticespacing
      double precision r(3,numresmax),u(3,3),t(3),rms,w(numresmax) !armsd is real
      data w /numresmax*1.0/

c      do i=1,nres
c         if(mv(i).lt.0)then
c            write(6,'i3,x,3(f8.3)'),i,ex(i)*0.87,ey(i)*0.87,ez(i)*0.87
c         else
c            write(6,'i3,x,3(f8.3)')i,x(i)*0.87,y(i)*0.87,z(i)*0.87
c         endif
c      enddo
c      write(6,*)'reject_accept'
c      call exit(1)

      latticespacing=0.87
      !write(6,'a,i2') 'irep=',irep
      rmsdtarget_irep=rmsdtarget(irep)
      rmsdtarget_irep2=rmsdtarget(irep+1)
      !write(6,*)rmsdtarget_irep,rmsdtarget_irep2
      id=3 !init as reject
c     calculate rmsd to native
c      write(6,*) 'latticespacing=',latticespacing
      do ires=1,nres
         if(mv(ires).lt.0)then
            !write(6,*)'juju',ires,ex(ires),ey(ires),ez(ires)
            r(1,ires)=ex(ires)*latticespacing
            r(2,ires)=ey(ires)*latticespacing
            r(3,ires)=ez(ires)*latticespacing
         else
            !write(6,*)ires,x(ires),y(ires),z(ires)
            r(1,ires)=x(ires)*latticespacing
            r(2,ires)=y(ires)*latticespacing
            r(3,ires)=z(ires)*latticespacing
         endif
      enddo
c      do ires=1,nres
c         write(6,'i3,x,3(f8.3)')ires,r(1,ires),r(2,ires),r(3,ires)
cc         write(6,*)ires,r(1,ires),r(2,ires),r(3,ires)
c      enddo
      call u3b(w,r,natxyz,nres,0,rms,u,t,ier) !u rotate r_1 to r_2
      rmsd=dsqrt(rms/nres)
c      write(6,*) 'rmsd=',rmsd
      rmsdnew=rmsd              !store for putative update of rmsdprev array
      if(irep.gt.zrep)then !no umbrella sampling above recorded conformations
         id=1
      else if(rmsd.lt.rmsdtarget_irep)then !smaller rmsd than target rmsd
         id=1 !accept
c     rmsd bigger than rmsdtarget but decreasing. This may happen if
c     prevrmsd_irep>rmsdtarget_irep because of a replica swap
      else if(rmsdtarget_irep.lt.rmsd.and.rmsd.lt.rmsdprev(irep))then
         id=1
      else if(rmsd.lt.rmsdtarget_irep2)then
         prob=aranzy() !uniform probability distribution generator in (0,1)
         if(prob.lt.(rmsd-rmsdtarget_irep)/
     &        (rmsdtarget_irep2-rmsdtarget_irep))then
            id=1
         endif
      else
         rmsdnew=-1 !don't accept the move
      endif
c      write(6,'a,i3,5(a,f5.1),a,i1')' irep=',irep,' rmsdtarget_irep=',
c     &     rmsdtarget_irep,' rmsd=',rmsd,' rmsdtarget_irep2=',
c     &     rmsdtarget_irep2,' rmsdprev=',rmsdprev(irep),' rmsdnew=',
c     &     rmsdnew,' id=',id
c      call exit(1)
      end


c     update rmsdprev array only if move is accepted
      subroutine update_rmsdprev(id)
      include 'umbrella.inc.f'
      if(id.eq.1) rmsdprev(irep)=rmsdnew
      end


c     retirm rmsd or current replica
      function armsdprev(k)
      include 'umbrella.inc.f'
      integer k
      armsdprev=rmsdprev(k)
      end


c     return lower rmsd target of the replica
      function armsdtarget(k)
      include 'umbrella.inc.f'
      integer k
      armsdtarget=rmsdtarget(k)
      end


c     exception handling function
      function report_exception(_ER_)
      include 'umbrella.inc.f'
      integer _ER_,report_exception,error_code
      character*512 mesg,mesg0
      mesg0='ERROR from umbrella::report_exception: '
      if(_ER_.eq._ER_NAT_)then !global variable _ER_ contains error code
         mesg=mesg0//'reading native file'
      else if(_ER_.eq._ER_TOOBIG_)then
         mesg=mesg0//'protein size exceeds parameter numresmax'
      else
         mesg=mesg0//'undefined error'
      endif
      write(6,*) mesg
      report_exception=_ER_
      end


c     ###########################################################
c     swap rmsd values for replicas
      subroutine umbrella_swap(irep1,irep2)
      include 'umbrella.inc.f'
      integer irep1,irep2,irep3
      irep3=rmsdprev(irep1)
      rmsdprev(irep1)=rmsdprev(irep2)
      rmsdprev(irep2)=irep3
      end


c     ###########################################################
c     return one particular coordinate (x, y, or z) at residue ii
      function umbrella_getNat(coord_number,ii)
      include 'umbrella.inc.f'
      real umbrella_getNat
      integer coord_number,ii
      umbrella_getNat=natR(coord_number,ii)
      end
