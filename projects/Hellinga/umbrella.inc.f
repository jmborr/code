c     set of parameters and common blocks for umbrella sampling in native rmsd
      parameter(numrepmax=100) !maximum number of replicas
      parameter(numresmax=300) !maximum number of residues
      common/umbrella/rmsdtarget,rmsdprev,rmsdnew,natxyz,nres,irep,zrep,
     &     natR
      integer nres,irep,zrep
      real rmsdtarget(numrepmax),rmsdprev(numrepmax),rmsdnew
      real natR(3,numresmax)
      double precision natxyz(3,numresmax)
c     file units
c     _fu_nat_=1
      common/umbrella_fu/_fu_nat_
      integer _fu_nat_

c     error codes
c     _ERR_NAT_=1 ,_ER_TOOBIG_=2
      common/umbrella_err/_ERR_NAT_,_ER_TOOBIG_
      integer _ERR_NAT_,_ER_TOOBIG_
