c #########################################################
c     this program reads in a pdb file extracts the coordinates
c     units are in lattice units(1lu=1.22 Angstroms=
      subroutine getCAcoords(pdbf,xyz,nres,aaa,chainID)
      parameter(nmax=2000)
      character*255 pdbf
      character *4 itype
      character*1 chainID,chain
      character*3 ares
      character aaa(3,nmax)
      character*6 atom
      integer nres
      DIMENSION X(3),xyz(3,nmax)
      do j=1,3      
         do i=1,nmax
            xyz(j,i)=0
         end do
      end do
c      write(6,*)pdbf//'.pdb is being analyzed'
      icount=0
      open(unit=10,file=pdbf,err=666)
      ntest=0
      nres=0
      inumold=-200
      ntest2=0
      icount=0	
 1    continue
      ntest=ntest+1
c     if(ntest .gt.nmax) go to 20
 100  read(10,10,err=1,end=20)atom,id,itype,
     &     ares,chain,inum,x(1),x(2),x(3)
 10   format(a6,i5,1x,a4,1x,a3,1x,a1,i4,1x,3x,3f8.3)
c     if(atom.eq. 'end' .or. atom .eq. 'ter')then
      if(atom.eq. 'END' )then
c         write(6,*)'end of file'
         go to 20
      end if
      if(chainID .ne.'_')then
         if(chain .ne. chainID)go to 1
         if(chain .eq.chainID)then
            if(atom .eq. 'ENDMDL')go to 20
            if(atom .eq.'TER')go to 20	
         end if
      else
         if(atom .eq.'TER')go to 20
         if(atom .eq. 'ENDMDL')go to 20
      end if
      if(atom .eq. 'HETATM')go to 1
      if(ares .eq. 'ACE')go to 1
      if(atom .ne. 'ATOM  ')go to 1
c     begin read in of new residue
      if(itype .eq. ' CA')then	
         icount=icount+1
         aaa(1:3,icount)=ares
         xyz(1,icount)=x(1)
         xyz(2,icount)=x(2)
         xyz(3,icount)=x(3)
      end if
      go to 100
 20   continue
 200  continue
      nres=icount
c      write(6,*)pdbf,' which has',nres,' residues is finished'
      close(10)
      return
 666  write(6,*)'ERROR from getCAcoords: cannot read '//pdbf
      END
