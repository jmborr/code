c     pgf77 -Mextend -O -s -fast -Wl,-static -o pdb2ss_tasserdef.x pdb2ss_tasserdef.f getCAcoords.f inputArgs.f f2kcli.f
      program pdb2ss_tasserdef
      parameter(nmax=1000) !maximum number of residues
      integer out_p,nres
      character pdbf*255,outf*255,aaa(3,nmax),a1
      real xyz(3,nmax)
      call read_command_line(pdbf,outf,out_p) !read input arguments
      call getCAcoords(pdbf,xyz,nres,aaa,a1)   !input CA coordinates from structure
      
      open(unit=out_p,file=outf) !write output
      
      close(out_p)
      end !end of program pdb2ss



c     ###########################################
      subroutine read_command_line(pdbf,outf,out_p)
c     these are the variables to be initialized from command line
      character*256 pdbf,outf
      integer out_p
c     variables intrinsic to the subroutine
      character inputs*255(26)
      character wellcome*255(0:26)
      integer c2i
c     the wellcome message
      wellcome(0)='./pdb2ss_tasserdef.x [options]'
      wellcome(1)='-a _RA_pdbf pdb file'
      wellcome(2)='-b __outf output file (def:stdout)'
c     pass command values from the command line
      call proccess_command_line(wellcome,inputs)
      read(inputs(c2i('a')), *) pdbf
      out_p=6 !defaults to standard output
      if(inputs(c2i('b'))(1:1).ne.'')then
         out_p=1
         read(inputs(c2i('b')), *) outf
      endif
      end
      

c     ##########################################
      subroutine output_ss(xyz,nres)
      paramter(nmax=1000)
      real xyz(3,nmax),ax(3,5),ar15,dot13,dot24,dot14
      integer nres,i,j,k,r14
      character ss,sss*nmax
c     artificially extend last five residues
      do i=1,3
         do j=1,5
            xyz(i,nres+j)=xyz(i,nres+j-5)
         enddo
      enddo
c     cycle through sequence
      do i=1,nres
         ss='C'   !init to coil
         do k=1,3 !init ax
            do j=1,5
               ax(k,j)=ayz(k,i+j-1)
            enddo
         enddo
         ar15=0
         do k=1,3
            ar15=ar15+(ax(k,1)-ax(k,5))**2
         enddo
         if(ar15.lt.57)then     !7.53A: helix
c ?????? WHO IS IBIN AND R14 ???????
            if(ibin(r14).gt.4.and.ibin(r14).lt.8)then !right-hand,3A->8A
               dot13=0.0
               do k=1,3
                  dot13=dot13+(ax(k,4)-ax(k,3))*(ax(k,2)-ax(k,1))
               enddo
               if(dot13.lt.0)then
                  dot24=0.0
                  do k=1,3
                     dot24=dot24+(ax(k,5)-ax(k,4))*(ax(k,3)-ax(k,2))
                  enddo
                  if(dot24.lt.0)then
                     dot14=0.0
                     do k=1,3
                        dot14=dot14+(ax(k,5)-ax(k,4))*(ax(k,2)-ax(k,1))
                     enddo
                     if(dot14.gt.0)then
                        ss='H'
                     endif
                  endif
               endif
            endif
         elseif(ar15.gt.121)then !11A, beta
            if(sec(i+1).ne.2.AND.sec(i+2).ne.2.AND.sec(i+3).ne.2)then
               if(mv(i+1).gt.0)then
                  bx2=hbx(ica(i),ica(i+1))
                  by2=hby(ica(i),ica(i+1))
                  bz2=hbz(ica(i),ica(i+1))
               else
                  bx2=ebx(i+1)  !Hb
                  by2=eby(i+1)
                  bz2=ebz(i+1)
               endif
               if(mv(i+3).gt.0)then
                  bx4=hbx(ica(i+2),ica(i+3))
                  by4=hby(ica(i+2),ica(i+3))
                  bz4=hbz(ica(i+2),ica(i+3))
               else
                  bx4=ebx(i+3)
                  by4=eby(i+3)
                  bz4=ebz(i+3)
               endif
               dot24=bx2*bx4+by2*by4+bz2*bz4 !cos(H1.H3)
               if(dot24.gt.0.71)then !angle<45
                  if(mv(i+2).gt.0)then
                     bx3=hbx(ica(i+1),ica(i+2))
                     by3=hby(ica(i+1),ica(i+2))
                     bz3=hbz(ica(i+1),ica(i+2))
                  else
                     bx3=ebx(i+2)
                     by3=eby(i+2)
                     bz3=ebz(i+2)
                  endif
                  dot23=bx2*bx3+by2*by3+bz2*bz3 !cos(H1.H2)
                  if(dot23.lt.-0.71)then !angle >135
                     ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
c                     ESHORT5c=ESHORT5c-2-ff
                     ESHORT5c=ESHORT5c-2-ff*2
                  endif
               endif
            endif
         endif
      enddo
      end
