************************************************************************
*     pgf77 -g -Mbounds -Mchkfpstk -Mextend  -Wl,-static -o TMalign_lists main_fragment.f main_initial.f main_tmsc.f TMalign_lists.f
*     pgf77 -Mextend -O3 -Bstatic -o TMalign_lists  main_fragment.f main_initial.f main_tmsc.f TMalign_lists.f
*
*                          FrTM-align                               ****
*                                                                   ****
*   Program to align two protein structures, with maximizing the    ****
*   TM-score. The program has been extensively modified from the    ****
*   previous version.                                               ****
*                                                                   ****        
*   Copyright (C) Shashi Pandit and Jeffrey Skolnick                ****
*                                                                   ****
************************************************************************  

      program FrTMalign

      parameter (maxres=3000)           ! no. of residues       
      parameter (maxlen=2*maxres)       ! for alignment
      parameter (maxaln=100)            ! Maximum no. of alignments
      parameter(lmax=16501)
      character*512 pdb(2),outfile,outname,buffer,fnam
      logical       ispresent

      character*3 aanam1(-2:20),resn(maxres,0:1)
      character*1 aanam2(-2:20),seqn(maxres,0:1)
      character aseq1(maxlen),aseq2(maxlen),aseq3(maxlen)

      dimension icon(maxres,20,0:1),ncon(maxres,0:1)    ! contact information
      dimension invmap0(maxres)
      dimension xtm1(maxres),ytm1(maxres),ztm1(maxres)
      dimension xtm2(maxres),ytm2(maxres),ztm2(maxres)
      dimension m1(maxres),m2(maxres)
      logical nglob
     
      data aanam1 /'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     &              'PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN',
     &              'ARG','HIS','PHE','TYR','TRP','CYX','MSE'/

      data aanam2 /'X','G','A','S','C','V','T','I','P','M','D','N',
     &     'L','K','E','Q','R','H','F','Y','W','C','M'/

      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /pdbinfo/ ires(maxres,0:1),resn,seqn
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure
      common /n1n2/ n1(maxres),n2(maxres)
      common /dinfo/ d8
ccc   RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),r_3(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /maxres*1.0/
ccc   additional variables for handling lists of files
      character seq1(0:maxres),seq2(0:maxres)
      real xas(3,maxres,0:1,lmax) !extension to XA
      integer nseqXs(lmax,0:1) !extension to nseq1,nseq2
      character seqXs(0:maxres,lmax,0:1)!extension to seq1,seq2
      character*3 ssXs(maxres,lmax,0:1) !extension to ss1,ss2
      character*200 structDir, targDir !directories of structures and targets
      character*200 listA, listB !file lists containing basenames of structures and targets
      character*200 stru(lmax),targ(lmax) !store basenames of structures and targets
      character*200 out_file !redirect standard output to this file
      integer nstru,ntarg !number of structures and targets
      integer mmXs(maxres,lmax,0:1) !residue number
      integer readPDBfile

cccc Arguments
c      outfile='STDOUT'
      narg=iargc()
      if (narg .eq. 0) call instru
      call getarg (1,fnam)
      if (fnam.eq.'?'.or.fnam.eq.'-h'.or.fnam.eq.'--'.or.
     &     fnam.eq.'--help'.or.fnam.eq.'-help') then
         call instru
      endif

      m_out=-1                  !decided output
      m_simpl=0                 !simplified output
      mout=-1                   ! output file
      mfix=-1                   ! Fixed length TM-score
      mmat=0                    ! Transformatrix output
      no=1                      ! Number of pdbs

      i=1
      do while (i.le.narg) 
         call getarg(i,fnam)
         if (fnam.eq.'-o') then
            mout=1
            i=i+1
            call getarg(i,outname)
         elseif (fnam.eq.'-l') then
            mfix=1
            i=i+1
            call getarg(i,fnam)
            read (fnam,*)Lfix
         elseif(fnam.eq.'-outf')then
            i=i+1
            call getarg(i,out_file)
            open(unit=6,file=out_file,access='append',status='unknown')
         elseif(fnam.eq.'-sdir')then
            i=i+1
            call getarg(i,structDir)
         elseif(fnam.eq.'-tdir')then
            i=i+1
            call getarg(i,targDir)
         elseif(fnam.eq.'-ss')then
            i=i+1
            call getarg(i,listA)
            call loadFileNames(listA,stru,nstru)
         elseif(fnam.eq.'-tt')then
            i=i+1
            call getarg(i,listB)
            call loadFileNames(listB,targ,ntarg)
         elseif(fnam.eq.'-simpl')then
             m_simpl=1
         else
            read(fnam,'(a)')pdb(no)
            no=no+1
         endif         
         i=i+1
      enddo

******************************************************************
****    main program starts
******************************************************************
cc reading pdb files 

c        do i=1,2
c         ip=i-1
c         call readpdb(pdb(i),ip,nlen,aanam1,aanam2)
c         if (i.eq.1) nseq1=nlen
c         if (i.eq.2) nseq2=nlen
c        enddo

c     read data from structures

      isD=lnblnk(structDir) !number of characters
      do istru=1,nstru
         j=lnblnk(stru(istru))
         write(fnam,'a')structDir(1:isD)//'/'//stru(istru)(1:j) !absolute file name
         nseqXs(istru,0)=readPDBfile(fnam,xas,seqXs,ssXs,mmXs,0,
     $        istru)
      enddo
c     read data from targets, if we are indeed passing any target files
      isD=lnblnk(targDir) !number of characters
      do itarg=1,ntarg
         j=lnblnk(targ(itarg))
         write(fnam,'a')targDir(1:isD)//'/'//targ(itarg)(1:j) !absolute file name
         nseqXs(itarg,1)=readPDBfile(fnam,xas,seqXs,ssXs,mmXs,1,itarg)
      enddo

c     iterate over structures (and targets)
      istru=1
      if(ntarg.eq.0) then !no targets, compare among structures
         endstru=nstru-1
      else
         endstru=nstru
      endif
 555  iM=0
      nseq1=nseqXs(istru,iM)    !initialize properties of structure
      pdb(1)=stru(istru)          
      do jres=1,nseq1
         seqn(jres,0)=seqXs(jres,istru,iM)
         resn(jres,0)=ssXs(jres,istru,iM)
         ires(jres,0)=mmXs(jres,istru,iM)
         do j=1,3
            xa(j,jres,0)=xas(j,jres,iM,istru)
         enddo
      enddo      
      if(ntarg.eq.0) then !no targets, compare among structures
         iN=0 !signals structure
         itarg=istru+1
         endtarg=nstru
      else !compare structure to other targets
         iN=1 !signals target
         itarg=1
         endtarg=ntarg
      endif
 556  nseq2=nseqXs(itarg,iN)   !initialize properties of target
      if(iN.eq.0) then
         pdb(2)=stru(itarg)
      else
         pdb(2)=targ(itarg)
      endif
      do jres=1,nseq2
         seqn(jres,1)=seqXs(jres,itarg,iN)
         resn(jres,1)=ssXs(jres,itarg,iN)
         ires(jres,1)=mmXs(jres,itarg,iN)
         do j=1,3
            xa(j,jres,1)=xas(j,jres,iN,itarg)
         enddo
      enddo

      write(0,*)pdb(1)(1:12),' ',pdb(2)(1:10)
c      do jres=1,nseq1
c         write(6,'3(f8.3,x)')xa(1,jres,1),xa(2,jres,1),xa(3,jres,1)
c      enddo
c      stop

cc fixing d0 for search ONLY with d0 fixed for small protein
      d0_min=0.5
      aminlen=min(nseq1,nseq2)
      d8=1.5*aminlen**0.3+3.5      !remove pairs with dis>d8 during search 
      if(aminlen.gt.15) then
        d0=1.24*(aminlen-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
        d0=d0_min
      endif
      nseq=max(nseq1,nseq2)
       do i=1,nseq
        n1(i)=i
        n2(i)=i
       enddo
      d0_search=d0
      nglob=.false.    !assumed smaller protein is globular
        
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.3.0)d0_search=3.0

      if (nseq1.lt.nseq2) then
        ipass=0
      else
        ipass=1
      endif

      call get_rg(ipass,nglob)

      if (nglob) then
        if (d0.lt.4.5)d0_search=4.5
      endif

c      write(6,*)d0,d0_search
      call super_align (d0,d0_search,invmap0)
      write(0,*)'super align finished'
ccc resuperose to find residues d < d8
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
            m1(n_al)=i          !for recording residue order
            m2(n_al)=j
         endif
      enddo
c      write(6,*)'n_al=',n_al
      d0_input=d0
      isearch=2 
      call TMsearch (d0_input,d0_search,n_al,xtm1,ytm1,ztm1,n1,n_al,
     &     xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm,isearch) !TM-score with dis<d8 only

cccc for Output TM-score set d0 
      d0_min=0.5                !for output
      anseq=nseq1               !length for defining final TMscore
      if(mfix.eq.1)anseq=L_fix  !input length
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.3.0)d0_search=3.0
      if (nglob) then
        if (d0.lt.4.5)d0_search=4.5
      endif

cccc remove dis>d8 in normal TM-score calculation for final report
      j=0
      n_eq=0
c      write(6,*)'n_al=',n_al
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)
            m1(j)=m1(i)
            m2(j)=m2(i)
            if(seqn(m1(i),0).eq.seqn(m2(i),1) )then
               n_eq=n_eq+1
            endif
         endif
      enddo
      n8_al=j
      d0_input=d0
      isearch=3
      call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &     xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm,isearch) !normal TMscore
      rmsd8_al=Rcomm
      TMfirstseq=TM8/anseq       !TM-score after cutoff
c      write(6,*) 'TMfirstseq=',TMfirstseq,d0
c      write(6,*) 'n8_al=',n8_al


*!!!  Output TM-score is based on the second protein------------------>
      anseq=nseq2               !length for defining final TMscore
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      seq_id=float(n_eq)/(n_al+0.00000001)
      d0_input=d0
      isearch=3
      call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &     xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm,isearch) !normal TMscore
      rmsd8_al=Rcomm
      TM8=TM8/anseq       !TM-score after cutoff

********* for output summary ******************************
      if(m_simpl.eq.1)then !simplifed output
         write(6,'axaxf6.4xf6.4') pdb(1)(1:lnblnk(pdb(1))),
     $        pdb(2)(1:lnblnk(pdb(2))),TMfirstseq,TM8
         goto 8999
      endif

cccc write rotation matrix
      if (mmat.eq.1) then
       open (unit=8,file='trf.mat',status='unknown')
       L=0
        do i=1,n8_al
          k=m1(i)
          L=L+1
          r_1(1,L)=xa(1,k,0)
          r_1(2,L)=xa(2,k,0)
          r_1(3,L)=xa(3,k,0)
          r_2(1,L)=xtm1(i)
          r_2(2,L)=ytm1(i)
          r_2(3,L)=ztm1(i)
        enddo
        if(L.gt.3)then
           call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2

           write(8,*)'Rotation matrix to rotate chain 1 to chain 2'
           write(8,*)'i          t(i)         u(i,1)         u(i,2) ',
     &               '        u(i,3)'
           do i=1,3
              write(8,204)i,t(i),u(i,1),u(i,2),u(i,3)
           enddo
  204     format(I2,f18.10,f15.10,f15.10,f15.10)
        endif
      endif
      close (8)
      
cccc write coordinates

      if (mout.eq.1) then
        open (unit=9,file=outname,status='unknown')

**  rasmol   script:
         write(9,900)'load inline'
         write(9,900)'select atomno<2000'
         write(9,900)'wireframe .45'
         write(9,900)'select none'
         write(9,900)'select atomno>2000'
         write(9,900)'wireframe .20'
         write(9,900)'color white'
 900     format(A)

         do i=1,n8_al
            dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &           (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
            if(dis2.le.5)then
               write(9,901)m1(i)
               write(9,900)'color red'
               write(9,901)2000+m2(i)
               write(9,900)'color red'
            endif
         enddo
 901     format('select atomno=',I4)
         write (9,900)'select all'
         write (9,900)'exit'
         write (9,102) pdb(1),nseq1
         write (9,103) pdb(2),nseq2
 102     format ('REMARK   Protein 1:',A10,' Length =',i4)
 103     format ('REMARK   Protein 2:',A10,' Length =',i4)
         write (9,104) n8_al,rmsd8_al,TM8
 104     format ('REMARK   Structural alignment summary: ',
     &     ' Aligned length=',i4,'; RMSD =',f6.2,'; TM-score=',f7.5)      

***   chain1:
         do i=1,n8_al
            write(9,1237)m1(i),resn(m1(i),0),ires(m1(i),0),
     &           xtm1(i),ytm1(i),ztm1(i)
         enddo
         write(9,1238)          !TER
         do i=2,n8_al
            write(9,1239)m1(i-1),m1(i) !connect atoms
         enddo
***   chain2:
         do i=1,n8_al
            write(9,1237)2000+m2(i),resn(m2(i),1),ires(m2(i),1),
     $           xtm2(i),ytm2(i),ztm2(i)
         enddo
         write(9,1238)
         do i=2,n8_al
            write(9,1239)2000+m2(i-1),2000+m2(i)
         enddo
      endif
      close(9)
 1237    format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
 1238    format('TER')
 1239    format('CONECT',I5,I5)
   
cccc write the structure based sequence alignment
      write(*,*)
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)'*                               FrTM-align           ',
     &     '                    *'
      write(*,*)'* A protein structural alignment algorithm based on T',
     &     'M-score             *'
      write(*,*)'* Reference: Y. Zhang and J. Skolnick, Nucl. Acids Re',
     &     's. 2005 33, 2302-9  *'
      write(*,*)'* Comments on the program, please email to: spandit3@',
     &     'mail.gatech.edu                *'
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)
      write(*,105)pdb(1),nseq1
 105  format('Chain 1:',A10,'  Size=',I4)
      write(*,106)pdb(2),nseq2,int(anseq)
 106  format('Chain 2:',A10,'  Size=',I4,
     &     ' (TM-score is normalized by ',I4,')')
      write(*,*)
      write(*,107)n8_al,rmsd8_al,TM8,seq_id
 107  format('Aligned length=',I4,', RMSD=',f6.2,
     &     ', TM-score=',f7.5,', ID=',f5.3)
      write(*,*)

       ii=0
       i1_old=1
       i2_old=1
       do i=1,n8_al
         do j=i1_old,m1(i)-1
            ii=ii+1
            aseq1(ii)=seqn(j,0)
            aseq2(ii)='-'
            aseq3(ii)=' '
         enddo
         do j=i2_old,m2(i)-1
            ii=ii+1
            aseq1(ii)='-'
            aseq2(ii)=seqn(j,1)
            aseq3(ii)=' '
         enddo

         ii=ii+1
         aseq1(ii)=seqn(m1(i),0)
         aseq2(ii)=seqn(m2(i),1)
         dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &     (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
         if(dis2.le.5)then
           aseq3(ii)=':'
         else
           aseq3(ii)='.'
         endif
         i1_old=m1(i)+1
         i2_old=m2(i)+1
       enddo

      do i=i1_old,nseq1
         ii=ii+1
         aseq1(ii)=seqn(i,0)
         aseq2(ii)='-'
         aseq3(ii)=' '
      enddo
      do i=i2_old,nseq2
         ii=ii+1
         aseq1(ii)='-'
         aseq2(ii)=seqn(i,1)
         aseq3(ii)=' '
      enddo
      write(*,50)
 50   format('(":" denotes the residue pairs of distance < 5.0 ',
     &     'Angstrom)')
      write(*,10)(aseq1(i),i=1,ii)
      write(*,10)(aseq3(i),i=1,ii)
      write(*,10)(aseq2(i),i=1,ii)
 10   format(10000A1)
      write(*,*)

 8999 continue
      itarg=itarg+1
      if(itarg.le.endtarg) goto 556 !iterate over remaining targets
      istru=istru+1
      if(istru.le.endstru)  goto 555 !iterate over remaining structures

 9999 end











******************************************************************
****            Subroutines
******************************************************************
cccc get radius of gyration
      subroutine get_rg (np,glob)
      parameter (maxres=3000)
      double precision xc,yc,zc,disx
      logical glob
      
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2

       glob=.false.
        xc=0.0
        yc=0.0
        zc=0.0
        disx=0.0
       if (np.eq.0)nlen=nseq1
       if (np.eq.1)nlen=nseq2

       do i=1,nlen
        xc=xa(1,i,np)+xc
        yc=xa(2,i,np)+yc
        zc=xa(3,i,np)+zc
       enddo

         xc=xc/float(nlen)
         yc=yc/float(nlen)
         zc=zc/float(nlen)

       do i=1,nlen
        disx=disx+( ((xa(1,i,np)-xc)**2)
     &       +((xa(2,i,np)-yc)**2)
     &       +((xa(3,i,np)-zc)**2) )
       enddo

       rg=sqrt(disx/float(nlen))
       predRG=2.2*exp(0.38*alog(float(nlen)))

       if (rg.gt.predRG) glob=.true.    ! protein is not globular

      return
      end

cccc making initial superposition

      subroutine super_align(dx,dxs,invmap0)
      parameter (maxres=3000)
      dimension score(maxres,maxres),invmap(maxres),invmap0(maxres)
      dimension invmap_i(maxres)
      
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres) !secondary structure
      
      call assignssp            ! secondary structure assignment
      call fragscan(dx,dxs,invmap0,atm)
      
      return
      end

cccc for getting the best alignment
      subroutine getbest(aTM,invmapi,aTMmax,invmapr)
      parameter (maxres=3000)
      dimension invmapi(maxres),invmapr(maxres)
      common /length/ nseq1,nseq2

        if (aTM.gt.aTMmax) then
            aTMmax=aTM
                do j=1,nseq2
                    invmapr(j)=invmapi(j)
                enddo
        endif
      return
      end

cccc making secondary structure assignment
      subroutine assignssp
      parameter (maxres=3000)
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure

       do i=1,nseq1
          isec(i)=1
          j1=i-2
          j2=i-1
          j3=i
          j4=i+1
          j5=i+2
          if(j1.ge.1.and.j5.le.nseq1)then
             dis13=diszy(0,j1,j3)
             dis14=diszy(0,j1,j4)
             dis15=diszy(0,j1,j5)
             dis24=diszy(0,j2,j4)
             dis25=diszy(0,j2,j5)
             dis35=diszy(0,j3,j5)
             isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
          endif
       enddo
       do i=1,nseq2
          jsec(i)=1
          j1=i-2
          j2=i-1
          j3=i
          j4=i+1
          j5=i+2
          if(j1.ge.1.and.j5.le.nseq2)then
             dis13=diszy(1,j1,j3)
             dis14=diszy(1,j1,j4)
             dis15=diszy(1,j1,j5)
             dis24=diszy(1,j2,j4)
             dis25=diszy(1,j2,j5)
             dis35=diszy(1,j3,j5)
             jsec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
          endif
       enddo
       call smooth               !smooth the assignment

      return
      end

cccc make secondary str
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end

cccc smooth the secondary structure assignment
      subroutine smooth
      parameter (maxres=3000)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure

***   smooth single -------------->
***   --x-- => -----
       do i=3,nseq1
          if(isec(i).eq.2.or.isec(i).eq.4)then
             j=isec(i)
             if(isec(i-2).ne.j)then
                if(isec(i-1).ne.j)then
                   if(isec(i+1).ne.j)then
                      if(isec(i+1).ne.j)then
                         isec(i)=1
                      endif
                   endif
                endif
             endif
          endif
       enddo
       do i=3,nseq2
          if(jsec(i).eq.2.or.jsec(i).eq.4)then
             j=jsec(i)
             if(jsec(i-2).ne.j)then
                if(jsec(i-1).ne.j)then
                   if(jsec(i+1).ne.j)then
                      if(jsec(i+1).ne.j)then
                         jsec(i)=1
                      endif
                   endif
                endif
             endif
          endif
       enddo

***   smooth double -------------->
***   --xx-- => ------
       do i=1,nseq1-5
          if(isec(i).ne.2)then
            if(isec(i+1).ne.2)then
                if(isec(i+2).eq.2)then
                    if(isec(i+3).eq.2)then
                        if(isec(i+4).ne.2)then
                            if(isec(i+5).ne.2)then
                                isec(i+2)=1
                                isec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
 
          if(isec(i).ne.4)then
            if(isec(i+1).ne.4)then
                if(isec(i+2).eq.4)then
                    if(isec(i+3).eq.4)then
                        if(isec(i+4).ne.4)then
                            if(isec(i+5).ne.4)then
                                isec(i+2)=1
                                isec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
       enddo
       do i=1,nseq2-5
          if(jsec(i).ne.2)then
            if(jsec(i+1).ne.2)then
                if(jsec(i+2).eq.2)then
                    if(jsec(i+3).eq.2)then
                        if(jsec(i+4).ne.2)then
                            if(jsec(i+5).ne.2)then
                                jsec(i+2)=1
                                jsec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
 
          if(jsec(i).ne.4)then
            if(jsec(i+1).ne.4)then
                if(jsec(i+2).eq.4)then
                    if(jsec(i+3).eq.4)then
                        if(jsec(i+4).ne.4)then
                            if(jsec(i+5).ne.4)then
                                jsec(i+2)=1
                                jsec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
       enddo

***   connect -------------->
***   x-x => xxx
       do i=1,nseq1-2
          if(isec(i).eq.2)then
            if(isec(i+1).ne.2)then
                if(isec(i+2).eq.2)then
                    isec(i+1)=2
                endif
            endif
          endif
 
          if(isec(i).eq.4)then
            if(isec(i+1).ne.4)then
                if(isec(i+2).eq.4)then
                    isec(i+1)=4
                endif
            endif
          endif
       enddo
       do i=1,nseq2-2
          if(jsec(i).eq.2)then
            if(jsec(i+1).ne.2)then
                if(jsec(i+2).eq.2)then
                    jsec(i+1)=2
                endif
            endif
          endif
 
          if(jsec(i).eq.4)then
            if(jsec(i+1).ne.4)then
                if(jsec(i+2).eq.4)then
                    jsec(i+1)=4
                endif
            endif
          endif
       enddo
 
      return
      end

      
cccc distance calculation

      function diszy(i,i1,i2)
       parameter (maxres=3000)
       common /coord/ xa(3,maxres,0:1)

        diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     &     +(xa(2,i1,i)-xa(2,i2,i))**2
     &     +(xa(3,i1,i)-xa(3,i2,i))**2)
       return
      end

cccc reading pdb file
       subroutine readpdb (filen,np,nlen,aanam1,aanam2) ! filen: file name, np: 0 or 1

        parameter (maxres=3000)           ! no. of residues       

        character*100 buffer,fnam,filen
        character*3 aanam1(-2:20),resn(maxres,0:1)
        character*1 aanam2(-2:20),seqn(maxres,0:1)

        common /coord/ xa(3,maxres,0:1)
        common /pdbinfo/ ires(maxres,0:1),resn,seqn

          open (unit=10,file=filen)
          i=0
          do while(.true.) 
            read(10,('a100'),end=101) buffer
                if (i.gt.0) then
                  if(buffer(1:6).eq.'ENDMDL') then
                    write(*,*) 'Warning: PDB contains multiple models, 
     &               only first model will be used'
                  endif
                  if(buffer(1:3).eq.'TER')goto 101
                endif

                if (buffer(1:4).eq.'ATOM'.or.
     &              (buffer(1:6).eq.'HETATM'.and.buffer(17:3).eq.'MSE')) then
                 if(buffer(13:16).eq.' CA '.or. buffer(13:16).eq.'  CA'
     &              .or.buffer(13:16).eq.'CA  ') then
                   if(buffer(17:17).eq.' '.or.buffer(17:17).eq.'A') then
                     i=i+1
                     read(buffer,90)fnam,resn(i,np),fnam,ires(i,np),fnam,
     &                              xa(1,i,np),xa(2,i,np),xa(3,i,np)
                     do j=-2,20
                       if (resn(i,np).eq.aanam1(j)) seqn(i,np)=aanam2(j)
                     enddo
                   endif
                 endif
                endif

          enddo
  101  continue
   90  format(a17,a3,a2,i4,a4,3f8.3)
       close (10)
        nlen=i
       return
       end

cccc Instruction for running of the program
       subroutine instru
       
       write(*,*)
       write(*,*) 'Instructions for running TM-align program:'
       write(*,*)   
       write(*,*) 'Required arguments:'
       write(*,*)   
       write(*,*) 'Two structures. Aligns structure.pdb to target.pdb'
       write(*,*) '(By default, TM-score is normalized by the length', 
     &  'of the target.pdb)'
       write(*,*) './TMalign structure.pdb target.pdb'
       write(*,*)   
       write(*,*) 'Optional arguments:'
       write(*,*)   
       write(*,*) '-o : Output the superposed coordinates in a file'
       write(*,*) '    ./TMalign structure.pdb target.pdb -o TM.sup'
       write(*,*) '    (To view the superimposed structures by rasmol:',
     &  './rasmol -script TM.sup)'
       write(*,*)   
       write(*,*) '-l : TM-score to be normalized by the given length'
       write(*,*) '    ./TMalign structure.pdb target.pdb -L 100'
       write(*,*)   
       write(*,*) '-m : <0 or 1 > If -m = 1, stores the rotation'
       write(*,*) '     matrix in trf.mat file. Default value is 0.'
       write(*,*)   
       write(*,*)'4. If you want TM-align of many structures '//
     $      'among themselves, pass a list and\n    directory'//
     $      ' of structures:'
       write(*,*)'  >TMalign -ss strList -sdir structDir'
       write(*,*)'   If you want TM-align of many structures',
     $      'to many targets:'
       write(*,*)'  >TMalign -ss strList -tt targList -sdir ',
     $      'structDir -tdir targDir'
       write(*,*)'   (structDir and targDir can be the same)'
       write(*,*)
       write(*,*)'5. Redirect output to file:'
       write(*,*)'  >TMalign ... -outf outputFile'
       write(*,*)
       write(*,*)'6. Simplifed output, "structure target TM-score":'
       write(*,*)'  >TMalign ... -simpl'
       write(*,*)
       write(*,*) '    ./TMalign structure.pdb target.pdb -m 1'
       write(*,*) 
       
       stop
       end

cccc 

      subroutine loadFileNames(list,files,nfiles)
c     load contents of file "list" into memory (files)
      PARAMETER(lmax=16501) !maximum number of model or native structures
      character*200 list,files(lmax)
      integer nfiles
      open(unit=10,file=list,status='old',err=666)
      nfiles=1
 10   read(10,end=101,err=101,*) files(nfiles)
      nfiles=nfiles+1
      goto 10
 101  close(10)
      nfiles=nfiles-1
      return
 666  call exit(1)
      end



      function readPDBfile(fnam,xas,seqXs,ssXs,mmXs,iM,iconf)
      PARAMETER(maxres=3000) !maximum number of residues
      PARAMETER(lmax=16501) !maximum number of model or native structures
      character*512 fnam          !pdb file name
      real xas(3,maxres,0:1,lmax)       !array to store coordinates
      character seqXs(0:maxres,lmax,0:1) !array to store sequence in numeric code
      character*3 ssXs(maxres,lmax,0:1)  !array to store sequence in three-letter code
      integer mmXs(maxres,lmax,0:1)      !array to store residue numbers
      integer iconf
      
      character*100 s,du
      character*3 anam
      integer readPDBfile,ires,jres
      dimension mm1(maxres)
      character*3 aa(-1:20)
      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     $     'PRO','MET','ASP','ASN','LEU',
     $     'LYS','GLU','GLN','ARG',
     $     'HIS','PHE','TYR','TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C','V','T','I',
     $     'P','M','D','N','L','K','E','Q','R',
     $     'H','F','Y','W','C'/      

 104  format(A100)
 9000 format(A17,A3,A2,i4,A4,3F8.3)
      ires=0
      open(unit=10,file=fnam,status='old',err=666)
 101  read(10,104,end=102) s
         if(s(1:3).eq.'TER' .or. s(1:3).eq.'END') goto 102
         if(s(1:4).eq.'ATOM')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &           eq.'  CA')then
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  ires=ires+1
                  if(ires.ge.maxres)goto 102 !too long sequence. Consider only first maxres res
                  read(s,9000)du,anam,du,mmXs(ires,iconf,iM),du,
     $                 xas(1,ires,iM,iconf),xas(2,ires,iM,iconf),
     $                 xas(3,ires,iM,iconf) 
                  do j=-1,20    !residue type in numeric code
                     if(anam.eq.aa(j))seqXs(ires,iconf,iM)=slc(j)
                  enddo
                  ssXs(ires,iconf,iM)=anam !residue type in three-letter code
               endif
            endif
         endif
      goto 101
 102  close(10)
      readPDBfile=ires
      return
 666  write(0,'a')'ERROR: reading',fnam
      call exit(1)
      end

