********************************************************************************
****               Release 1.0
****
**** IMPORTANT NOTICE:   
**** This routine is part of development. Hence, there are many redundant parts, 
**** and also parts, which are not finally used. This will be cleaned
**** in next version of the program. Due to this issue, progam virtual
**** memory usuage is higher. However, the usuage of RESIDENT memory is < 100MB (typically)
**** The JOBS can be submitted with memory requirement of 100mb
****
**** Subroutine for scanning fragment based for structural alignment
**** Split each sequence into set length. Now, using the fragments of
**** second sequence scan the first sequence fragments. Join the
**** fragments with DP (using GL/SSP matrices) and make the initial alignment.
**** This will be extened using normal iteration to obtain best TM-score
********************************************************************************

      subroutine fragscan (dx,dxs,mapr,TMmax)
      parameter (maxres=3000)
      parameter (maxfr=1000)    ! maximum no. of fragments
      parameter (nfraln=6000)   ! maximum no. of frag. combination used for
                                ! dpsc*, ssp
      parameter (maln=200)      ! maximum of frag. combination allowed 
                                ! as set by maxgl
      parameter (mfrali=100)   ! max. TOTAL no. of alignments of FRAGMENTS

      dimension invmap(maxres),mapA(maxres),mapr(maxres)
      dimension score(maxres,maxres),scorei(maxres,maxres),gap(10)

      dimension ista(maxfr,2),iend(maxfr,2) ! start and end residues information
      dimension iali(maln,2),iiali(maln,2),tmal(maln),inum(maln,2)
      dimension iinum(maln,2)
      
      dimension dpsc(nfraln,nfraln),dpsccp(nfraln,nfraln),ssp(nfraln,nfraln)
      dimension smtx(nfraln,nfraln)
      dimension ifrali(maxfr),jali(mfrali,maxfr)
     
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)

      imin=min(nseq1,nseq2)
      imax=max(nseq1,nseq2)
      fr=float(imin)/float(imax)

      call assp (asspt1,asspt2,nlen)

      if (nlen.lt.6) then
        call getlen (fr,nlen,asspt1,asspt2) ! get length of the fragme
      endif

      nshift=4
      if (fr.gt.0.70)nshift=3
      nfr1=nseq1/nlen
      nfr2=nseq2/nlen
      
      do kk=1,2
        if (kk.eq.1) ntot=nfr1
        if (kk.eq.2) ntot=nfr2
          do j=1,ntot
            ista(j,kk)=(nlen)*(j-1)+1
            iend(j,kk)=ista(j,kk)+(nlen)-1
          enddo
      enddo
      nfrmax=max(nfr1,nfr2)

      call initialize (dpsc,dpsccp,ssp,smtx,nfrmax) !useful when doing TMalign for many proteins

      maxgl=60
      ntotal=0
      icou=0

      nfr1_mod=0
      nfr2_mod=0
      asec_max=0    ! maximum sec. str. matching content in the fragment.

      do i=1,nfr2
       nfr2_mod=i

        do j=1,nfr1

            do 111 kno=1,nshift ! two shifts
              call fillinvmap(invmap)

              if (kno.eq.1) then 
                istart1=ista(j,1)
                nfr1_mod=(j-1)*2+1
                nfr2_mod=(i-1)*2+1
                istart2=ista(i,2)
              elseif (kno.eq.2) then
                istart1=ista(j,1)+(nlen/2)
                istart2=ista(i,2)
                if ((istart1+nlen).gt.nseq1) then 
                  goto 111
                endif
                nfr1_mod=(j-1)*2+2
                nfr2_mod=(i-1)*2+1
              elseif(kno.eq.3) then
                istart1=ista(j,1)
                istart2=ista(i,2)+(nlen/2)
                if ((istart2+nlen).gt.nseq2) then
                    goto 11
                endif
                nfr1_mod=(j-1)*2+1
                nfr2_mod=(i-1)*2+2
              elseif (kno.eq.4) then
                istart1=ista(j,1)+(nlen/2)
                istart2=ista(i,2)+(nlen/2)

                if ((istart1+nlen).gt.nseq1) goto 111
                if ((istart2+nlen).gt.nseq2) goto 11

                nfr1_mod=(j-1)*2+2
                nfr2_mod=(i-1)*2+2
              endif

                itmpst=istart1
                iss=0
                do kk=istart2,(istart2+nlen-1)
                    invmap(kk)=itmpst
                        if (isec(itmpst).eq.jsec(kk))iss=iss+1
                    itmpst=itmpst+1
                enddo
                ntotal=ntotal+1
               
                call get_GL (dxs,dx,invmap,GL)

                dpsc(nfr1_mod,nfr2_mod)=GL
                asec=(float(iss)/float(nlen))  ! sec. str.
                ssp(nfr1_mod,nfr2_mod)=asec    ! sec. str.
                if (asec.gt.asec_max) asec_max=asec    ! sec. str.

                 if (ntotal.lt.(maxgl+1)) then
                    tmal(ntotal)=GL
                    iali(ntotal,1)=nfr1_mod ! fragment no. for protein 1
                    iali(ntotal,2)=nfr2_mod ! fragment no. for protein 2
                    inum(ntotal,1)=istart
                    inum(ntotal,2)=istart2
                        if (ntotal.eq.maxgl) then
                           call sortn (maxgl,tmal,iali,inum,iiali,iinum)
                        endif
                 else
                        do ii=maxgl,1,-1
                            if(GL.gt.tmal(ii)) then
                                do jj=1,ii-1
                                    iiali(jj,1)=iiali(jj+1,1)
                                    iiali(jj,2)=iiali(jj+1,2)
                                    iinum(jj,1)=iinum(jj+1,1)
                                    iinum(jj,2)=iinum(jj+1,2)
                                    tmal(jj)=tmal(jj+1)
                                enddo
                                iiali(ii,1)=nfr1_mod
                                iiali(ii,2)=nfr2_mod
                                iinum(ii,1)=istart    ! residue no. for protein 1
                                iinum(ii,2)=istart2   ! residue no. for protein 2
                                tmal(ii)=GL
                              goto 12
                            endif
                        enddo
                 endif
   12                          continue                 
                 
 111        continue
   11           continue
        enddo
        if (icou.eq.0) then
            nfr1_tot=nfr1_mod
            icou=1
        endif
      enddo
      nfr2_tot=nfr2_mod

       if (ntotal.lt.maxgl) then
           maxgl=ntotal
           call sortn (maxgl,tmal,iali,inum,iiali,iinum)
       endif

       do ii=1,nfr1_tot
           do jj=1,nfr2_tot
               dpsc(ii,jj)=(dpsc(ii,jj)/tmal(maxgl))
           enddo
       enddo
        
cccc  Scanning the paths generated by DPSC or SSP scores using DP with
cccc  three gap opening penalty. and each alignment is passed for
cccc  extension with iteration.

       TMmax=0.0

       do idp=1,2

        if (idp.eq.1) then
            call cpdpsc (nfr1_tot,nfr2_tot,dpsc,dpsccp)
            niter=6      !change from 5 to 2
        elseif (idp.eq.2) then
            call cpdpsc (nfr1_tot,nfr2_tot,ssp,dpsccp)
            niter=6
            if (fr.lt.0.26.and.(imin.ge.45.and.imin.lt.80)) niter=8
        endif

         nali=0
         igp_st=0
         if (idp.eq.2) igp_st=-1

          do  igp=-1,1
            call cpdpsc (nfr1_tot,nfr2_tot,dpsccp,smtx)
            gap_open=igp
            
            do iter=1,niter
             call dpnew (nfr1_tot,nfr2_tot,smtx,gap_open,ifrali,dpout)

              if (nali.eq.0) then
                  nali=nali+1
                  do kk=1,nfr2_tot
                    jali(nali,kk)=ifrali(kk)
                  enddo
              else
                do ialn=1,nali
                    inos=0
                    do kk=1,nfr2_tot
                      if (ifrali(kk).eq.jali(ialn,kk)) inos=inos+1
                    enddo
                    if (inos.eq.nfr2_tot) then
                        goto 18
                    endif
                enddo

                   nali=nali+1
                   do kk=1,nfr2_tot
                      jali(nali,kk)=ifrali(kk)
                   enddo
              endif

              do kk=1,nfr2_tot
                ipp=ifrali(kk)
                if (ipp.gt.0) smtx(ipp,kk)=0.0
              enddo
            enddo
   18     continue
          enddo
         call calbesttm (dx,dxs,nlen,nali,jali,nfr1_tot,nfr2_tot,dpsccp,fr,mapA,rTM)

         if (rTM.gt.TMmax) then
            TMmax=rTM
            do jj=1,nseq2
                mapr(jj)=mapA(jj)
            enddo
         endif

       enddo

      return
      end

***** routine for sorting the alignment with GLscore

      subroutine sortn (naln,gscore,ialix,inumx,iialix,iinumx)
      parameter (maxaln=200)    ! maximum no. of alignments
      dimension ialix(maxaln,2),gscore(maxaln) ! store alignments and scores
      dimension iswap(maxaln)   ! for swapping alignments
      dimension iialix(maxaln,2)   ! output alignment
      dimension iinumx(maxaln,2),inumx(maxaln,2)  ! output alignment

        do i=1,naln
            iswap(i)=i
        enddo

        do i=1,naln
           do j=1,naln-i
            if (gscore(j+1).lt.gscore(j)) then
              tmp=gscore(j)
              gscore(j)=gscore(j+1)
              gscore(j+1)=tmp
              ii=iswap(j)
              iswap(j)=iswap(j+1)
              iswap(j+1)=ii
            endif
           enddo
        enddo

cc swapping the alignment
        do i=1,naln
          iialix(i,1)=ialix(iswap(i),1)
          iialix(i,2)=ialix(iswap(i),2)
          iinumx(i,1)=inumx(iswap(i),1)
          iinumx(i,2)=inumx(iswap(i),2)
        enddo

      return
      end

***** filling the invmap
       subroutine fillinvmap (invmapr)
       parameter (maxres=3000)
       dimension invmapr(maxres)
       common /length/ nseq1,nseq2

       do jj=1,nseq2
        invmapr(jj)=-1
       enddo

       return
       end

**** rotuine for doing DP on the fragments
      subroutine dpnew (nfr1,nfr2,score,gap_open,ifralir,dpout)
      parameter (maxfr=1000)    ! maximum no. of fragments
      parameter (maxa=6000)             ! no. of fragments
      dimension score(maxa,maxa),ifralir(maxfr),ifralio(maxfr),scorea(maxa,maxa)
      dimension scorei(maxa,maxa)
      logical dir
      dimension dir(0:maxa,0:maxa),val(0:maxa,0:maxa)
      real h,v

cccc initialize matrix

      istate=0

      if (nfr2.gt.nfr1) then    ! forces the alignment to be done with nfr2 being the smallest 
                                ! sequence to avoid many altenative paths when doing DP other way around

        istate=1   ! swap the alignments in the end

        call cpdpsc(nfr1,nfr2,score,scorei)

        do ii=1,nfr1
          do jj=1,nfr2
            scorea(jj,ii)=score(ii,jj)
          enddo
        enddo

        ntemp=nfr2
        nfr2=nfr1
        nfr1=ntemp
        do kk=1,nfr1
            ifralio(kk)=0
        enddo

        do ii=1,nfr1
          do jj=1,nfr2
            score(ii,jj)=scorea(ii,jj)
          enddo
        enddo
      endif
      
       val(0,0)=0.0
       do i=1,nfr1
         dir(i,0)=.false.
         val(i,0)=0.0
       enddo
       do j=1,nfr2
         dir(0,j)=.false.
         val(0,j)=0.0
         ifralir(j)=-1
       enddo

cc fill the matrix and path

       do j=1,nfr2
        do i=1,nfr1
           d=val(i-1,j-1)+score(i,j)
           h=val(i-1,j)
           v=val(i,j-1)

            if(dir(i-1,j))h=h+gap_open
            if(dir(i,j-1))v=v+gap_open

                if((d.ge.h).and.(d.ge.v)) then
                    dir(i,j)=.true.
                    val(i,j)=d
                else
                    dir(i,j)=.false.
                        if(v.ge.h)then
                            val(i,j)=v
                        else
                            val(i,j)=h
                        end if
                endif
        enddo
       enddo

cccc   extract the alignment:
       i=nfr1
       j=nfr2

       do while((i.gt.0).and.(j.gt.0))

         if ( dir(i,j) ) then
           ifralir(j)=i
           ifralio(i)=j
           i=i-1
           j=j-1
         else
           h=val(i-1,j)
           v=val(i,j-1)

           if(dir(i-1,j))h=h+gap_open
           if(dir(i,j-1))v=v+gap_open

           if(v.ge.h) then
             j=j-1
           else
             i=i-1
           endif
         endif
       enddo

       j=nfr2
       do while(ifralir(j).lt.1)
        j=j-1
       enddo
       dpout=val(ifralir(j),j)
       

       if (istate.eq.1) then
            do jj=1,nfr1
                if (ifralio(jj).eq.0) then
                    ifralir(jj)=-1
                else
                    ifralir(jj)=ifralio(jj)
                endif
            enddo
        ntemp=nfr2
        nfr2=nfr1
        nfr1=ntemp
        call cpdpsc(nfr1,nfr2,scorei,score) 
       endif
       
       return
       end
cccc

***** Note: This will be obsolete in the newer version of the program
cccc Subroutine to get the length of the fragment 
cccc such that fragments covers most of the length of the protein with
cccc lots of optimization to acheive minimum difference between old and new
cccc TMalign programs

      subroutine getlen (fri,nlen,asspt1i,asspt2i)
      common /length/ nseq1,nseq2

      minlen=min(nseq1,nseq2)
      maxlen=max(nseq1,nseq2)
      isum=nseq1+nseq2
      imin0=100
      nlen=0
      imax=0
      imax1=0

       if (maxlen.le.100) then

        do ii=8,14,2
            ipp=isum-((isum/ii)*ii)
            if ((ipp).gt.(ii/2)) ipp=ipp-(ii/2)
            if (ipp.le.imin0) then
                imin0=ipp
                nlen=ii
            endif
        enddo

       else
        istep=2

        if (minlen.le.50) then
          if (fri.gt.0.10) then
            ist=8   
            jj=14
          else
            ist=8  
            jj=12
          endif
        elseif (minlen.gt.50.and.minlen.le.90) then
          if (fri.gt.0.30) then
            ist=8
            jj=14
          else
            ist=8
            jj=18
          endif
        elseif (minlen.gt.90.and.minlen.le.200) then
            ist=10
            jj=18
        else
          if (fr.lt.0.20) then
            jj=10
            ist=minlen/3
            if (ist.gt.60)ist=60
            if (mod(ist,2).ne.0) ist=ist+1
            istep=-2
          else
            ist=8
            jj=18
          endif
        endif
        
        do ii=ist,jj,2
            iret1=(nseq1/ii)*ii
            iret2=(nseq2/ii)*ii
            idiff=(nseq1+nseq2)-(iret1+iret2)
            if ((nseq1-iret1).ge.(ii/2))then
                idiff=idiff-(ii/2)
                iret1=iret1+(ii/2)
            endif
            if ((nseq2-iret2).ge.(ii/2)) then
                idiff=idiff-(ii/2)
                iret2=iret2+(ii/2)
            endif
                iret1=nseq1-iret1
                iret2=nseq2-iret2

            if (minlen.le.50)then
              if (fri.lt.0.10.and.asspt1i.lt.0.10) then
               if (idiff.le.imin0) then
                  imin0=idiff
                  nlen=ii
               endif
              else
              if (idiff.le.imin0.and.iret1.ne.0.and.iret2.ne.0) then
                 imin0=idiff
                 nlen=ii
              endif
              endif
            else 
              if (fri.le.0.10.and.asspt1i.lt.0.10) then
              if (idiff.ge.imax)then
                 imax=idiff
                 nlen=ii
              endif
              else
              if (idiff.le.imin0)then
                 imin0=idiff
                 nlen=ii
              endif
              endif
           endif
        enddo

       endif

       if (asspt1i.lt.0.05.and.minlen.le.50) then
        if (fri.le.0.05) then
         nlen=8
        elseif (fri.gt.0.05.and.fri.lt.0.10) then
         nlen=6
        elseif (fri.ge.0.16) then
            ist=8
            jj=14
            istep=2

            do ii=ist,jj,istep
                ipp=isum-((isum/ii)*ii)
                if ((ipp).gt.(ii/2)) ipp=ipp-(ii/2)
                if (ipp.le.imin0) then
                    imin0=ipp
                    nlen=ii
                endif
            enddo
        endif
       endif

       if (nlen.eq.0)nlen=8

       return
       end
cccc 

**** subroutine to store the original score.
**** needed for the one fragment scanning.

       subroutine cpscore (scorei,scorer)
       parameter (maxres=3000)
       dimension scorei(maxres,maxres),scorer(maxres,maxres)
       common /length/ nseq1,nseq2

       do ii=1,nseq1
        do jj=1,nseq2
            scorer(ii,jj)=scorei(ii,jj)
        enddo
       enddo

       return
       end
cccc

****  subroutine to store the original dpsc matrix
**** needed for the joining of the fragments

       subroutine cpdpsc (nfr1,nfr2,smtxi,smtxr)
       parameter (mcom=6000)
       dimension smtxi(mcom,mcom),smtxr(mcom,mcom)
        
       do ii=1,nfr1
        do jj=1,nfr2
            smtxr(ii,jj)=smtxi(ii,jj)
        enddo
       enddo

       return
       end
cccc 

**** subroutine to find the SEQUENCE alignment, given the alignment of fragments
**** needs: DPSC, JALI (alignment), NALI (No. of alignments, 
****
**** returns: TMmax: best TM-score, MAP (alignment of sequence)

      subroutine calbesttm (dx,dxs,nlen,nali,jali,nfr1,nfr2,dpsci,fri,mapr,b1TMmax)
      parameter (maxres=3000)
      parameter (maxfr=1000)    ! maximum no. of fragments
      parameter (nfraln=6000)   ! maximum no. of frag. combination used for
      parameter (mfrali=100)   ! max. TOTAL no. of alignments of FRAGMENTS

      dimension dpsci(nfraln,nfraln),jali(mfrali,maxfr),ifrali(maxfr)
      dimension invmap(maxres),invmap_i(maxres),gap(10),mapr(maxres)
      dimension score(maxres,maxres),scorei(maxres,maxres),invmapr(maxres)

      common /length/ nseq1,nseq2

      ngap=2
      gap(1)=-0.6
      gap(2)=0.0
      b1TMmax=0.0
        
      do 114 iter=1,nali
        do jj=1,nfr2
            ifrali(jj)=jali(iter,jj)
        enddo

           call fillinvmap(invmap)
           kk=0

           do jj=1,nfr2

            if (ifrali(jj).gt.0) then
                ist1=ifrali(jj)
                ist2=jj

                if (mod(ist1,2) .eq. 0) then 
                    istart1=(ist1-1)*(nlen/2)+1
                else
                    istart1=( ((ist1+1)/2) -1)*nlen+1
                endif
                if (mod(ist2,2) .eq. 0) then 
                    istart2=(ist2-1)*(nlen/2)+1
                else
                    istart2=( ((ist2+1)/2) -1)*nlen+1
                endif
                dpinf=dpsci(ist1,ist2)

                  if (kk.eq.0) then
                   itmpst=istart1
                   do ii=istart2,(istart2+nlen-1)
                     invmap(ii)=itmpst
                     itmpst=itmpst+1
                   enddo
                   kk=kk+1
                   istart1_old=istart1
                   istart2_old=istart2
                   dpold=dpinf
                  else
                    imove=0 ! accept second alignment
                    iol=0   ! second non overlapping
                    iend1=istart1_old+nlen-1
                    ires=istart2-istart2_old

                    if (dpinf.lt.dpold) then
                        imove=1
                    endif
                    if (ires.lt.nlen)iol=1

                    if (iend1.gt.istart1) then
                        if (imove.eq.0.and.iol.eq.0) then
                            iend2_old=istart2_old+nlen
                            itmp2=iend2_old-(nlen/2)

                            do ii=itmp2,(itmp2+(nlen/2)-1)
                                invmap(ii)=-1
                            enddo
                            itmpst=istart1
                            iend2=istart2+nlen-1
                            if (iend2.gt.nseq2) iend2=nseq2

                            do ii=istart2,iend2
                                if (itmpst.le.nseq1) then
                                    invmap(ii)=itmpst
                                    itmpst=itmpst+1
                                endif
                            enddo
                        else
                            if(imove.eq.0.and.iol.eq.1) then
                                itmpst1=istart1
                                itmpst2=istart2
                                iend2=itmpst2+nlen-1
                            elseif(imove.eq.1.and.iol.eq.1)then
                                itmpst1=istart1+nlen/2
                                itmpst2=istart2+nlen/2
                                iend2=itmpst2+nlen/2-1
                            else
                                itmpst1=istart1+nlen/2
                                itmpst2=istart2+nlen/2
                                iend2=itmpst2+nlen/2-1
                            endif

                            if (iend2.gt.nseq2) iend2=nseq2

                            do ii=itmpst2,iend2
                                if (itmpst1.le.nseq1) then
                                    invmap(ii)=itmpst1
                                    itmpst1=itmpst1+1
                                endif
                            enddo
                        endif
                    else
                        if (iol.eq.0) then
                            itmpst=istart1
                            iend2=istart2+nlen-1
                            if (iend2.gt.nseq2) iend2=nseq2

                            do ii=istart2,iend2
                                if (itmpst.le.nseq1) then
                                    invmap(ii)=itmpst
                                    itmpst=itmpst+1
                                endif
                            enddo
                        else
                            if (imove.eq.0) then
                                itmpst1=istart1
                                itmpst2=istart2
                                iend2=itmpst2+nlen-1
                            else
                                itmpst1=istart1+nlen/2
                                itmpst2=istart2+nlen/2
                                iend2=itmpst2+nlen/2-1
                            endif

                            if (iend2.gt.nseq2) iend2=nseq2

                            do ii=itmpst2,iend2
                                if (itmpst1.le.nseq1) then
                                    invmap(ii)=itmpst1
                                    itmpst1=itmpst1+1
                                endif
                            enddo
                        endif
                    endif
                    istart1_old=istart1
                    istart2_old=istart2
                    dpold=dpinf
                  endif
            endif
           enddo

            gap_open=0.0
            call get_score1 (dx,invmap,score)
            call dp(score,gap_open,invmap_i,dpout)
            call get_score (dx,dxs,invmap_i,score,TM)
            aTM=TM
              if (fri.le.0.20) then
                if (aTM.le.0.20) then
                    add=0.3     !change from 0.3
                else
                    add=0.2 ! change from 0.2
                endif
              endif
              if (fri.gt.0.20) then
                if (aTM.ge.0.20) then
                    add=0.5
                else
                    add=0.2
                endif
              endif
            call fixscore (add,score)

            call get_score (dx,dxs,invmap,scorei,aTM)
                if (fri.gt.0.25) then
                   if(aTM.lt.0.20) then
                        add=0.4
                   else
                        add=0.1
                   endif
                   call fixscore (add,scorei)
                else
                    if (aTM.lt.0.20) then
                         add=0.2    !change from 0.3
                    else
                         add=0.1
                    endif
                    call fixscore (add,scorei)
                endif

            do igap=1,2
                if (igap.eq.2) then 
                    call cpscore (scorei,score) ! copying scorei in score, since with dp score changes
                endif

                call dp (score,gap_open,invmap_i,dpout)
                call get_score (dx,dxs,invmap_i,score,TM)
                tTM=TM
                call make_iter(dx,dxs,score,15,0,invmap_i,tTM,invmapr) ! changing from 20 to 10
                call getbest (tTM,invmapr,b1TMmax,mapr)
            enddo      
 114  continue

      return
      end

cccc 

      subroutine assp (asspt1,asspt2,nlenrep) 
      parameter (maxres=3000)
      common /secstr/ isec(maxres),jsec(maxres)
      common /length/ nseq1,nseq2

      imin=min(nseq1,nseq2)
      imax=max(nseq1,nseq2)
      fr=float(imin)/float(imax)

      isspt1=0
      isspt2=0
      if (imin.eq.nseq1) then
        do jj=1,nseq1
            if (isec(jj).eq.2.or.isec(jj).eq.4) isspt1=isspt1+1
        enddo
        asspt1=float(isspt1)/float(nseq1)

        do jj=1,nseq2
            if (jsec(jj).eq.2.or.jsec(jj).eq.4) isspt2=isspt2+1
        enddo
        asspt2=float(isspt2)/float(nseq2)
      else
        do jj=1,nseq2
            if (jsec(jj).eq.2.or.jsec(jj).eq.4) isspt1=isspt1+1
        enddo
        asspt1=float(isspt1)/float(nseq2)

        do jj=1,nseq1
            if (isec(jj).eq.2.or.isec(jj).eq.4) isspt2=isspt2+1
        enddo
        asspt2=float(isspt2)/float(nseq1)
      endif

      nlenrep1a=0
      nlenrep1b=0
      nlenrep2a=0
      nlenrep2b=0
      ilena=0
      ilenb=0

      do jj=1,nseq1-1
        if (isec(jj).eq.2.and.isec(jj+1).eq.2) then
            ilena=ilena+1
        elseif (isec(jj).eq.2.and.isec(jj+1).ne.2) then
            ilena=ilena+1
            if (ilena.gt.nlenrep1a) nlenrep1a=ilena
            ilena=0
        endif

        if (isec(jj).eq.4.and.isec(jj+1).eq.4) then
            ilenb=ilenb+1
        elseif (isec(jj).eq.4.and.isec(jj+1).ne.4) then
            ilenb=ilenb+1
            if (ilenb.gt.nlenrep1b) nlenrep1b=ilenb
            ilenb=0
        endif
      enddo

      ilena=0
      ilenb=0

      do jj=1,nseq2-1
        if (jsec(jj).eq.2.and.jsec(jj+1).eq.2) then
            ilena=ilena+1
        elseif (jsec(jj).eq.2.and.jsec(jj+1).ne.2) then
            ilena=ilena+1
            if (ilena.gt.nlenrep2a) nlenrep2a=ilena
            ilena=0
        endif

        if (jsec(jj).eq.4.and.jsec(jj+1).eq.4) then
            ilenb=ilenb+1
        elseif (jsec(jj).eq.4.and.jsec(jj+1).ne.4) then
            ilenb=ilenb+1
            if (ilenb.gt.nlenrep2b) nlenrep2b=ilenb
            ilenb=0
        endif
      enddo
        
      nlenrep1=max(nlenrep1a,nlenrep1b)
      nlenrep2=max(nlenrep2a,nlenrep2b)

      if (min(nlenrep1,nlenrep2).ne.0) then
        if (asspt1.gt.0.40.and.asspt2.gt.0.40) then
            nlenrep=max(nlenrep1,nlenrep2)
        else
            nlenrep=min(nlenrep1,nlenrep2)
        endif
      else
        nlenrep=max(nlenrep1,nlenrep2)
      endif

      idec=21
      if (min(nseq1,nseq2).lt.100) idec=17
      
      do while (nlenrep.gt.idec)
        nlenrep=nlenrep/2
      enddo

      if (mod(nlenrep,2).ne.0) nlenrep=nlenrep-1
      
      if (nlenrep.le.6) nlenrep=8

        nlen=0
        imin0=100
        imax0=0
        imax1=0

        ist=nlenrep-6
        jj=nlenrep+8
        if (ist.le.6) ist=6
        if (jj.ge.22)jj=22

        if (imin.le.60) then
          jj=16
        elseif (imin.gt.60.and.imin.le.100) then
          jj=18
        endif
        
        if (imax.lt.100) then
          jj=14
        endif

        icou2=0
        icou0=0

        do ii=ist,jj,2
            ip=0
            iret1=(nseq1/ii)*ii
            iret2=(nseq2/ii)*ii
            idiff=(nseq1+nseq2)-(iret1+iret2)
           if ((nseq1-iret1).ge.(ii/2))then
                idiff=idiff-(ii/2)
                iret1=iret1+(ii/2)
                ip=ip+1
           endif
           if ((nseq2-iret2).ge.(ii/2)) then
                idiff=idiff-(ii/2)
                iret2=iret2+(ii/2)
                ip=ip+1
           endif
                iret1=nseq1-iret1
                iret2=nseq2-iret2
                 if (idiff.ge.imax0.and.(ip.eq.2)) then
                     imax0=idiff
                     nlenrep=ii
                 endif

                 if (idiff.lt.imin0.and.(ip.eq.2)) then
                    imin0=idiff
                    nlenp1=ii
                 endif

                 if (idiff.ge.imax1.and.(ip.eq.0)) then
                     imax1=idiff
                     nlenp2=ii
                 endif

                 if (ip.eq.2) icou2=icou2+1
                 if (ip.eq.0) icou0=icou0+1
        enddo
        
        if (icou2.le.2) then
         imina=100
         ist=8
         jj=18

         if (imin.lt.55) then
            ist=6
            jj=14
            if (fr.gt.0.10.and.fr.lt.0.30) ist=8
            if (fr.ge.0.35) ist=8    !changed from 0.35
         elseif (imin.gt.55.and.imin.lt.60) then
            jj=14
         elseif (imin.ge.60.and.imin.lt.100) then
            ist=10
         endif

        do ii=ist,jj,2
            ip=0
            iret1=(nseq1/ii)*ii
            iret2=(nseq2/ii)*ii
            idiff=(nseq1+nseq2)-(iret1+iret2)
           if ((nseq1-iret1).ge.(ii/2))then
                idiff=idiff-(ii/2)
                iret1=iret1+(ii/2)
                ip=ip+1
           endif
           if ((nseq2-iret2).ge.(ii/2)) then
                idiff=idiff-(ii/2)
                iret2=iret2+(ii/2)
                ip=ip+1
           endif

            if (idiff.le.imina.and.ip.eq.1) then
                imina=idiff
                nlenp3=ii
            endif
        enddo
        endif

        if ( (icou2.lt.2.and.icou0.lt.2).or.(icou2.eq.2.and.icou0.lt.2)) then
           nlenrep=nlenp3 
        elseif (icou2.lt.2.and.icou0.ge.2) then
          nlenrep=nlenp2
          if (imin.le.50.and.fr.gt.0.20) then
           if (nlenp3.ne.0.and.nlenp3.lt.nlenp2) nlenrep=nlenp3
          endif
        elseif (icou2.eq.2.and.icou0.ge.2) then
          if (fr.lt.0.08) then
           if (imax1.lt.imin0) then
            nlenrep=nlenp2
           else
            nlenrep=nlenp1
           endif
          else
            if (nlenp2.gt.nlenp1) then
                nlenrep=nlenp2
            else
                nlenrep=nlenp1
            endif
          endif
           if (imax1.eq.imin0) nlenrep=nlenp2
           if (imin.gt.100.and.nlenp3.ne.0) then
                nlenrep=nlenp3
           endif
        endif

 151  continue   

      end

      subroutine initialize (dpsci,dpsccpi,sspi,smtxi,nfrmax)
      parameter (nfraln=6000) 
      dimension dpsci(nfraln,nfraln),dpsccpi(nfraln,nfraln),sspi(nfraln,nfraln)
      dimension smtxi(nfraln,nfraln)
      
       nfrmax=nfrmax*2

        do ii=1,nfrmax
          do jj=1,nfrmax
           dpsci(ii,jj)=0.0
           dpsccpi(ii,jj)=0.0
           sspi(ii,jj)=0.0
           smtxi(ii,jj)=0.0
          enddo
        enddo

      return
      end
