c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c     Read FASTA search output
c     2005_02_03
c
c***********************************************************************

c
c PARAMETER DECLARATION
c
      integer    inlinelen
      parameter (inlinelen =    50)
      integer    outlinelen
      parameter (outlinelen =   50)
      integer    mxseqlen
      parameter (mxseqlen =  20000)
      integer    mxseqn
      parameter (mxseqn =     2001)
      integer    mxfnlen
      parameter (mxfnlen =     512)
      integer    mxseqnamelen
      parameter (mxseqnamelen = 40)


c
c VARIABLE DECLARATION
c
      logical fileok
      logical dok(mxseqn)
      logical cok(mxseqn)

      integer iunit01, iunit02
      integer ounit01, ounit02
      integer spcpos1, slhpos1
      integer spcpos2, slhpos2
      integer spcpos3
      integer spcposq
      integer i, j, k
      integer iline
      integer iq, qlen
      integer isn1, isn2, isnq
      integer sq_len1, sq_len2
      integer al_start1, al_start2
      integer al_stop1, al_stop2
      integer al_display_start1, al_display_start2
      integer iseq1, iseq2
      integer fseq1, fseq2
      integer start
      integer albegin1, alend1
      integer albegin2, alend2
      integer aligned, match
      integer is, inogap, ns
      integer isn(mxseqn)
      integer idok, icok
      integer nlin

      real maxinter
      real mindist
      real minclose
      real maxclose
      real si(mxseqn,mxseqn)

      character*003 version
      character*006 ifn01
      character*003 ofn01
      character*012 ofn02
      character*002 xn, xl

      character*(mxfnlen) a(6)
      character*(mxfnlen) pname1, pname2, header
      character*(inlinelen) line
      character*(mxseqnamelen) sn1, sn2, snq
      character*(mxseqlen) seq1, seq2, seqline
      character*(mxseqlen)     s (mxseqn) 
      character*(mxseqnamelen) sn(mxseqn)

c
c DATA INITIALIZATION
c

c
c VARIABLE INITIALIZATION
c
      version = '002'
      ifn01   = 'infile'
      ofn01   = 'aln'
      ofn02   = 'fastaout.log'
      iunit01 = 10
      iunit02 = 11
      ounit01 = 12
      ounit02 = 13
      write(xn,'I2') mxseqnamelen
      write(xl,'I2') outlinelen

c
c MAIN
c

c Read parameters
      if (iargc() .ne. 7) then
       write(unit=6,*)'Enter the following 7 arguments:'
       write(unit=6,*)'1) path/name of the input sequence file'
       write(unit=6,*)'2) path/name of the output file'
       write(unit=6,*)'3) header to prepend to the output file names'
       write(unit=6,*)'4) max. seq. ident. for redundancy filter'
       write(unit=6,*)'5) min. seq. ident. to parent for distant set'
       write(unit=6,*)'6) min. seq. ident. to parent for close set'
       write(unit=6,*)'7) max. seq. ident. to parent for both sets'
       stop
      end if
      do i=1, 7
       call getarg(i,a(i))
      end do
c Input sequence file
      read (a(1), fmt='A') pname1
c Alignment output file
      read (a(2), fmt='A') pname2
c Header to prepend files
      read (a(3), fmt='A') header
c Read sequence identity parameters 
      read (a(4), fmt='F') maxinter
      read (a(5), fmt='F') mindist
      read (a(6), fmt='F') minclose
      read (a(7), fmt='F') maxclose
      close (unit=iunit01)

c Open files
      open  (unit=ounit02, file=ofn02, status='UNKNOWN')

c Check presence of name_input file
c     inquire (file=ifn01,exist=fileok)
c     if (.not. fileok) then
c      write (unit=ounit02,fmt='A,A,A')
c    $  'File ',ifn01,' does not exist'
c      stop
c     end if
c Read name_input file
c     open  (unit=iunit01, file=ifn01, status='OLD')
c Input sequence file
c     read (unit=iunit01, fmt='A', end=5) pname1
c Alignment output file
c     read (unit=iunit01, fmt='A', end=5) pname2
c Read sequence identity parameters 
c     read (unit=iunit01, fmt='F', end=5) maxinter
c     read (unit=iunit01, fmt='F', end=5) mindist
c     read (unit=iunit01, fmt='F', end=5) minclose
c     read (unit=iunit01, fmt='F', end=5) maxclose
c     close (unit=iunit01)
c End parameter reading
c     go to 7
c Bad input file format
c5     continue
c     write (unit=ounit02, fmt='A,A,A') 
c    $ 'File ',ifn01,' has missing parameter lines'
c     stop

c7     continue
c Detect first space in path&file string [file1]
      spcpos1=0
      do i=mxfnlen, 1, -1
       if (pname1(i:i) .eq. ' ') spcpos1=i
      end do
c Detect last slash in path&file string  [file1]
      slhpos1=0
      do i=1, mxfnlen
       if (pname1(i:i) .eq. '/') slhpos1=i
      end do
c Check presence of the file             [file1]
      inquire (file=pname1(1:spcpos1-1),exist=fileok)
      if (.not. fileok) then
       write(unit=ounit02,'A,A,A') 
     $  'File ',pname1(1:spcpos1-1),' does not exist'
       stop
      end if
c Detect first space in path&file string [file2]
      spcpos2=0
      do i=mxfnlen, 1, -1
       if (pname2(i:i) .eq. ' ') spcpos2=i
      end do
c Detect last slash in path&file string [file2]
      slhpos2=0
      do i=1, mxfnlen
       if (pname2(i:i) .eq. '/') slhpos2=i
      end do
c Check presence of the file            [file2]
      inquire (file=pname2(1:spcpos2-1),exist=fileok)
      if (.not. fileok) then
       write(unit=ounit02,'A,A,A') 
     $  'File ',pname2(1:spcpos2-1),' does not exist'
       stop
      end if
c Detect first space in header string
      spcpos3=0
      do i=mxfnlen, 1, -1
       if (header(i:i) .eq. ' ') spcpos3=i
      end do

c
c Read from input sequence file         [file1]
c
      write (unit=ounit02,fmt='A,A') 
     $ 'Reading: ',pname1(1:slhpos1)//pname1(slhpos1+1:spcpos1-1)
      open   (unit=iunit01,
     $        file=pname1(1:spcpos1-1),
     $        status='OLD')
      read  (unit=iunit01, fmt='A') seqline
      if (seqline(1:1) .ne. '>') then
       write (unit=ounit02,fmt='A,A') 
     $  'Input sequence in file: ',pname1(1:spcpos1-1),
     $  ' does not have FASTA format'
       stop
      end if
c Detect first space
      spcposq=0
      do i=mxseqnamelen, 1, -1
       if (seqline(i:i) .eq. ' ') spcposq=i
      end do
      isnq=spcposq-1      
c Extract fasta title      
      snq(2:isnq)=seqline(2:isnq)
       sn(1)=snq
      isn(1)=isnq
      write (unit=ounit02,fmt='A,A') 'Input seq. title : ', snq(2:isnq)
c Read input sequence
      iq=0
10    continue
      read  (unit=iunit01, fmt='A', end=50) seqline
      do i=1, mxseqlen
       if ( lge(seqline(i:i),'A') .and. lle(seqline(i:i),'Z') ) then
        iq=iq+1
        s(1)(iq:iq)=seqline(i:i)
       end if
      end do
      go to 10
50    continue
      qlen=iq
      write (unit=ounit02,fmt='A,I') 'Input seq. length: ', qlen

c                                                          
c Read from alignment output file       [file2]
c
      is=1
      write (unit=ounit02,fmt='A,A') 
     $ 'Reading: ',pname2(1:slhpos2)//pname2(slhpos2+1:spcpos2-1)
      open   (unit=iunit02,
     $        file=pname2(1:spcpos2-1),
     $        status='OLD')
      iline=0
100   continue
      read  (unit=iunit02, fmt='A', end=200) line
      iline=iline+1
      if (line(1:2) .eq. '>>' .and. line(3:3) .ne. '>') then
c Alignment record starts [ >> ]
       sq_len1=0
       al_start1=0
       al_stop1=0
       al_display_start1=0
       sq_len2=0
       al_start2=0
       al_stop2=0
       al_display_start2=0
110    continue

c      write (6,*) sq_len1, al_start1, al_stop1, al_display_start1
c      write (6,*) sq_len2, al_start2, al_stop2, al_display_start2

       read  (unit=iunit02, fmt='A', end=200) line
       iline=iline+1
c Ignore parameter lines
       if (line(1:1) .eq. ';') go to 110
       if (line(1:1) .ne. '>') then
        write (unit=ounit02,fmt='A,A') 
     $   'Unexpected format in line: ', iline
       end if
c Query sequence starts [ > ]
       iseq1=0
       sn1=line(1:mxseqnamelen)
       isn1=1
       do i=2, mxseqnamelen-2
        if (sn1(i:i+2) .eq. ' ..') then
         isn1=i-1
        end if
       end do
120    continue
       read  (unit=iunit02, fmt='A', end=200) line
       iline=iline+1
c Read parameter lines
       if (line(1:1) .eq. ';') then
        if (line(1:10) .eq. '; sq_len: ') then
         read (line(11:),'I') sq_len1
        end if
        if (line(1:12) .eq. '; al_start: ') then
         read (line(13:),'I') al_start1
        end if
        if (line(1:11) .eq. '; al_stop: ') then
         read (line(12:),'I') al_stop1
        end if
        if (line(1:20) .eq. '; al_display_start: ') then
         read (line(21:),'I') al_display_start1
        end if
        go to 120
       end if
       if (line(1:1) .ne. '>') then
c Query sequence line
        do i=1, inlinelen 
         iseq1=iseq1+1
         seq1(iseq1:iseq1)=line(i:i)
        end do
        go to 120
       end if
c Remove initial gaps
       start=iseq1
       do i=1, iseq1       
        if (seq1(i:i) .ne. '-') then
         start=i
         go to 125
        end if
       end do
125    continue
       do i=start, iseq1
        seq1(i-start+1:i-start+1)=seq1(i:i)
       end do
       iseq1=iseq1-start+1
       fseq1=iseq1
c Find alignment indexes
       albegin1=al_start1-al_display_start1+1
       j=al_start1-1
       do i=albegin1, fseq1 
        if (lge(seq1(i:i),'A') .and. lle(seq1(i:i),'Z')) j=j+1
        if (j .eq. al_stop1) then
         alend1=i
         go to 128
        end if
       end do
128    continue
c Library sequence starts [ > ]
       iseq2=0
       sn2=line(1:mxseqnamelen)
       isn2=1
       do i=2, mxseqnamelen-2
        if (sn2(i:i+2) .eq. ' ..') then
         isn2=i-1
        end if
       end do
130    continue
       read  (unit=iunit02, fmt='A', end=200) line
       iline=iline+1
c Read parameter lines
       if (line(1:1) .eq. ';') then
        if (line(1:10) .eq. '; sq_len: ') then
         read (line(11:),'I') sq_len2
        end if
        if (line(1:12) .eq. '; al_start: ') then
         read (line(13:),'I') al_start2
        end if
        if (line(1:11) .eq. '; al_stop: ') then
         read (line(12:),'I') al_stop2
        end if
        if (line(1:20) .eq. '; al_display_start: ') then
         read (line(21:),'I') al_display_start2
        end if
        go to 130
       end if
       if (line(1:1) .ne. '>') then
c Library sequence line
        do i=1, inlinelen 
         iseq2=iseq2+1
         seq2(iseq2:iseq2)=line(i:i)
        end do
        go to 130
       end if
       fseq2=iseq2
c Remove initial gaps
       start=iseq2
       do i=1, iseq2       
        if (seq2(i:i) .ne. '-') then
         start=i
         go to 135
        end if
       end do
135    continue
       do i=start, iseq2
        seq2(i-start+1:i-start+1)=seq2(i:i)
       end do
       iseq2=iseq2-start+1
       fseq2=fseq2-start+1
c Find alignment indexes
       albegin2=al_start2-al_display_start2+1
       j=al_start2-1
       do i=albegin2, fseq2
        if (lge(seq2(i:i),'A') .and. lle(seq2(i:i),'Z')) j=j+1
        if (j .eq. al_stop2) then
         alend2=i
         go to 138
        end if
       end do
138    continue
c Verify consistency of the alignments
       aligned=alend1-albegin1+1
       if (alend2-albegin2+1 .ne. aligned) then
        write (unit=ounit02,fmt='3(1XA)') 
     $  'Unexpected alignment format: ', sn1(2:isn1), sn2(2:isn2)
        stop
       end if
c Save aligned sequence and name
       is=is+1
c Name
        sn(is)=sn2
       isn(is)=isn2
c Sequence
       do i=1,qlen
        s(is)(i:i)='-'
       end do
       inogap=0
       do i=1, aligned
        if (lge(seq1(albegin1+i-1:albegin1+i-1),'A') .and.
     $      lle(seq1(albegin1+i-1:albegin1+i-1),'Z')      ) then
         inogap=inogap+1
         k=inogap-1+al_start1
         s(is)(k:k)=seq2(albegin2+i-1:albegin2+i-1)
        end if
       end do
c Sequences output
c      write (unit=ounit02,fmt='A,A') 
c    $  sn1, s(1)(1:qlen)
c      write (unit=ounit02,fmt='A,A') 
c    $  sn2, s(is)(1:qlen)

c Check end of aligments [ >>><<< ]
       if (line(1:6) .ne. '>>><<<') go to 110 
      end if

      go to 100

200   continue

c     Total number of sequences (is=1 is parent sequence)
      ns = is

c Sequence identity matrix
      do i=1, ns
       dok(i)=.true.
       do j=1, ns
        match=0
        do k=1, qlen
         if (s(i)(k:k) .eq. s(j)(k:k)) match=match+1
        end do
        si(i,j)=(1.*match)/(1.*qlen)
c       write (unit=ounit02,fmt='I4,1X,I4,1X,A'//xn//',1X,F5.3') 
c    $   i, j, sn(j)(2:isn(j)), si(i,j)
       end do
      end do

c Selection of sequences
      do i=1, ns
       dok(i)=.true.
       cok(i)=.true.
      end do
c No sequence less similar than mindist to parent in distant set
      do j=2, ns
       if (dok(j)) then
        if (si(1,j) .lt. mindist) then
         dok(j)=.false.
         cok(j)=.false.
        end if
       end if
      end do
c No sequence more similar than maxclose to parent in both sets
      do j=2, ns
       if (dok(j)) then
        if (si(1,j) .gt. maxclose) then
         dok(j)=.false.
         cok(j)=.false.
        end if
       end if
      end do
c No pair more similar than maxinter in both sets
      do i=1, ns-1
       if (dok(i)) then
        do j=i+1, ns
         if (dok(j)) then
          if (si(i,j) .gt. maxinter) then
           dok(j)=.false.
           cok(j)=.false.
          end if
         end if
        end do
       end if 
      end do
c No sequence less similar than minclose to parent in close set
      do j=2, ns
       if (dok(j)) then
        if (si(1,j) .lt. minclose) cok(j)=.false.
       end if
      end do
c Number of sequences is alignments
      idok=0
      icok=0
      do i=1, ns
       if (dok(i)) idok=idok+1
       if (cok(i)) icok=icok+1
      end do

c Output alignments
      nlin=1+int((qlen-1)/outlinelen)

c Distant alignment
      open  (unit=ounit01, 
     $       file=header(1:spcpos3-1)//'.dist.'//ofn01, 
     $       status='UNKNOWN')
      write (unit=ounit01, fmt='I4,1X,I5') idok, qlen
      do i=1, ns
       if (dok(i)) then
        if (nlin .eq. 1) then
         write  (unit=ounit01, fmt='A') s(i)(1:qlen)
                         else
         do j=1, nlin-1
          write (unit=ounit01, fmt='A') 
     $     s(i)(1+outlinelen*(j-1):outlinelen*j)
         end do
         write (unit=ounit01, fmt='A') 
     $    s(i)(1+outlinelen*(nlin-1):qlen)
        end if
       end if
      end do
      close  (unit=ounit01)

c Close alignment
      open  (unit=ounit01, 
     $       file=header(1:spcpos3-1)//'.close.'//ofn01, 
     $       status='UNKNOWN')
      write (unit=ounit01, fmt='I4,1X,I5') icok, qlen
      do i=1, ns
       if (cok(i)) then
        if (nlin .eq. 1) then
         write  (unit=ounit01, fmt='A') s(i)(1:qlen)
                         else
         do j=1, nlin-1
          write (unit=ounit01, fmt='A') 
     $     s(i)(1+outlinelen*(j-1):outlinelen*j)
         end do
         write (unit=ounit01, fmt='A') 
     $    s(i)(1+outlinelen*(nlin-1):qlen)
        end if
       end if
      end do

c Verify
      idok=0
      icok=0
      do i=1, ns
       if (dok(i)) idok=idok+1
       if (cok(i)) icok=icok+1
       write (unit=ounit02,
     $  fmt='I4,1X,A'//xn//',1X,F5.3,2(1X,L1,1X,I4)') 
     $  i,sn(i)(2:isn(i)), si(1,i), dok(i), idok, cok(i), icok
      end do

      close  (unit=ounit02)


      end
