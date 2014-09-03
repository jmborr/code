c     pgf77 -Mextend -O -s -fast -Wl,-static -o inputArgs_example.x inputArgs_example.f inputArgs.f
      program inputArgs_example
      integer one,two
      character*259 three
      real four
      character five*200
      call read_command_line(one,two,three,four,five)
      end

      subroutine read_command_line(one,two,three,four,five)
c     these are the variables to be initialized from command line
      integer one,two
      character*259 three
      real four
      character five*200
c     variables intrinsic to the subroutine
      character inputs*255(26)
      character wellcome*255(0:26)
      integer c2i
c     the wellcome message
      wellcome(0)='./junk.x [options]'
      wellcome(1)='-a _RA_one one argument'
      wellcome(2)='-b __two two is optional(def:3)'
      wellcome(3)='-d __three three is optional too(def:seq.txt)'
      wellcome(4)='-c _A_four four is optional too(def:1.0)'
      wellcome(5)='-e _R_five five is required'
c     init default variables
      two=3
      three='seq.txt'
      four=1.0
c     pass command values from the command line
      call proccess_command_line(wellcome,inputs)
      read(inputs(c2i('a')), *) one
      if(inputs(c2i('b'))(1:1).ne.'') read(inputs(c2i('b')), *) two
      if(inputs(c2i('c'))(1:1).ne.' ') read(inputs(c2i('c')), *) three
      if(inputs(c2i('d'))(1:1).ne.' ') read(inputs(c2i('d')),*) four
      read(inputs(c2i('e')), *) five
      end
