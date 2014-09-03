c     pgf77 -Mextend -O -s -fast -Wl,-static -o inputArgs.x inputArgs.f
c      program inputArgs
c      integer one,two
c      character*259 three
c      real four
c      character five*200
c      call read_command_line(one,two,three,four,five)


cc     ############################################################
cc     this is a sample subroutine, the only subroutine that has to be
cc     hard-coded for every program that you write
cc     read input arguments from command line
c      subroutine read_command_line(one,two,three,four,five)
cc     these are the variables to be initialized from command line
c      integer one,two
c      character*259 three
c      real four
c      character five*200
cc     variables intrinsic to the subroutine
c      character inputs*255(26)
c      character wellcome*255(0:26)
c      integer c2i
cc     the wellcome message
c      wellcome(0)='./junk.x [options]'
c      wellcome(1)='-a _RA_one one argument'
c      wellcome(2)='-b __two two is optional(def:3)'
c      wellcome(3)='-d __three three is optional too(def:seq.txt)'
c      wellcome(4)='-c _A_four four is optional too(def:1.0)'
c      wellcome(5)='-e _R_five five is required'
cc     init default variables
c      two=3
c      three='seq.txt'
c      four=1.0
cc     pass command values from the command line
c      call proccess_command_line(wellcome,inputs)
c      read(inputs(c2i('a')), *) one
c      if(inputs(c2i('b'))(1:1).ne.'') read(inputs(c2i('b')), *) two
c      if(inputs(c2i('c'))(1:1).ne.' ') read(inputs(c2i('c')), *) three
c      if(inputs(c2i('d'))(1:1).ne.' ') read(inputs(c2i('d')),*) four
c      read(inputs(c2i('e')), *) five
c      end


c     ############################################################
      subroutine proccess_command_line(wellcome,inputs)
      parameter(nargsmax=26)
      common/inputArgs/message,is_valid,nvalid,is_required,nreq
      character inputs*255(nargsmax),wellcome*255(0:nargsmax)
      character message*7140 !255*nargsmax
      integer is_valid(nargsmax),nvalid,is_required(nargsmax),nreq
      integer i,ireq,narg,iarg,length,istatus,option,c2i
      character onearg*255
      integer  command_argument_count
      external command_argument_count
      call parse_wellcome(wellcome)
      do i=1,nargsmax
         inputs(i)=''
      enddo
c      call print_message(message)
c      call exit(1)

      narg=command_argument_count() !number of passed arguments
c      write(6,*)narg
      if(narg.eq.0)then
         if(nreq.gt.0)  call print_message(message)
         return
      endif
      ireq=0
      iarg=0
 5    continue
      iarg=iarg+1
      call get_command_argument(iarg,onearg,length,istatus)
c      write(6,*)'onearg=',onearg(2:2)
      option=c2i(onearg(2:2))
c      write(6,*)option
c      call exit(1)
      if(is_valid(option).eq.0) call print_message(message)
      if(is_required(option).eq.1) ireq=ireq+1
      iarg=iarg+1
      call get_command_argument(iarg,onearg,length,istatus)
c      write(6,'a,a')'onearg=',onearg
c      call exit(1)
      read(onearg,*) inputs(option) !store argument in inputs
c      write(6,*)option,'inputs(option)=',inputs(option)
c      call exit(1)
c     check if we finished processing the arguments
      if(iarg.lt.narg) then
         goto 5
      endif
c      write(6,*)iarg,narg,ireq,nreq
c      call exit(1)
      if(ireq.lt.nreq)then
         call print_message(message)   !not enough number of required arguments
      endif
      end


c     ############################################################
      subroutine parse_wellcome(wellcome)
      parameter(nargsmax=26)
      common/inputArgs/message,is_valid,nvalid,is_required,nreq
      integer i,j,k,is_valid(nargsmax),nvalid,is_required(nargsmax),nreq
      integer size,msize,rsize,osize,optbegin,optend,we_require,c2i
      character message*7140 !255*nargsmax
      character line*255,letter,requiredline*7140,optionline*7140
      character wellcome*255(0:nargsmax)

c     initialized arrays is_valid and is_required
      do i=1,nargsmax
         is_valid(i)=0
         is_required(i)=0
      enddo
c     initialize message with wellcome(0), which should be the greeting
      line=wellcome(0)
      do i=255,1,-1
         if(line(i:i) .ne. ' ')then
            size=i
            goto 101
         endif
      enddo
 101  message='Usage: '//wellcome(0)(1:size)//'\n'
      msize=size+8 !current size of message
      rsize=0 !size of message related to required arguments
      osize=0 !size of message related to optional arguments
      do i=1,nargsmax
         line=wellcome(i)
         do j=255,1,-1
            if(line(j:j) .ne. ' ')then
               size=j
               goto 111
            endif
         enddo
 111     if(size.eq.255)then    !no more input lines
            goto 121
         endif
         k=c2i(line(2:2))  !1 for a, 2 for b, 3 for c,...
         is_valid(k)=1
         if( we_require(line,optbegin,optend).eq.1 )then
            is_required(k)=1
            if(rsize.eq.0)then
               requiredline='  '//line(1:optbegin-1)//
     &              line(1+optend:size)//'\n'
               rsize=2+optbegin+size-optend
            else
               requiredline=requiredline(1:rsize)//'  '//
     &              line(1:optbegin-1)//line(1+optend:size)
     &              //'\n'
               rsize=2+rsize+optbegin+size-optend
            endif
            nreq=nreq+1
         else
            if(osize.eq.0)then
               optionline='  '//line(1:optbegin-1)//
     &              line(1+optend:size)//'\n'
               osize=2+optbegin+size-optend
            else
               optionline=optionline(1:osize)//'  '//
     &              line(1:optbegin-1)//line(1+optend:size)/
     &              /'\n'
               osize=2+osize+optbegin+size-optend
            endif   
         endif
         nvalid=nvalid+1
      enddo
 121  if(rsize.gt.0)then
         message=message(1:msize)//' Required arguments:\n'//
     &     requiredline(1:rsize)
         msize=msize+21+rsize
      endif
      if(osize.gt.0)then
         message=message(1:msize)//' Optional arguments:\n'//
     &     optionline(1:osize)
         msize=msize+21+osize
      endif
      end


c     ##########################################################
      function we_require(line,optbegin,optend)
      integer we_require,i,j,optbegin,optend
      character line*255
      we_require=0 !init as not required
      optbegin=0
      do i=1,255
         if(line(i:i).eq.'_')then
            if(optbegin.eq.0) then
               optbegin=i
            else
               optend=i
               goto 101
            endif
         endif
      enddo
 101  if(optend.eq.1+optbegin) goto 111
      do i=optbegin,optend
         if(line(i:i).eq.'R') we_require=1
      enddo
 111  do i=optend,255
         if(line(i:i).eq.' ')then
            optend=i
            goto 121
         endif
      enddo
 121  end 


c     ##########################################################
c     print message and exit
      subroutine print_message(message)
      character message*7140
      integer i,size
      size=1
      do i=7140,1,-1
         if(message(i:i).ne.' ')then
            size=i
            goto 101
         endif
      enddo
 101  call system('clear')
      write(6,'a') message(1:size)
      call exit(1)
      end


c     ##########################################################
      function c2i(letter)
      character letter
      integer c2i
      c2i=1+ichar(letter)-ichar('a')
      end

C F2KCLI : Fortran 200x Command Line Interface
C copyright Interactive Software Services Ltd. 2001
C For conditions of use see manual.txt
C
C Platform    : Linux/Solaris
C Compiler    : PGI Fortran 77
C To compile  : pgf77 -c f2kcli.f
C Implementer : Lawson B. Wakefield, I.S.S. Ltd.
C Date        : October 2005
C
      INTEGER FUNCTION COMMAND_ARGUMENT_COUNT()
C
C Description. Returns the number of command arguments.
C
C Class. Inquiry function
C
C Arguments. None.
C
C Result Characteristics. Scalar default integer.
C
C Result Value. The result value is equal to the number of command
C   arguments available. If there are no command arguments available
C   or if the processor does not support command arguments, then
C   the result value is 0. If the processor has a concept of a command
C   name, the command name does not count as one of the command
C   arguments.
C
      INTEGER  IARGC
      EXTERNAL IARGC
C
      COMMAND_ARGUMENT_COUNT = IARGC()
      RETURN
      END
C
      SUBROUTINE GET_COMMAND(COMMAND,LENGTH,STATUS)
C
C Description. Returns the entire command by which the program was
C   invoked.
C
C Class. Subroutine.
C
C Arguments.
C COMMAND (optional) shall be scalar and of type default character.
C   It is an INTENT(OUT) argument. It is assigned the entire command
C   by which the program was invoked. If the command cannot be
C   determined, COMMAND is assigned all blanks.
C LENGTH (optional) shall be scalar and of type default integer. It is
C   an INTENT(OUT) argument. It is assigned the significant length
C   of the command by which the program was invoked. The significant
C   length may include trailing blanks if the processor allows commands
C   with significant trailing blanks. This length does not consider any
C   possible truncation or padding in assigning the command to the
C   COMMAND argument; in fact the COMMAND argument need not even be
C   present. If the command length cannot be determined, a length of
C   0 is assigned.
C STATUS (optional) shall be scalar and of type default integer. It is
C   an INTENT(OUT) argument. It is assigned the value 0 if the
C   command retrieval is sucessful. It is assigned a processor-dependent
C   non-zero value if the command retrieval fails.
C
C NOTE
C (1) The Fortran 77 implementation of this routine does not support
C     optional arguments, so all arguments must be specified by the
C     caller.
C
      CHARACTER*(*)  COMMAND
      INTEGER        LENGTH
      INTEGER        STATUS
C
      INTEGER        IARG,NARG,IPOS,LENARG
      CHARACTER*2000 ARGSTR
      LOGICAL        GETCMD
C
      SAVE           ARGSTR
      SAVE           LENARG
C
      INTEGER        IARGC
      EXTERNAL       IARGC
C
      DATA           GETCMD/.TRUE./
C
C Under Unix we must reconstruct the command line from its constituent
C parts. This will not be the original command line. Rather it will be
C the expanded command line as generated by the shell.
C
      IF (GETCMD) THEN
          NARG = IARGC()
          IF (NARG.GT.0) THEN
              IPOS = 1
              DO 100 IARG = 1,NARG
                CALL GETARG(IARG,ARGSTR(IPOS:))
                LENARG = LEN_TRIMF2K(ARGSTR)
                IPOS   = LENARG + 2
                IF (IPOS.GT.LEN(ARGSTR)) GOTO 200
  100         CONTINUE
          ELSE
              ARGSTR = ' '
              LENARG = 0
          ENDIF
  200     GETCMD = .FALSE.
      ENDIF
      COMMAND = ARGSTR
      LENGTH  = LENARG
      STATUS  = 0
      RETURN
      END
C
      SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
C
C Description. Returns a command argument.
C
C Class. Subroutine.
C
C Arguments.
C NUMBER shall be scalar and of type default integer. It is an
C   INTENT(IN) argument. It specifies the number of the command
C   argument that the other arguments give information about. Useful
C   values of NUMBER are those between 0 and the argument count
C   returned by the COMMAND_ARGUMENT_COUNT intrinsic.
C   Other values are allowed, but will result in error status return
C   (see below).  Command argument 0 is defined to be the command
C   name by which the program was invoked if the processor has such
C   a concept. It is allowed to call the GET_COMMAND_ARGUMENT
C   procedure for command argument number 0, even if the processor
C   does not define command names or other command arguments.
C   The remaining command arguments are numbered consecutively from
C   1 to the argument count in an order determined by the processor.
C VALUE (optional) shall be scalar and of type default character.
C   It is an INTENT(OUT) argument. It is assigned the value of the
C   command argument specified by NUMBER. If the command argument value
C   cannot be determined, VALUE is assigned all blanks.
C LENGTH (optional) shall be scalar and of type default integer.
C   It is an INTENT(OUT) argument. It is assigned the significant length
C   of the command argument specified by NUMBER. The significant
C   length may include trailing blanks if the processor allows command
C   arguments with significant trailing blanks. This length does not
C   consider any possible truncation or padding in assigning the
C   command argument value to the VALUE argument; in fact the
C   VALUE argument need not even be present. If the command
C   argument length cannot be determined, a length of 0 is assigned.
C STATUS (optional) shall be scalar and of type default integer.
C   It is an INTENT(OUT) argument. It is assigned the value 0 if
C   the argument retrieval is sucessful. It is assigned a
C   processor-dependent non-zero value if the argument retrieval fails.
C
C NOTE
C (1) One possible reason for failure is that NUMBER is negative or
C     greater than COMMAND_ARGUMENT_COUNT().
C (2) The Fortran 77 implementation of this routine does not support
C     optional arguments, so all arguments must be specified by the
C     caller.
C
      INTEGER       NUMBER
      CHARACTER*(*) VALUE
      INTEGER       LENGTH
      INTEGER       STATUS
C
      INTEGER       IARGC
      EXTERNAL      IARGC
C
C Possible error codes:
C 1 = Argument number is less than minimum
C 2 = Argument number exceeds maximum
C
      IF (NUMBER.LT.0) THEN
          VALUE  = ' '
          LENGTH = 0
          STATUS = 1
          RETURN
      ELSE IF (NUMBER.GT.IARGC()) THEN
          VALUE  = ' '
          LENGTH = 0
          STATUS = 2
          RETURN
      ENDIF
C
C Get the argument
C
      CALL GETARG(NUMBER,VALUE)
C
C The LENGTH option is fairly pointless under Unix.
C Trailing spaces can only be specified using quotes.
C Since the command line has already been processed by the
C shell before the application sees it, we have no way of
C knowing the true length of any quoted arguments.
C Just find last non-blank character in string.
C
      LENGTH = LEN_TRIMF2K(VALUE)
C
C Since GETARG does not return a result code, assume success
C
      STATUS = 0
      RETURN
      END
C
      INTEGER FUNCTION LEN_TRIMF2K(STRING)
C
C Return actual length of supplied string,
C excluding trailing blanks or zero if blank.
C
C STRING = String to search
C
      CHARACTER*(*) STRING
C
      INTEGER       IPOS
C
      DO 100 IPOS = LEN(STRING),1,-1
        IF (STRING(IPOS:IPOS).NE.' ') THEN
            LEN_TRIMF2K = IPOS
            RETURN
        ENDIF
  100 CONTINUE
C
C String is blank
C
      LEN_TRIMF2K = 0
      RETURN
      END
