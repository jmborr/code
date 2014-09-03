C F2KCLI : Fortran 200x Command Line Interface
C copyright Interactive Software Services Ltd. 2001
C For conditions of use see manual.txt
C
C Platform    : DOS (16-bit real mode)
C Compiler    : Lahey F77L
C               Microsoft Fortran 4.x/5.x
C               Prospero Fortran
C               Watcom Fortran 77/16
C To compile  : f77l   f2kcli.for   (Lahey)
C               fl /c  f2kcli.for   (MS)
C               profor f2kcli/h2/b1 (Prospero)
C               wfc    f2kcli.for   (Watcom)
C Implementer : Lawson B. Wakefield, I.S.S. Ltd.
C Date        : February 2001
C
      SUBROUTINE F2KSUBSTR(STRING,ISTART,IEND)
C
C Locate start and end of first sub-string in supplied string
C
C STRING = Character string to search
C ISTART = Start position of delimited string
C          Returned as zero if string is blank
C IEND   = End position of delimited string
C          Returned as zero if string is blank
C
      CHARACTER*(*) STRING
      INTEGER       ISTART,IEND
C
C Find start of sub-string
C
      DO 100 IPOS = 1,LEN(STRING)
        IF (STRING(IPOS:IPOS).NE.' ') THEN
            ISTART = IPOS
C
C Find the end of the sub-string
C
            IPOS2 = INDEX(STRING(IPOS:),' ')
            IF (IPOS2.EQ.0) THEN
                IEND = LEN(STRING)
            ELSE
                IEND = IPOS2 + IPOS - 2
            ENDIF
            RETURN
        ENDIF
  100 CONTINUE
C
C String is blank
C
      ISTART = 0
      IEND   = 0
      RETURN
      END
C
      INTEGER FUNCTION LEN_TRIM(STRING)
C
C Return actual length of supplied string,
C excluding trailing blanks or zero if blank.
C
C STRING = String to search
C
      CHARACTER*(*) STRING
C
      DO 100 IPOS = LEN(STRING),1,-1
        IF (STRING(IPOS:IPOS).NE.' ') THEN
            LEN_TRIM = IPOS
            RETURN
        ENDIF
  100 CONTINUE
C
C String is blank
C
      LEN_TRIM = 0
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
      CHARACTER*(*) COMMAND
      INTEGER       LENGTH
      INTEGER       STATUS
C
      CHARACTER*127 ARGSTR
      EQUIVALENCE  (ARGSTR,IARGSTR)
      LOGICAL       GETCMD
C
      SAVE          ARGSTR
C
      DATA          GETCMD/.TRUE./
C
      IF (GETCMD) THEN
          ARGSTR = ' '
          CALL F2KGETCL(IARGSTR)
          GETCMD = .FALSE.
      END IF
      COMMAND = ARGSTR
      LENGTH  = LEN_TRIM(COMMAND)
      STATUS  = 0
      RETURN
      END
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
      INTEGER       IPOS,ISTART,IEND,NARGS,IPOS1,LENARG,ISTAT
      CHARACTER*127 ARGSTR
C
      DATA          NARGS/-1/
C
      IF (NARGS.EQ.-1) THEN
C
C Get whole command line
C
          CALL GET_COMMAND(ARGSTR,LENARG,ISTAT)
C
C Count command line arguments
C
          NARGS = 0
          IPOS  = 1
  100     CALL F2KSUBSTR(ARGSTR(IPOS:),ISTART,IEND)
          IF (ISTART.GT.0) THEN
              ISTART = ISTART + IPOS - 1
              IEND   = IEND   + IPOS - 1
              IPOS   = IEND   + 2
C
C Is argument quoted ?
C
              IF (ARGSTR(ISTART:ISTART).NE.'"') THEN
C
C No - increment arg count
C
                  NARGS = NARGS + 1
              ELSE IF (ISTART.LT.LENARG) THEN
C
C Yes it is quoted and quote isn't at end of string
C
                  ISTART = ISTART + 1
                  CALL F2KSUBSTR(ARGSTR(ISTART:),IPOS1,IEND)
                  IF (IPOS1.GT.0) THEN
                      IEND = INDEX(ARGSTR(ISTART:),'"')
C
C Ignore null quotes
C
                      IF (IEND.NE.1) THEN
                          IF (IEND.EQ.0) THEN
                              IEND = LENARG
                          ELSE
                              IEND = ISTART + IEND - 2
                          ENDIF
                          NARGS = NARGS + 1
                          IPOS  = IEND  + 3
                      ENDIF
                  ENDIF
              ENDIF
C
              IF (IPOS.LE.LENARG) GOTO 100
          ENDIF
      ENDIF
C
      COMMAND_ARGUMENT_COUNT = NARGS
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
      INTEGER       IPOS,ISTART,IEND,NARGS,IPOS1,IENV,LENARG,ISTAT
      INTEGER*2     ICUR,IPREV
      CHARACTER*127 ARGSTR
C
C If the argument number is negative, return an error code of 1.
C
      IF (NUMBER.LT.0) THEN
          VALUE  = ' '
          LENGTH = 0
          STATUS = 1
          RETURN
      ENDIF
C
C Get executable name if requested
C
      IF (NUMBER.EQ.0) THEN
C
C Get address of environment
C
          CALL F2KENVAD(IENV)
C
C Skip past environment strings table (terminated by double null)
C
          ICUR  = -1
          IPREV = -1
   10     CALL F2KPEEK(IENV,ICUR)
          IENV  = IENV + 1
          IF (ICUR.NE.0.OR.IPREV.NE.0) THEN
              IPREV = ICUR
              GOTO 10
          ENDIF
C
C Skip past next two bytes after environment strings
C
          CALL F2KPEEK(IENV,ICUR)
          IENV  = IENV + 1
          CALL F2KPEEK(IENV,ICUR)
          IENV  = IENV + 1
C
C Collect null terminated executable name
C
          LENGTH = 0
          VALUE  = ' '
   20     CALL F2KPEEK(IENV,ICUR)
          IENV   = IENV + 1
          IF (ICUR.NE.0) THEN
              IF (LENGTH.EQ.LEN(VALUE)) THEN
                  STATUS = 2
                  RETURN
              ENDIF
              LENGTH = LENGTH + 1
              VALUE(LENGTH:LENGTH) = CHAR(ICUR)
              GOTO 20
          ENDIF
          STATUS = 0
          RETURN
      ENDIF
C
C Get whole command line
C
      CALL GET_COMMAND(ARGSTR,LENARG,ISTAT)
C
C Find required command line argument
C
      NARGS = 0
      IPOS  = 1
  100 CALL F2KSUBSTR(ARGSTR(IPOS:),ISTART,IEND)
      IF (ISTART.GT.0) THEN
          ISTART = ISTART + IPOS - 1
          IEND   = IEND   + IPOS - 1
          IPOS   = IEND   + 2
C
C Is argument quoted ?
C
          IF (ARGSTR(ISTART:ISTART).NE.'"') THEN
C
C No - increment arg count
C
              NARGS = NARGS + 1
          ELSE IF (ISTART.LT.LENARG) THEN
C
C Yes it is quoted and quote isn't at end of string
C
              ISTART = ISTART + 1
              CALL F2KSUBSTR(ARGSTR(ISTART:),IPOS1,IEND)
              IF (IPOS1.GT.0) THEN
                  IEND = INDEX(ARGSTR(ISTART:),'"')
C
C Ignore null quotes
C
                  IF (IEND.NE.1) THEN
                      IF (IEND.EQ.0) THEN
                          IEND = LENARG
                      ELSE
                          IEND = ISTART + IEND - 2
                      ENDIF
                      NARGS = NARGS + 1
                      IPOS  = IEND  + 3
                  ENDIF
              ENDIF
          ENDIF
C
C If this is the required command line argument, return value
C and exit otherwise continue if not at end of command line
C
          IF (NUMBER.EQ.NARGS) THEN
              VALUE  = ARGSTR(ISTART:IEND)
              LENGTH = IEND - ISTART + 1
              STATUS = 0
              RETURN
          ELSE IF (IPOS.LE.LENARG) THEN
              GOTO 100
          ENDIF
      ENDIF
C
C Error code = 2 : NUMBER too large
C
      VALUE  = ' '
      LENGTH = 0
      STATUS = 2
      RETURN
      END
