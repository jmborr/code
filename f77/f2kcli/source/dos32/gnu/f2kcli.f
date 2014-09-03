C F2KCLI : Fortran 200x Command Line Interface
C copyright Interactive Software Services Ltd. 2001
C For conditions of use see manual.txt
C
C Platform    : DOS (EMX 32-bit DOS extender)
C Compiler    : GNU g77
C To compile  : g77 -c f2kcli.f
C Implementer : Lawson B. Wakefield, I.S.S. Ltd.
C Date        : February 2001
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
      CHARACTER*(*) COMMAND
      INTEGER       LENGTH
      INTEGER       STATUS
C
      INTEGER       IARG,NARG,IPOS,LENARG
      CHARACTER*127 ARGSTR
      LOGICAL       GETCMD
C
      SAVE          ARGSTR
      SAVE          LENARG
C
      DATA          GETCMD/.TRUE./
C
C g77 provides Unix-style command line access, so we must
C reconstruct the command line from its constituent parts.
C
      IF (GETCMD) THEN
          NARG = IARGC()
          IF (NARG.GT.0) THEN
              IPOS = 1
              DO 100 IARG = 1,NARG
                CALL GETARG(IARG,ARGSTR(IPOS:))
                LENARG = LEN_TRIM(ARGSTR)
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
      INTEGER       IPOS
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
      LENGTH = LEN_TRIMF2K(VALUE)
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