! F2KCLI : Fortran 200x Command Line Interface
! copyright Interactive Software Services Ltd. 2001
! For conditions of use see manual.txt
!
! Platform    : DOS (Phar Lap 32-bit DOS extender)
! Compiler    : Lahey LF90
! To compile  : lf90 -c f2kcli.f90
! Implementer : Lawson B. Wakefield, I.S.S. Ltd.
! Date        : February 2001
!
      MODULE F2KCLI
!
      CONTAINS
!
      SUBROUTINE F2KSUBSTR(STRING,ISTART,IEND)
!
! Locate start and end of first sub-string in supplied string
!
! STRING = Character string to search
! ISTART = Start position of delimited string
!          Returned as zero if string is blank
! IEND   = End position of delimited string
!          Returned as zero if string is blank
!
      CHARACTER(LEN=*), INTENT(IN)  :: STRING
      INTEGER         , INTENT(OUT) :: ISTART,IEND
!
      INTEGER :: IPOS,IPOS2
!
! Find start of sub-string
!
      DO IPOS = 1,LEN(STRING)
        IF (STRING(IPOS:IPOS) /= ' ') THEN
            ISTART = IPOS
!
! Find the end of the sub-string
!
            IPOS2 = INDEX(STRING(IPOS:),' ')
            IF (IPOS2 == 0) THEN
                IEND = LEN(STRING)
            ELSE
                IEND = IPOS2 + IPOS - 2
            END IF
            RETURN
        END IF
      END DO
!
! String was blank
!
      ISTART = 0
      IEND   = 0
      RETURN
      END SUBROUTINE F2KSUBSTR
!
      SUBROUTINE GET_COMMAND(COMMAND,LENGTH,STATUS)
!
! Description. Returns the entire command by which the program was
!   invoked.
!
! Class. Subroutine.
!
! Arguments.
! COMMAND (optional) shall be scalar and of type default character.
!   It is an INTENT(OUT) argument. It is assigned the entire command
!   by which the program was invoked. If the command cannot be
!   determined, COMMAND is assigned all blanks.
! LENGTH (optional) shall be scalar and of type default integer. It is
!   an INTENT(OUT) argument. It is assigned the significant length
!   of the command by which the program was invoked. The significant
!   length may include trailing blanks if the processor allows commands
!   with significant trailing blanks. This length does not consider any
!   possible truncation or padding in assigning the command to the
!   COMMAND argument; in fact the COMMAND argument need not even be
!   present. If the command length cannot be determined, a length of
!   0 is assigned.
! STATUS (optional) shall be scalar and of type default integer. It is
!   an INTENT(OUT) argument. It is assigned the value 0 if the
!   command retrieval is sucessful. It is assigned a processor-dependent
!   non-zero value if the command retrieval fails.
!
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: COMMAND
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
!
      CHARACTER(LEN=127), SAVE :: ARGSTR
      LOGICAL           , SAVE :: GETCMD = .TRUE.
!
      IF (GETCMD) THEN
          CALL GETCL(ARGSTR)
          GETCMD = .FALSE.
      END IF
!
      IF (PRESENT(COMMAND)) COMMAND = ARGSTR
      IF (PRESENT(LENGTH )) LENGTH  = LEN_TRIM(ARGSTR)
      IF (PRESENT(STATUS )) STATUS  = 0
      RETURN
      END SUBROUTINE GET_COMMAND
!
      INTEGER FUNCTION COMMAND_ARGUMENT_COUNT()
!
! Description. Returns the number of command arguments.
!
! Class. Inquiry function
!
! Arguments. None.
!
! Result Characteristics. Scalar default integer.
!
! Result Value. The result value is equal to the number of command
!   arguments available. If there are no command arguments available
!   or if the processor does not support command arguments, then
!   the result value is 0. If the processor has a concept of a command
!   name, the command name does not count as one of the command
!   arguments.
!
      INTEGER            :: IPOS,ISTART,IEND,IPOS1,LENARG
      CHARACTER(LEN=127) :: ARGSTR
      INTEGER, SAVE      :: NARGS = -1
!
      IF (NARGS == -1) THEN
!
! Get whole command line
!
          CALL GET_COMMAND(ARGSTR,LENARG)
!
! Count command line arguments
!
          NARGS = 0
          IPOS  = 1
  100     CALL F2KSUBSTR(ARGSTR(IPOS:),ISTART,IEND)
          IF (ISTART > 0) THEN
              ISTART = ISTART + IPOS - 1
              IEND   = IEND   + IPOS - 1
              IPOS   = IEND   + 2
!
! Is argument quoted ?
!
              IF (ARGSTR(ISTART:ISTART) /= '"') THEN
!
! No - increment arg count
!
                  NARGS = NARGS + 1
              ELSE IF (ISTART < LENARG) THEN
!
! Yes it is quoted and quote isn't at end of string
!
                  ISTART = ISTART + 1
                  CALL F2KSUBSTR(ARGSTR(ISTART:),IPOS1,IEND)
                  IF (IPOS1 > 0) THEN
                      IEND = INDEX(ARGSTR(ISTART:),'"')
!
! Ignore null quotes
!
                      IF (IEND /= 1) THEN
                          IF (IEND == 0) THEN
                              IEND = LENARG
                          ELSE
                              IEND = ISTART + IEND - 2
                          END IF
                          NARGS = NARGS + 1
                          IPOS  = IEND  + 3
                      END IF
                  END IF
              END IF
!
              IF (IPOS <= LENARG) GO TO 100
          END IF
      END IF
!
      COMMAND_ARGUMENT_COUNT = NARGS
      RETURN
      END FUNCTION COMMAND_ARGUMENT_COUNT
!
      SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
!
! Description. Returns a command argument.
!
! Class. Subroutine.
!
! Arguments.
! NUMBER shall be scalar and of type default integer. It is an
!   INTENT(IN) argument. It specifies the number of the command
!   argument that the other arguments give information about. Useful
!   values of NUMBER are those between 0 and the argument count
!   returned by the COMMAND_ARGUMENT_COUNT intrinsic.
!   Other values are allowed, but will result in error status return
!   (see below).  Command argument 0 is defined to be the command
!   name by which the program was invoked if the processor has such
!   a concept. It is allowed to call the GET_COMMAND_ARGUMENT
!   procedure for command argument number 0, even if the processor
!   does not define command names or other command arguments.
!   The remaining command arguments are numbered consecutively from
!   1 to the argument count in an order determined by the processor.
! VALUE (optional) shall be scalar and of type default character.
!   It is an INTENT(OUT) argument. It is assigned the value of the
!   command argument specified by NUMBER. If the command argument value
!   cannot be determined, VALUE is assigned all blanks.
! LENGTH (optional) shall be scalar and of type default integer.
!   It is an INTENT(OUT) argument. It is assigned the significant length
!   of the command argument specified by NUMBER. The significant
!   length may include trailing blanks if the processor allows command
!   arguments with significant trailing blanks. This length does not
!   consider any possible truncation or padding in assigning the
!   command argument value to the VALUE argument; in fact the
!   VALUE argument need not even be present. If the command
!   argument length cannot be determined, a length of 0 is assigned.
! STATUS (optional) shall be scalar and of type default integer.
!   It is an INTENT(OUT) argument. It is assigned the value 0 if
!   the argument retrieval is sucessful. It is assigned a
!   processor-dependent non-zero value if the argument retrieval fails.
!
! NOTE
!   One possible reason for failure is that NUMBER is negative or
!   greater than COMMAND_ARGUMENT_COUNT().
!
      INTEGER         , INTENT(IN)            :: NUMBER
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VALUE
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
!
      INTEGER            :: IPOS,ISTART,IEND,NARGS,IPOS1,LENARG,IENV
      CHARACTER(LEN=127) :: ARGSTR
      CHARACTER(LEN=260) :: EXENAME
      INTEGER*2             ICUR,IPREV
!
! If the argument number is negative, return an error code of 1.
!
      IF (NUMBER < 0) THEN
          IF (PRESENT(VALUE))  VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 1
          RETURN
      END IF
!
! Get executable name if requested
!
      IF (NUMBER == 0) THEN
!
! IENV = Address of environment in the Phar Lap DOS extender (2C0000h)
!
          IENV  = 2883584
!
! Skip past environment strings table (terminated by double null)
!
          ICUR  = -1
          IPREV = -1
   10     CALL F2KPEEK(IENV,ICUR)
          IENV  = IENV + 1
          IF (ICUR /= 0.OR.IPREV /= 0) THEN
              IPREV = ICUR
              GO TO 10
          END IF
!
! Skip past next two bytes after environment strings
!
          CALL F2KPEEK(IENV,ICUR)
          IENV  = IENV + 1
          CALL F2KPEEK(IENV,ICUR)
          IENV  = IENV + 1
!
! Collect null terminated executable name
!
          IPOS    = 0
          EXENAME = ' '
   20     CALL F2KPEEK(IENV,ICUR)
          IENV   = IENV + 1
          IF (ICUR /= 0) THEN
              IPOS = IPOS + 1
              IF (IPOS <= LEN(EXENAME)) EXENAME(IPOS:IPOS) = CHAR(ICUR)
              GO TO 20
          END IF
          IF (PRESENT(VALUE))  VALUE  = EXENAME
          IF (PRESENT(LENGTH)) LENGTH = IPOS
          IF (PRESENT(STATUS)) STATUS = 0
          RETURN
      ENDIF
!
! Get whole command line
!
      CALL GET_COMMAND(ARGSTR,LENARG)
!
! Find required command line argument
!
      NARGS = 0
      IPOS  = 1
  100 CALL F2KSUBSTR(ARGSTR(IPOS:),ISTART,IEND)
      IF (ISTART > 0) THEN
          ISTART = ISTART + IPOS - 1
          IEND   = IEND   + IPOS - 1
          IPOS   = IEND   + 2
!
! Is argument quoted ?
!
          IF (ARGSTR(ISTART:ISTART) /= '"') THEN
!
! No - increment arg count
!
              NARGS = NARGS + 1
          ELSE IF (ISTART < LENARG) THEN
!
! Yes it is quoted and quote isn't at end of string
!
              ISTART = ISTART + 1
              CALL F2KSUBSTR(ARGSTR(ISTART:),IPOS1,IEND)
              IF (IPOS1 > 0) THEN
                  IEND = INDEX(ARGSTR(ISTART:),'"')
!
! Ignore null quotes
!
                  IF (IEND /= 1) THEN
                      IF (IEND == 0) THEN
                          IEND = LENARG
                      ELSE
                          IEND = ISTART + IEND - 2
                      END IF
                      NARGS = NARGS + 1
                      IPOS  = IEND  + 3
                  END IF
              END IF
          END IF
!
! If this is the required command line argument, return value
! and exit otherwise continue if not at end of command line
!
          IF (NUMBER == NARGS) THEN
              IF (PRESENT(VALUE))  VALUE  = ARGSTR(ISTART:IEND)
              IF (PRESENT(LENGTH)) LENGTH = IEND - ISTART + 1
              IF (PRESENT(STATUS)) STATUS = 0
              RETURN
          ELSE IF (IPOS <= LENARG) THEN
              GO TO 100
          END IF
      END IF
!
! Error code = 2 : NUMBER too large
!
      IF (PRESENT(VALUE))  VALUE  = ' '
      IF (PRESENT(LENGTH)) LENGTH = 0
      IF (PRESENT(STATUS)) STATUS = 2
      RETURN
      END SUBROUTINE GET_COMMAND_ARGUMENT
!
      END MODULE F2KCLI
