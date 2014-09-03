! F2KCLI : Fortran 200x Command Line Interface
! copyright Interactive Software Services Ltd. 2001-2003
! For conditions of use see manual.txt
!
! Platform    : Win32
! Compiler    : F
! To compile  : f -c f2kcli.f90
! Implementer : Lawson B. Wakefield, I.S.S. Ltd.
! Date        : February 2003
!
      MODULE F2KCLI
!
      PRIVATE :: F2KSUBSTR
      PRIVATE :: f2kgetcl
      PRIVATE :: f2kgetexe
!
      PUBLIC  :: GET_COMMAND
      PUBLIC  :: COMMAND_ARGUMENT_COUNT
      PUBLIC  :: GET_COMMAND_ARGUMENT
!
! Interface definitions for C binding to GetCommandLine and GetModuleFilename
!
      INTERFACE
        PURE SUBROUTINE f2kgetcl(STR)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(OUT) :: STR
        END SUBROUTINE f2kgetcl
        PURE SUBROUTINE f2kgetexe(STR)
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(OUT) :: STR
        END SUBROUTINE f2kgetexe
      END INTERFACE
!
      CONTAINS
!
      PURE SUBROUTINE F2KSUBSTR(STRING,ISTART,IEND)
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
        IF (STRING(IPOS:IPOS) /= " ") THEN
            ISTART = IPOS
!
! Find the end of the sub-string
!
            IPOS2 = INDEX(STRING(IPOS:)," ")
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
      CHARACTER(LEN=2000), SAVE :: ARGSTR
      INTEGER            , SAVE :: LENARG,ISTART
      LOGICAL            , SAVE :: GETCMD = .TRUE.
      INTEGER                   :: INULL,IEND
!
!  Get whole command line including executable name and terminating
!  null. Remove executable and null before returning result.
!
      IF (GETCMD) THEN
          CALL f2kgetcl(ARGSTR)
          INULL = INDEX(ARGSTR,CHAR(0))
          IF (INULL > 0) THEN
              ARGSTR(INULL:) = " "
              CALL F2KSUBSTR(ARGSTR,ISTART,IEND)
              ISTART = IEND + 2
              LENARG = INULL - ISTART
          ELSE
              ARGSTR = " "
              ISTART = 1
              LENARG = 0
          END IF
          GETCMD = .FALSE.
      END IF
      IF (PRESENT(COMMAND)) COMMAND = ARGSTR(ISTART:)
      IF (PRESENT(LENGTH )) LENGTH  = LENARG
      IF (PRESENT(STATUS )) STATUS  = 0
      RETURN
      END SUBROUTINE GET_COMMAND
!
      FUNCTION COMMAND_ARGUMENT_COUNT() RESULT (RESCOMMAND_ARGUMENT_COUNT)
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
      INTEGER             :: RESCOMMAND_ARGUMENT_COUNT
      INTEGER             :: IPOS,ISTART,IEND,IPOS1,INULL
      CHARACTER(LEN=2000) :: ARGSTR
!
      INTEGER :: NARGS
!
!
! Get whole command line
! (null terminated and including executable name)
!
      CALL f2kgetcl(ARGSTR)
      INULL = INDEX(ARGSTR,CHAR(0))
      IF (INULL > 0) ARGSTR(INULL:) = " "
!
! Count command line arguments
!
      NARGS = -1
      IPOS  = 1
      DO
        CALL F2KSUBSTR(ARGSTR(IPOS:),ISTART,IEND)
        IF (ISTART > 0) THEN
            ISTART = ISTART + IPOS - 1
            IEND   = IEND   + IPOS - 1
            IPOS   = IEND   + 2
!
! Is argument quoted ?
!
            IF (ARGSTR(ISTART:ISTART) /= """") THEN
!
! No - increment arg count
!
                NARGS = NARGS + 1
            ELSE IF (ISTART < LEN(ARGSTR)) THEN
!
! Yes it is quoted and quote isn't at end of string
!
                ISTART = ISTART + 1
                CALL F2KSUBSTR(ARGSTR(ISTART:),IPOS1,IEND)
                IF (IPOS1 > 0) THEN
                    IEND = INDEX(ARGSTR(ISTART:),"""")
!
! Ignore null quotes
!
                    IF (IEND /= 1) THEN
                        IF (IEND == 0) THEN
                            IEND = LEN(ARGSTR)
                        ELSE
                            IEND = ISTART + IEND - 2
                        END IF
                        NARGS = NARGS + 1
                        IPOS  = IEND  + 3
                    END IF
                END IF
            END IF
!
            IF (IPOS <= LEN(ARGSTR)) CYCLE
        END IF
        EXIT
      END DO
!
      RESCOMMAND_ARGUMENT_COUNT = NARGS
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
      INTEGER             :: IPOS,ISTART,IEND,NARG,IPOS1,INULL
      CHARACTER(LEN=2000) :: ARGSTR
!
! If the argument number is negative, return an error code of 1.
!
      IF (NUMBER < 0) THEN
          IF (PRESENT(VALUE))  VALUE  = " "
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 1
          RETURN
      ELSE IF (NUMBER == 0) THEN
!
! If argument zero requested, get the executable name via the
! GetModuleFilename API function. The name at the start of the
! command line is not consistent across Windows 9x/Me and NT/2K.
! The former give the more useful full pathname whereas the latter
! only give the entered name. GetModuleFilename gives the full
! pathname regardless of Windows version.
!
          CALL f2kgetexe(ARGSTR)
          INULL = INDEX(ARGSTR,CHAR(0))
          IF (INULL > 0) ARGSTR(INULL:) = " "
          IF (PRESENT(VALUE))  VALUE  = ARGSTR
          IF (PRESENT(LENGTH)) LENGTH = MAX(INULL-1,0)
          IF (PRESENT(STATUS)) STATUS = 0
          RETURN
      END IF
!
! Get whole command line - remove terminating null
!
      CALL f2kgetcl(ARGSTR)
      INULL = INDEX(ARGSTR,CHAR(0))
      IF (INULL > 0) ARGSTR(INULL:) = " "
!
! Find required command line argument - skip executable name
!
      NARG = -1
      IPOS = 1
      DO
        CALL F2KSUBSTR(ARGSTR(IPOS:),ISTART,IEND)
        IF (ISTART > 0) THEN
            ISTART = ISTART + IPOS - 1
            IEND   = IEND   + IPOS - 1
            IPOS   = IEND   + 2
!
! Is argument quoted ?
!
            IF (ARGSTR(ISTART:ISTART) /= """") THEN
!
! No - increment arg count
!
                NARG = NARG + 1
            ELSE IF (ISTART < LEN(ARGSTR)) THEN
!
! Yes it is quoted and quote isn't at end of string
!
                ISTART = ISTART + 1
                CALL F2KSUBSTR(ARGSTR(ISTART:),IPOS1,IEND)
                IF (IPOS1 > 0) THEN
                    IEND = INDEX(ARGSTR(ISTART:),"""")
!
! Ignore null quotes
!
                    IF (IEND /= 1) THEN
                        IF (IEND == 0) THEN
                            IEND = LEN(ARGSTR)
                        ELSE
                            IEND = ISTART + IEND - 2
                        END IF
                        NARG = NARG + 1
                        IPOS = IEND + 3
                    END IF
                END IF
            END IF
!
! If this is the required command line argument, return value
! and exit otherwise continue if not at end of command line
!
            IF (NUMBER == NARG) THEN
                IF (PRESENT(VALUE))  VALUE  = ARGSTR(ISTART:IEND)
                IF (PRESENT(LENGTH)) LENGTH = IEND - ISTART + 1
                IF (PRESENT(STATUS)) STATUS = 0
                RETURN
            ELSE IF (IPOS <= LEN(ARGSTR)) THEN
                CYCLE
            END IF
        END IF
        EXIT
      END DO
!
! Error code = 2 : NUMBER too large
!
      IF (PRESENT(VALUE))  VALUE  = " "
      IF (PRESENT(LENGTH)) LENGTH = 0
      IF (PRESENT(STATUS)) STATUS = 2
      RETURN
      END SUBROUTINE GET_COMMAND_ARGUMENT
!
      END MODULE F2KCLI