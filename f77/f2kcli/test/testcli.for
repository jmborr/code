C F2KCLI : Fortran 200x Command Line Interface
C copyright Interactive Software Services Ltd. 2001
C For conditions of use see manual.txt
C
C Fixed format Fortran 77 test program
C
      PROGRAM TESTCLI
C
      CHARACTER*256 LINE,EXE
      CHARACTER*40  CMD
      INTEGER       NARG,IARG
C
      INCLUDE 'f2kcli.inc'
C
      NARG = COMMAND_ARGUMENT_COUNT()
      WRITE(*,*) 'Arg count=', NARG
C
      CALL GET_COMMAND(LINE,LENGTH,ISTATUS)
      WRITE(*,*) 'Line=',LINE
C
      CALL GET_COMMAND_ARGUMENT(0,EXE,LENGTH,ISTATUS)
      WRITE(*,*) 'Program=',EXE
C
      DO 100 IARG = 1,NARG
        CALL GET_COMMAND_ARGUMENT(IARG,CMD,LENGTH,ISTATUS)
        WRITE(*,*) 'Arg ',IARG,'=',CMD
  100 CONTINUE
      END

