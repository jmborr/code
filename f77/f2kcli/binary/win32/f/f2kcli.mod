V9 0 0 3 0 0
MODULE F2KCLI,0 0
FILE 0,f2kcli.f90
PROC F2KSUBSTR,3,8,0: 8,0,0,0,0,60000,2,400800,0
VAR STRING,3,1,,: 3,0,a,-1,0,100010b,0,0,0
VAR ISTART,3,2,,: 1,0,3,0,0,83,0,0,0
VAR IEND,3,2,,: 1,0,3,0,0,83,0,0,0
ENDPROC
PROC F2KGETCL,1,10,0: 8,0,0,0,0,60000,2,800,0
VAR STR,3,2,,: 3,0,a,-1,0,1000003,0,0,0
ENDPROC
PROC GET_COMMAND_ARGUMENT,4,8,0: 8,0,0,0,0,40000,1,400000,0
VAR NUMBER,3,1,,: 1,0,3,0,0,103,0,0,0
VAR VALUE,3,2,,: 3,0,a,-1,0,100009b,0,0,0
VAR LENGTH,3,2,,: 1,0,3,0,0,9b,0,0,0
VAR STATUS,3,2,,: 1,0,3,0,0,9b,0,0,0
ENDPROC
PROC F2KGETEXE,1,10,0: 8,0,0,0,0,60000,2,800,0
VAR STR,3,2,,: 3,0,a,-1,0,1000003,0,0,0
ENDPROC
PROC GET_COMMAND,3,8,0: 8,0,0,0,0,40000,1,400000,0
VAR COMMAND,3,2,,: 3,0,a,-1,0,100009b,0,0,0
VAR LENGTH,3,2,,: 1,0,3,0,0,9b,0,0,0
VAR STATUS,3,2,,: 1,0,3,0,0,9b,0,0,0
ENDPROC
PROC COMMAND_ARGUMENT_COUNT,0,8,0: 1,0,3,0,0,40081,1,10400000,0
RESULTVAR RESCOMMAND_ARGUMENT_COUNT,4,0,,COMMAND_ARGUMENT_COUNT: 1,0,3,0,0,81,0,0,0
ENDPROC
END
