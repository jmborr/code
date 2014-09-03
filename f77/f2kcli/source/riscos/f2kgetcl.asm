    TTL "f2kgetcl"

; F2KCLI : Fortran 200x Command Line Interface
; copyright Interactive Software Services Ltd. 2001
; For conditions of use see manual.txt
;
; Platform    : RISC OS
; Compiler    : Acornsoft Fortran 77
; To compile  : objasm asm.f2kgetcl aof.f2kgetcl -quit -stamp
; Implementer : Lawson B. Wakefield, I.S.S. Ltd.
; Date        : February 2001
;

; Registers

R0  RN   0
R1  RN   1
R2  RN   2
R3  RN   3
R4  RN   4
R5  RN   5
R6  RN   6
R7  RN   7
R8  RN   8
R9  RN   9

FP  RN  10
SP  RN  12
SB  RN  13
R14 RN  14
PC  RN  15

F0  FN   0
F1  FN   1
F2  FN   2
F3  FN   3
F4  FN   4
F5  FN   5
F6  FN   6
F7  FN   7

OS_GetEnv              * &10

; Data

        AREA   |F77$$Data|,DATA
RESULT  %      4               ; Result location

; Code

        AREA   |F77$$CODE|,CODE,READONLY
DATAPTR DCD    |F77$$Data|

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

NAME01  DCB   "F2KGETCL    "
        EXPORT F2KGETCL

;  SUBROUTINE F2KGETCL(ARGSTR)
;
;  Fortran callable binding to OS_GetEnv
;
;  ARGSTR = Command line string
;

F2KGETCL       ADR   R1,NAME01+12
               STMFD SP!,{R1,FP,SB,R14}
               LDR   SB,DATAPTR      ; Address data area
               MOV   FP,R0           ; copy argument list
;
               SWI   OS_GetEnv       ; get address of command line string
                                     ; into R0
               LDR   R1,[FP]         ; get return string descriptor address
               LDMIA R1,{R1,R2}      ; return string address/length in R1/R2
;
COPYARG        LDRB  R3,[R0],#1      ; get character & increment pointer
               STRB  R3,[R1],#1      ; store character & increment pointer
               SUBS  R2,R2,#1        ; decrement loop counter
               BGT   COPYARG         ; loop till end of return string
;
               MOV   R0,SB
               LDMIA SP!,{R1,FP,SB,PC}

               END
