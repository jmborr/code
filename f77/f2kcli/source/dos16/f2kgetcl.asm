; F2KCLI : Fortran 200x Command Line Interface
; copyright Interactive Software Services Ltd. 2001
; For conditions of use see manual.txt
;
; Platform    : DOS (16-bit real mode)
; Compiler    : Lahey F77L
;               Microsoft Fortran 4.x/5.x
;               Prospero Fortran
;               Watcom Fortran 77/16
; To compile  : masm /DLAHEY f2kgetcl;   (Lahey)
;               masm f2kgetcl;           (MS)
;               masm f2kgetcl;           (Prospero)
;               masm /DWATCOM f2kgetcl;  (Watcom)
; Implementer : Lawson B. Wakefield, I.S.S. Ltd.
; Date        : February 2001
;
IFDEF WATCOM
arg1_1   EQU    2
ELSE
arg1_1   EQU    6
ENDIF
;
data     SEGMENT WORD PUBLIC 'DATA'
data     ENDS
DGROUP   GROUP data
codeseg  SEGMENT BYTE PUBLIC 'CODE'
         ASSUME CS:codeseg,DS:DGROUP,SS:DGROUP
;
;  SUBROUTINE F2KGETCL(STRING)
;
;  Get program command line - method depends on version of DOS
;
;  DOS 3.0 or later : Use documented DOS function 62h to get PSP address
;
;  DOS 2.x          : Use undocumented DOS function 51h
;
;  STRING = Receiving string for command line
;
;  INTEGER STRING
; (should be equivalenced in calling routine to a CHARACTER variable
;  - non standard but all target compilers allow this)
;
         PUBLIC   f2kgetcl
f2kgetcl PROC FAR
;
IFDEF WATCOM
         PUSH CX
         PUSH BX
         PUSH DX
         PUSH AX
ENDIF
         PUSH BP
         MOV  BP,SP
;
         PUSH DS
         MOV  AH,30h
         INT  21h          ; Get DOS version
         MOV  AH,62h       ; By default assume documented DOS 3 function
         CMP  AL,02h       ; Is it DOS 2 ?
         JA   DOS_1        ; If it's 3 or more call DOS
         MOV  AH,51h       ; DOS 2 : use undocumented call
DOS_1:   INT  21h          ; call DOS to get original address of PSP
         LES  DI,[BP+arg1_1] ; make ES:DI point to receiving string
         MOV  DS,BX        ; PSP is returned in BX
         MOV  SI,0080h     ; point SI to byte containing length of cmd line
         CLD
         LODSB             ; get length of command line into AX & increment SI
         AND  AX,007Fh     ; return up to 127 chars
         MOV  CX,AX        ; transfer number of bytes to copy into CX
         CMP  CL,0
         JE   DONE
         REPZ MOVSB        ; move the string from DS:SI to ES:DI
                           ; (note : passed string is assumed to be
                           ;         already blank filled by calling prog)
DONE:    POP  DS
;
         POP  BP
IFDEF WATCOM
         POP  AX
         POP  DX
         POP  BX
         POP  CX
         RET
ELSE
IFDEF LAHEY
         RET
ELSE
         RET 4
ENDIF
ENDIF

f2kgetcl ENDP
;
;  SUBROUTINE F2KENVAD(IADDR)
;
;  Get address of environment
;
;  IADDR = Environment address
;
         PUBLIC   f2kenvad
f2kenvad PROC FAR
;
IFDEF WATCOM
         PUSH CX
         PUSH BX
         PUSH DX
         PUSH AX
ENDIF
         PUSH BP
         MOV  BP,SP
;
         PUSH  DS
         MOV   AH,30h
         INT   21h          ; Get DOS version
         MOV   AH,62h       ; By default assume documented DOS 3 function
         CMP   AL,02h       ; Is it DOS 2 ?
         JA    DOS_2          ; If it's 3 or more call DOS
         MOV   AH,51h       ; DOS 2 : use undocumented call
DOS_2:   INT   21h          ; call DOS to get original address of PSP
         MOV   DS,BX        ; PSP address is returned in BX
         MOV   SI,002Ch     ; DS:SI points to env address in PSP
         LODSW              ; get environment address
         POP   DS
         XOR   CX,CX        ; environment starts at offset zero
         LES   BX,[BP+arg1_1]
         MOV   ES:[BX],CX   ; return address offset and segment
         MOV   ES:[BX+2],AX
;
         POP  BP
IFDEF WATCOM
         POP  AX
         POP  DX
         POP  BX
         POP  CX
         RET
ELSE
IFDEF LAHEY
         RET
ELSE
         RET 4
ENDIF
ENDIF

f2kenvad ENDP
;
;  SUBROUTINE F2KPEEK(IADDR,IBYTE)
;
;  Peek a byte from a specified address & return it in an INT*2 variable
;
;  (I ) IADDR = Byte address
;  ( O) IBYTE = Returned byte value
;
;  INTEGER*4 IADDR
;  INTEGER*2 IBYTE
;
IFDEF WATCOM
arg2_1   EQU    2
arg2_2   EQU    6
ELSE
IFDEF LAHEY
arg2_1   EQU    6
arg2_2   EQU    10
ELSE
arg2_1   EQU    10
arg2_2   EQU    6
ENDIF
ENDIF

         PUBLIC f2kpeek
f2kpeek  PROC  FAR
;
IFDEF WATCOM
         PUSH  CX
         PUSH  BX
         PUSH  DX
         PUSH  AX
ENDIF
         PUSH  BP
         MOV   BP,SP
         PUSH  SI
         PUSH  DI
;
         PUSH  DS
         LES   BX,[BP+arg2_1]
         MOV   AX,ES:[BX]
         MOV   DS,ES:[BX+2]
         MOV   SI,AX        ; DS:SI points to byte to get
         LODSB              ; load byte into AL
         POP   DS
         XOR   AH,AH        ; clear top half of AX
         LES   BX,[BP+arg2_2]
         MOV   ES:[BX],AX   ; return 2 byte word
;
         POP   DI
         POP   SI
         POP   BP
IFDEF WATCOM
         POP   AX
         POP   DX
         POP   BX
         POP   CX
         RET
ELSE
IFDEF LAHEY 
         RET
ELSE
         RET 8
ENDIF
ENDIF

f2kpeek  ENDP
;
codeseg  ENDS
         END
