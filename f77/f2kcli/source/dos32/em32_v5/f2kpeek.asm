; F2KCLI : Fortran 200x Command Line Interface
; copyright Interactive Software Services Ltd. 2001
; For conditions of use see manual.txt
;
; Platform    : DOS (Phar Lap 32-bit DOS extender)
; Compiler    : Lahey F77L-EM/32 v5.x
; To compile  : masm f2kgetcl;
; Implementer : Lawson B. Wakefield, I.S.S. Ltd.
; Date        : February 2001
;
; SUBROUTINE F2KPEEK(IADDR,IBYTE)
;
; Peek a byte from a specified address & return it in an INT*2 variable
;
; IADDR = Byte address
; IBYTE = Returned byte value
;
; INTEGER*4 IADDR
; INTEGER*2 IBYTE
;
arg1     EQU    8
arg2     EQU    12

.386
data     SEGMENT WORD PUBLIC 'DATA'
data     ENDS

DGROUP   GROUP data

codeseg  SEGMENT BYTE PUBLIC 'CODE'
         ASSUME  CS:codeseg,DS:DGROUP,SS:DGROUP
         PUBLIC  f2kpeek
;
f2kpeek  PROC  NEAR
         PUSH  EBP
         MOV   EBP,ESP
         PUSH  EBX
         PUSH  ESI
         PUSH  EDI
;
         PUSH  DS
         MOV   EBX,DWORD PTR[EBP+arg1]
         MOV   AX,WORD PTR[EBX]
         MOV   DS,WORD PTR[EBX+2]
         XOR   ESI,ESI
         MOV   SI,AX                    ; DS:SI points to byte to get
         LODSB                          ; load byte into AL
         POP   DS
;
         XOR   AH,AH                    ; clear top half of AX
         MOV   EBX,DWORD PTR[EBP+arg2]  ; return 2 byte word
         MOV   WORD PTR[EBX],AX
;
         POP   EDI
         POP   ESI
         POP   EBX
         POP   EBP
         RET
f2kpeek  ENDP
;
codeseg  ENDS
         END
