F2KCLI for WATFOR77
-------------------

The following contribution from E. P. Chandler describes how to implement
F2KCLI for use with WATFOR77. Since we are not in a position to test the
described implementation, it is reproduced here on an "as is" basis, in
case it is useful to other users of that compiler.

Interactive Software Services Ltd. Oct 2005

Original submission from E.P. Chandler:
===================================================

Greetings:

Congratulations to Lawson B. Wakefield for a very nice package.

I noticed that this package is listed to work with WATCOM F77 but not the
older WATFOR77. I've put together a list of changes needed to the source
file F2KGETCLI.ASM in the 16 bit dos sources and to F2KCLI.INC, F2KCLI.FOR
and TESTCLI.FOR in order for these programs to run under WATFOR77
(when compiled to .exe files). So far I have had success with versions
1.4, 2.0 and 3.7 ('77 & '87).

Please feel free to include the following in your distribution, if you
choose. (I have done _some_ limited testing but the user is at his own risk.)

There are no restrictions on its use as long as my name and e-mail address
remain in the source code.

---- start text ----
edit testcli.for:

change INCLUDE 'f2kcli.inc' to *$include f2kcli.inc
after main program add *$include f2kcli.for

edit f2kgetcl.asm:

insert before first IFDEF WATCOM in f2kgetcl & f2kenvad

IFDEF WATFOR77
; create a MSF stack frame
         POP AX
         POP DX            ;save return address
         LDS SI,ES:[BX]    ;get ptr to 1st arg
         PUSH DS
         PUSH SI           ;push it
         PUSH DX
         PUSH AX           ;push return address
ENDIF

insert before first IFDEF WATCOM in f2kpeek

IFDEF WATFOR77
; create a MSF stack frame
         POP AX
         POP DX            ;save return address
         LDS SI,ES:[BX]    ;get ptr to 1st arg
         PUSH DS
         PUSH SI           ;push it
         LDS SI,ES:[BX+4]  ;get ptr to 2nd arg
         PUSH DS
         PUSH SI           ;push it
         PUSH DX
         PUSH AX           ;push return address
ENDIF

to compile:

masm /DWATFOR77 f2kgetcl;

to set up:

lib f2kcli +f2kgetcl;
set library=f2kcli;watfor
copy f2kcli.inc f2kcli.for f2kcli.lib & watfor.lib
  to working directory
watfor77/exe/noedit testcli

Additional changes are needed in version 1.4:

1. change _ to $ in variable names in f2kcli.for, f2kcli.inc & testcli.for
2. remove ARGSTR after SAVE in f2kcli.for

Notes:

1. Watfor77 passes a far pointer to an array of far pointers
   to arguments (or string descriptors) in ES:BX.
2. The string argument to the first subroutine is actually passed
   as an integer (via equivalence) so the argument points directly
   to the string, not to a string descriptor.
3. I have added a prolog to each subroutine which changes from
   Watfor77's calling convention to MS-FORTRAN. The code then falls
   into the MS-FORTRAN routine.
4. F2KCLI.FOR is included as source code since WATFOR77 does not
   create OBJ files.


E. P. Chandler
e-mail: epc8 at juno.com
---- end text -----

===================================================
