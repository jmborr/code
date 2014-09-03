/*
 F2KCLI : Fortran 200x Command Line Interface
 copyright Interactive Software Services Ltd. 2001-2003
 For conditions of use see manual.txt

 Platform    : Win32
 Compiler    : Most Fortran-compatible Win32 C compilers
 To compile  : see below
 Implementer : Lawson B. Wakefield, I.S.S. Ltd.
 Date        : February 2001
               Updated May 2002

 C binding to API GetCommandLine and GetModuleFilename functions

 Fortran             C         C compiler command line
 -------             -         -----------------------
 Absoft f77/f9x      Absoft    acc -c -windefs f2kgetcl.c
 Absoft f77/f9x      Absoft    acc -c -windefs -DUPPER f2kgetcl.c
 GNU g77             GNU       gcc -c -DUS f2kgetcl.c
 Intel               MS VC++   cl /c /DUPPER f2kgetcl.c
 Lahey Elf90         Borland   bcc32 -c -W -X -DUPPER -u- f2kgetcl.c
 Lahey LF90          Borland   bcc32 -c -W -X -DUPPER -u- f2kgetcl.c
 Lahey LF95 5.0-5.6  Borland   bcc32 -c -W -X -DUS f2kgetcl.c
 Lahey LF95 5.7/7.x  MS VC++   cl /c /DUS /Gd f2kgetcl.c
 MS PowerStation     MS VC++   cl /c /Gz /DUPPER f2kgetcl.c
 Salford FTNxx       Salford   scc /ansi_c /def UPPER f2kgetcl.c
 VF(x86)             MS VC++   cl /c /DUPPER /Gz f2kgetcl.c
 VF(AXP)             MS VC++   cl /c /DUPPER f2kgetcl.c
 Watcom              Watcom    wcc386 /mf /ei /zp4 /bt=nt_win /i=%watcom%\h
                                /i=%watcom%\h\nt /dUPPER /dWATCOM f2kgetcl.c
*/

#include <windows.h>

#ifdef US
#define f2kgetcl  f2kgetcl_
#define f2kgetexe f2kgetexe_
#else
#ifdef UPPER
#define f2kgetcl  F2KGETCL
#define f2kgetexe F2KGETEXE
#endif
#endif

#ifdef WATCOM

struct character {char *ptr; int len;};

void f2kgetcl(struct character *argstr)
{
int ncopy;
ncopy = lstrlen(GetCommandLine()) + 1;
if (ncopy > argstr->len) ncopy = argstr->len;
lstrcpyn(argstr->ptr,GetCommandLine(),ncopy);
}

void f2kgetexe(struct character *exestr)
{
GetModuleFileName(NULL,exestr->ptr,exestr->len);
}

#else

void f2kgetcl(char *argstrptr,int argstrlen)
{
int ncopy;
ncopy = lstrlen(GetCommandLine()) + 1;
if (ncopy > argstrlen) ncopy = argstrlen;
lstrcpyn(argstrptr,GetCommandLine(),ncopy);
}

void f2kgetexe(char *exestrptr,int exestrlen)
{
GetModuleFileName(NULL,exestrptr,exestrlen);
}

#endif
