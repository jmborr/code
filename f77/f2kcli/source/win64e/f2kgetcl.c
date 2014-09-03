/*
 F2KCLI : Fortran 200x Command Line Interface
 copyright Interactive Software Services Ltd. 2001-2005
 For conditions of use see manual.txt

 Platform    : Win64 (AMD64/EM64T)
 Compiler    : Microsoft C
 To compile  : cl /c /DUPPER /Zl f2kgetcl.c
 Implementer : Lawson B. Wakefield, I.S.S. Ltd.
 Date        : October 2005

 C binding to API GetCommandLine and GetModuleFilename functions

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

void f2kgetcl(char *argstrptr,unsigned long argstrlen)
{
int ncopy;
ncopy = lstrlen(GetCommandLine()) + 1;
if (ncopy > argstrlen) ncopy = argstrlen;
lstrcpyn(argstrptr,GetCommandLine(),ncopy);
}

void f2kgetexe(char *exestrptr,unsigned long exestrlen)
{
GetModuleFileName(NULL,exestrptr,exestrlen);
}
