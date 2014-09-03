/*
 F2KCLI : Fortran 200x Command Line Interface
 copyright Interactive Software Services Ltd. 2001-2003
 For conditions of use see manual.txt

 Platform    : Win32
 Compiler    : F and gcc (mingw32)
 To compile  : gcc -c f2kgetcl.c
 Implementer : Lawson B. Wakefield, I.S.S. Ltd.
 Date        : February 2003

 C binding to API GetCommandLine and GetModuleFilename functions

*/

#include <windows.h>

void f2kgetcl_(char *argstrptr,int argstrlen)
{
int ncopy;
ncopy = lstrlen(GetCommandLine()) + 1;
if (ncopy > argstrlen) ncopy = argstrlen;
lstrcpyn(argstrptr,GetCommandLine(),ncopy);
}

void f2kgetexe_(char *exestrptr,int exestrlen)
{
GetModuleFileName(NULL,exestrptr,exestrlen);
}
