$! make_vms.com  
$! OpenVMS build procedure for FASTA 3
$!
$! from:
$! David Mathog
$! mathog@seqaxp.bio.caltech.edu
$! Manager, sequence analysis facility, biology division, Caltech 
$!
$!mycc was gcc -O -DUNIX -DTIMES -DHZ=100 -DBIGMEM \
$! -DSFCHAR="':'" -c -DMAX_WORKERS=4 -DTHR_EXIT=pthread_exit \
$! -DPROGRESS -DLINUX
$!
$! Hey Bill, how about writing to a SINGLE C standard :-), this sucker 
$! requires both POSIX and XOpen extensions
$!    added _POSIX_C_SOURCE=2 to pick up GETOPT, isatty(), SIGHUP, and others
$!    defined _XOPEN_SOURCE_EXTENDED to pick up MAXFLOAT and others
$! no TIMES for OpenVMS <7, at 7 or up it should be ok as tms will be present
$!
$! disable ARGLISGTR255 warning - not a problem for Alpha, might be
$! on VAX because of the way it keeps track of stack arguments, and
$! some of the routines here pass big structures, which put them over the
$! 255 byte limit.
$!
$ docdebug = "/debug/noopt"
$ doldebug = "/debug"
$ if (P1 .eqs. "" )
$ then
$   docdebug = ""
$   doldebug = ""
$ else
$   docdebug = "/debug/noopt"
$   doldebug = "/debug"
$ endif
$!
$!  /reentrancy=multithread required for _t variants, doesn't seem
$!  to hurt the others.  Without it, _t variants only run up to 100%
$!  of one CPU.
$!
$ mycc :== cc/standard=ansi89/prefix=all/include=[] 'docdebug' -
  /warn=(disable=ARGLISGTR255)/reentrancy=multithread -
  /define=(M10_CONS,M10_CONS_L,UNIX, -
  HZ=60,BIGMEM,SFCHAR="""':'""",STAR_X,MAX_WORKERS=4,-
  THR_EXIT="""pthread_exit""",PROGRESS,_POSIX_C_SOURCE=2,_XOPEN_SOURCE_EXTENDED)
$ mycc02 == mycc - "EXTENDED)" + "EXTENDED,TFASTA)"
$ mycc03 == mycc - "EXTENDED)" + "EXTENDED,TFASTX)"
$ mycc04 == mycc - "EXTENDED)" + "EXTENDED,FASTX)"
$ mycc05 == mycc - "EXTENDED)" + "EXTENDED,FASTA)"
$ mycc06 == mycc - "EXTENDED)" + "EXTENDED,PRSS)"
$ mycc07 == mycc - "EXTENDED)" + "EXTENDED,NOLIB)"
$ mycc08 == mycc - "EXTENDED)" + "EXTENDED,FASTF)"
$ mycc09 == mycc - "EXTENDED)" + "EXTENDED,TFASTF)"
$ mycc10 == mycc - "EXTENDED)" + "EXTENDED,FASTS)"
$ mycc11 == mycc - "EXTENDED)" + "EXTENDED,TFASTS)"
$!
$! stolen from OSU link_server.com
$!
$ if f$search("sys$share:pthread$rtl.exe") .eqs. ""
$ then
$!   Classic DECthreads (draft 4).
$   create link.opt
sys$share:cma$lib_shr/share,cma$rtl/share
sys$share:cma$open_lib_shr/share,cma$open_rtl/share
$ kt_enable = ""
$ else
$!   POSIX-based DECthreads (VMS V7).  Link with flag to enable upcalls, but
$!   only if multithreading upcalls turned on.
$   create link.opt
!   nothing here
$   kt_enable = "/THREADS_ENABLE=(UPCALLS,MULTIPLE)"
$ endif
$
$ mylink :== link 'doldebug' /nomap'kt_enable'
$!
$!
$ mycc  url_subs.c 
$ mycc  complib.c 
$ rename complib.obj comp_lib.obj
$ mycc  compacc.c 
$ mycc  showbest.c 
$ mycc  showalign.c 
$ rename showalign.obj showalig.obj
$ mycc  htime.c
$ mycc  apam.c
$ mycc  doinit.c 
$ mycc05 initfa.c 
$ rename initfa.obj init_fa.obj
$ mycc  scaleswn.c 
$ mycc  dropnfa.c 
$ rename dropnfa.obj drop_nfa.obj
$ mycc  nxgetaa.c 
$ rename nxgetaa.obj lgetaa.obj
$ mycc  c_dispn.c 
$ mycc  ncbl_lib.c 
$ mycc  lib_sel.c 
$ mycc  nrand48.c 
$ mylink/exe=fasta3.exe -
comp_lib.obj,compacc.obj,showbest.obj,showalig.obj,-
htime.obj,apam.obj,doinit.obj,init_fa.obj,drop_nfa.obj,scaleswn.obj,-
lgetaa.obj,c_dispn.obj,ncbl_lib.obj,lib_sel.obj,url_subs.obj,nrand48.obj
$ mycc  initsw.c 
$ rename initsw.obj init_sw.obj
$ mycc  dropgsw.c 
$ mylink/exe=ssearch3.exe -
comp_lib.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,doinit.obj,init_sw.obj,-
dropgsw.obj,scaleswn.obj,lgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,url_subs.obj,nrand48.obj
$ mycc06 complib.c 
$ rename complib.obj rcomplib.obj
$ mycc06 initsw.c 
$ rename initsw.obj init_rss.obj
$ mycc  scaleswe.c 
$ mycc  dropnsw.c 
$ mycc07 llgetaa.c 
$ mycc  showrss.c 
$ mylink/exe=prss3.exe -
rcomplib.obj,compacc.obj,htime.obj,apam.obj,-
doinit.obj,init_rss.obj,dropnsw.obj,scaleswe.obj,llgetaa.obj,-
showrss.obj,lib_sel.obj,nrand48.obj
$ mycc08 initfa.c 
$ rename initfa.obj init_ff.obj
$ mycc  scaleswg.c 
$ mycc  dropffa.c 
$ rename dropffa.obj drop_ff.obj
$ mylink/exe=fastf3.exe -
comp_lib.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,doinit.obj,init_ff.obj,-
drop_ff.obj,scaleswg.obj,lgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,nrand48.obj,url_subs.obj
$ mycc10 initfa.c 
$ rename initfa.obj init_fs.obj
$ mycc10 dropfs.c 
$ rename dropfs.obj drop_fs.obj
$ mylink/exe=fasts3.exe -
comp_lib.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,doinit.obj,init_fs.obj,-
drop_fs.obj,scaleswn.obj,lgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,nrand48.obj,url_subs.obj
$ mycc02 complib.c 
$ rename complib.obj tcomplib.obj
$ mycc09 initfa.c 
$ rename initfa.obj init_tf.obj
$ mycc09 dropffa.c 
$ rename dropffa.obj drop_tff.obj
$ mycc02 nxgetaa.c 
$ rename nxgetaa.obj tgetaa.obj
$ mycc  faatran.c 
$ mylink/exe=tfastf3.exe -
tcomplib.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,doinit.obj,init_tf.obj,-
drop_tff.obj,scaleswg.obj,tgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,nrand48.obj,url_subs.obj
$ mycc02 initfa.c 
$ rename initfa.obj init_tfa.obj
$ mycc02 dropnfa.c 
$ rename dropnfa.obj drop_tfa.obj
$ mylink/exe=tfasta3.exe -
tcomplib.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,doinit.obj,init_tfa.obj,-
drop_tfa.obj,scaleswn.obj,tgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj
$ mycc04 complib.c 
$ rename complib.obj complibx.obj
$ mycc04 showalign.c 
$ rename showalign.obj showaligx.obj
$ mycc04 initfa.c 
$ rename initfa.obj initfx.obj
$ mycc04 dropfx.c 
$ rename dropfx.obj drop_fx.obj
$ mylink/exe=fastx3.exe -
complibx.obj,compacc.obj,showbest.obj,-
showaligx.obj,htime.obj,apam.obj,doinit.obj,initfx.obj,-
drop_fx.obj,scaleswn.obj,lgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj
$ mycc03 complib.c 
$ rename complib.obj tcomplibx.obj
$ mycc03 showalign.c 
$ rename showalign.obj showalitx.obj
$ mycc03 initfa.c 
$ rename initfa.obj inittfx.obj
$ mycc03 dropfx.c 
$ rename dropfx.obj tdropfx.obj
$ mylink/exe=tfastx3.exe -
tcomplibx.obj,compacc.obj,showbest.obj,-
showalitx.obj,htime.obj,apam.obj,doinit.obj,inittfx.obj,-
tdropfx.obj,scaleswn.obj,tgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj
$ mycc04 dropfz.c 
$ rename dropfz.obj drop_fz.obj
$ mylink/exe=fasty3.exe -
complibx.obj,compacc.obj,showbest.obj,-
showaligx.obj,htime.obj,apam.obj,doinit.obj,initfx.obj,-
drop_fz.obj,scaleswn.obj,lgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj
$ mycc03 dropfz.c 
$ rename dropfz.obj tdropfz.obj
$ mylink/exe=tfasty3.exe -
tcomplibx.obj,compacc.obj,showbest.obj,-
showalitx.obj,htime.obj,apam.obj,doinit.obj,inittfx.obj,-
tdropfz.obj,scaleswn.obj,tgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj
$ mycc11 initfa.c 
$ rename initfa.obj init_ts.obj
$ mycc11 dropfs.c 
$ rename dropfs.obj drop_tfs.obj
$ mylink/exe=tfasts3.exe -
tcomplib.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,doinit.obj,init_ts.obj,-
drop_tfs.obj,scaleswn.obj,tgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,nrand48.obj,url_subs.obj
$ mycc  comp_thr.c 
$ mycc  work_thr.c 
$ mycc  pthr_subs.c 
$ mylink/exe=ssearch3_t.exe -
comp_thr.obj,work_thr.obj,pthr_subs.obj,-
compacc.obj,showbest.obj,showalig.obj,htime.obj,apam.obj,doinit.obj,-
init_sw.obj,dropgsw.obj,scaleswn.obj,-
lgetaa.obj,c_dispn.obj,ncbl_lib.obj,lib_sel.obj,url_subs.obj,nrand48.obj 
$ mylink/exe=fasta3_t.exe -
comp_thr.obj,work_thr.obj,pthr_subs.obj,compacc.obj,-
showbest.obj,showalig.obj,htime.obj,apam.obj,doinit.obj,-
init_fa.obj,drop_nfa.obj,scaleswn.obj,lgetaa.obj,c_dispn.obj,-
ncbl_lib.obj,lib_sel.obj,url_subs.obj,nrand48.obj 
$ mylink/exe=fastf3_t.exe -
comp_thr.obj,work_thr.obj,pthr_subs.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,-
doinit.obj,init_ff.obj,drop_ff.obj,scaleswg.obj,lgetaa.obj,-
c_dispn.obj,ncbl_lib.obj,lib_sel.obj,-
nrand48.obj,url_subs.obj 
$ mylink/exe=fasts3_t.exe -
comp_thr.obj,work_thr.obj,pthr_subs.obj,compacc.obj,showbest.obj,showalig.obj,-
htime.obj,apam.obj,doinit.obj,init_fs.obj,drop_fs.obj,scaleswn.obj,lgetaa.obj,-
c_dispn.obj,ncbl_lib.obj,lib_sel.obj,nrand48.obj,url_subs.obj 
$ mycc02 comp_thr.c 
$ rename comp_thr.obj tcomp_thr.obj
$ mylink/exe=tfastf3_t.exe -
 tcomp_thr.obj,work_thr.obj,-
pthr_subs.obj,compacc.obj,showbest.obj,-
showalig.obj,htime.obj,apam.obj,doinit.obj,-
init_tf.obj,drop_tff.obj,scaleswg.obj,tgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,nrand48.obj,url_subs.obj 
$ mylink/exe=tfasta3_t.exe -
 tcomp_thr.obj,work_thr.obj,-
pthr_subs.obj,compacc.obj,showbest.obj,showalig.obj,-
htime.obj,apam.obj,doinit.obj,-
init_tfa.obj,drop_tfa.obj,scaleswn.obj,tgetaa.obj,-
c_dispn.obj,ncbl_lib.obj,lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj 
$ mycc04 comp_thr.c 
$ rename comp_thr.obj comp_thrx.obj
$ mylink/exe=fastx3_t.exe-
 comp_thrx.obj,work_thr.obj,pthr_subs.obj,-
compacc.obj,showbest.obj,showaligx.obj,htime.obj,apam.obj,doinit.obj,-
initfx.obj,drop_fx.obj,faatran.obj,scaleswn.obj,lgetaa.obj,-
c_dispn.obj,ncbl_lib.obj,lib_sel.obj,url_subs.obj,nrand48.obj 
$ mycc03 comp_thr.c 
$ rename comp_thr.obj tcomp_thrx.obj
$ mylink/exe=tfastx3_t.exe -
tcomp_thrx.obj,work_thr.obj,-
pthr_subs.obj,compacc.obj,showbest.obj,showalitx.obj,htime.obj,apam.obj,-
doinit.obj,inittfx.obj,tdropfx.obj,scaleswn.obj,tgetaa.obj,c_dispn.obj,-
ncbl_lib.obj,lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj 
$ mylink/exe=fasty3_t.exe -
 comp_thrx.obj,work_thr.obj,-
pthr_subs.obj,compacc.obj,showbest.obj,showaligx.obj,-
htime.obj,apam.obj,doinit.obj,-
initfx.obj,drop_fz.obj,faatran.obj,scaleswn.obj,lgetaa.obj,-
c_dispn.obj,ncbl_lib.obj,lib_sel.obj,url_subs.obj,nrand48.obj 
$ mylink/exe=tfasty3_t.exe -
 tcomp_thrx.obj,work_thr.obj,pthr_subs.obj,-
compacc.obj,showbest.obj,showalitx.obj,htime.obj,-
apam.obj,doinit.obj,inittfx.obj,tdropfz.obj,scaleswn.obj,tgetaa.obj,-
c_dispn.obj,ncbl_lib.obj,lib_sel.obj,faatran.obj,url_subs.obj,nrand48.obj 
$ mylink/exe=tfasts3_t.exe -
 tcomp_thr.obj,work_thr.obj,pthr_subs.obj,compacc.obj,-
showbest.obj,showalig.obj,htime.obj,apam.obj,doinit.obj,init_ts.obj,-
drop_tfs.obj,scaleswn.obj,tgetaa.obj,c_dispn.obj,ncbl_lib.obj,-
lib_sel.obj,faatran.obj,nrand48.obj,url_subs.obj 
