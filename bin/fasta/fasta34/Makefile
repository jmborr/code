# $Name: fa34t20b3 $ - $Id: Makefile,v 1.52 2002/08/28 21:10:15 wrp Exp $
#
# makefile for fasta3, fasta3_t.  Use makefile.pvm for pvcompxx.
# this file works for DEC Alphas
#
# this file supports mmap()'ed databases in BLAST2 format use -DUSE_MMAP
# for mmap()ed BLAST2 format.

# the -DDEBUG option provides additional debugging information, particularly
# with -D on the command line.

# use -DBIG_LIB64 to generate and use 64-bit offsets in map_db .xin
# files

# for Tru64 4.0F, no "<inttypes.h>" 4.0G has inttypes.h

CC= cc -O -std1
#CC= insure -g -DDEBUG
#CC= cc -g -DDEBUG -std1

#CC= gcc -g -Wall
#
# standard line for normal searching
CFLAGS= -DM10_CONS -DUNIX -DTIMES -DHZ=60 -DBIGMEM -DMAX_WORKERS=4 -DTHR_EXIT=pthread_exit -DPROGRESS -DFASTA_HOST='"fasta.bioch.virginia.edu/fasta/cgi"' -DUSE_MMAP -DIS_LITTLE_ENDIAN -I/usr/local/include/mysql -DGAP_OPEN -DBIG_LIB64 -DMYSQL_DB -DSAMP_STATS
#
#(-DMYSQL_DB for mySQL databases)  (also requires change to Makefile34.common)

# special options for SUPERFAMLIES
#CFLAGS= -DM10_CONS -DUNIX -DTIMES -DHZ=60 -DBIGMEM -DSFCHAR="'|'" -DMAX_WORKERS=4 -DTHR_EXIT=pthread_exit -DPROGRESS -DSUPERFAMNUM -DIS_LITTLE_ENDIAN -DUSE_MMAP -DMAXBEST=200000 -DGAP_OPEN
 
#LIB_M = -lm
LIB_M = -L/usr/local/lib/mysql -lmysqlclient -lm
# for mySQL databases

HFLAGS= -o
NFLAGS= -o

#for DEC Unix V4.0 
THR_SUBS = pthr_subs2
THR_LIBS = -threads
THR_CC =

#for Sun
#THR_SUBS = uthr_subs
#THR_LIBS = -lthread
#THR_CC =
#
# for SGI with current pthreads
#THR_SUBS = pthr_subs
#THR_LIBS = -lpthreads
#THR_CC = 
#
# for IBM with current pthreads
#CC= xlc_r -v -g
#THR_SUBS = ibm_pthr_subs
#THR_LIBS = -lpthreads
#THR_CC =

XDIR = /seqprg/slib/bin

# renamed (fasta34)  programs
include Makefile34m.common_sql

