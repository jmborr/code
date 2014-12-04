/*	ncbl2_lib.c	functions to read ncbi-blast format files from
			formatdb (blast2.0 format files)

		copyright (c) 1999 William R. Pearson
*/

/* $Name: fa34t20b3 $ - $Id: ncbl2_mlib.c,v 1.24 2002/08/02 00:49:33 wrp Exp $ */

/* to turn on mmap()ing for Blast2 files: */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* ****************************************************************
This code reads NCBI Blast2 format databases from formatdb version 3

(From NCBI) This section describes the format of the databases.

Formatdb creates three main files for proteins containing indices,
sequences, and headers with the extensions, respectively, of pin, psq,
and phr (for nucleotides these are nin, nsq, and nhr).  A number of
other ISAM indices are created, but these are described elsewhere.

FORMAT OF THE INDEX FILE
------------------------

1.) formatdb version number 	[4 bytes].

2.) protein dump flag (1 for a protein database, 0 for a nucleotide
    database) [4 bytes].

3.) length of the database title in bytes	[4 bytes].
4.) the database title		[length given in 3.)].
5.) length of the date/time string	[4 bytes].
6.) the date/time string	[length given in 5.)].
7.) the number of sequences in the database	[4 bytes].
8.) the total length of the database in residues/basepairs	[4 bytes].
9.) the length of the longest sequence in the database 		[4 bytes].

10.) a list of the offsets for definitions (one for each sequence) in
the header file.  There are num_of_seq+1 of these, where num_of_seq is
the number of sequences given in 7.).

11.) a list of the offsets for sequences (one for each sequence) in
the sequence file.  There are num_of_seq+1 of these, where num_of_seq
is the number of sequences given in 7.).

12.) a list of the offsets for the ambiguity characters (one for each
sequence) in the sequence file.  This list is only present for
nucleotide databases and, since the database is compressed 4/1 for
nucleotides, allows the ambiguity characters to be restored when the
sequence is generated.  There are num_of_seq+1 of these, where
num_of_seq is the number of sequences given in 7.).


FORMAT OF THE SEQUENCE FILE
---------------------------

There are different formats for the protein and nucleotide sequence files.

The protein sequence files is quite simple.  The first byte in the
file is a NULL byte, followed by the sequence in ncbistdaa format
(described in the NCBI Software Development Toolkit documentation).
Following the sequence is another NULL byte, followed by the next
sequence.  The file ends with a NULL byte, following the last
sequence.

The nucleotide sequence file contains the nucleotide sequence, with
four basepairs compressed into one byte.  The format used is NCBI2na,
documented in the NCBI Software Development Toolkit manual.  Any
ambiguity characters present in the original sequence are replaced at
random by A, C, G or T.  The true value of ambiguity characters are
stored at the end of each sequence to allow true reproduction of the
original sequence.

FORMAT OF THE HEADER FILE
-------------------------

The format of the header file depends on whether or not the identifiers in the
original file were parsed or not.  For the case that they were not, then each
entry has the format:

gnl|BL_ORD_ID|entry_number my favorite yeast sequence...

Here entry_number gives the ordinal number of the sequence in the
database (with zero offset).  The identifier
gnl|BL_ORD_ID|entry_number is used by the BLAST software to identify
the entry, if the user has not provided another identifier.  If the
identifier was parsed, then gnl|BL_ORD_ID|entry_number is replaced by
the correct identifier, as described in
ftp://ncbi.nlm.nih.gov/blast/db/README .

There are no separators between these deflines.

**************************************************************** */

#ifdef USE_MMAP
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#ifdef IBM_AIX
#include <fcntl.h>
#else
#include <sys/fcntl.h>
#endif
#endif

#ifdef USE_MMAP
#ifndef MAP_FILE
#define MAP_FILE 0
#endif
#endif

#ifdef UNIX
#define RBSTR "r"
#else
#define RBSTR "rb"
#endif

#define XTERNAL
#include "uascii.h"

#define XTERNAL
#include "upam.h"
#include "ncbl2_head.h"

#include "defs.h"
#include "mm_file.h"

unsigned int bl2_uint4_cvt(unsigned int);
unsigned int bl2_long4_cvt(long);
int64_t bl2_long8_cvt(int64_t);
void src_int4_read(FILE *fd,  int *valp);
void src_uint4_read(FILE *fd,  unsigned int *valp);
void src_long4_read(FILE *fd,  long *valp);
void src_char_read(FILE *fd,  char *valp);

/* nt_btoa maps  from blast 2bit  format to ascii  characters */
static char nt_btoa[5] = {"ACGT"};

static char aa_b2toa[27]= {"-ABCDEFGHIKLMNPQRSTVWXYZU*"};

static int aa_btof[32];	/* maps to fasta alphabet */

static int dbtype, dbformat, amb_cnt;

#define NCBIBL20 12

int ncbl2_getliba(char *, int, char *, int, FSEEK_T *, int *, struct lmf_str *, long *);
int ncbl2_getlibn(char *, int, char *, int, FSEEK_T *, int *, struct lmf_str *, long *);
void newname(char *, char *, char *, int);
 
void ncbl2_ranlib(char *, int, FSEEK_T, char *, struct lmf_str *m_fd);

/* ncbl2_openlib() is used to open (and memory map) a BLAST2.0 format
   file.  Ifdef USE_MMAP, then ncbl2_openlib returns a structure that can
   be used to read the database. */
   
struct lmf_str *
ncbl2_openlib(char *name, int ldnaseq)
{
  char hname[256];
  char sname[256];
  char tname[256];
  int title_len;
  char *title_str=NULL;
  int date_len;
  char *date_str=NULL;
  long ltmp;
  int tmp;
  int i;
  unsigned long line_len, c_len, clean_count;
#ifdef USE_MMAP
  struct stat statbuf;
#endif
  FILE *ifile;	/* index offsets, also DB info */
  unsigned int *f_pos_arr;
  struct lmf_str *m_fptr;

  if (!ldnaseq) {	/* read a protein database */
    newname(tname,name,AA_INDEX_EXT,(int)sizeof(tname));
    newname(hname,name,AA_HEADER_EXT,(int)sizeof(hname));
    newname(sname,name,AA_SEARCHSEQ_EXT,(int)sizeof(sname));

    /* initialize map of BLAST2 amino acids to FASTA amino acids */
    for (i=0; i<sizeof(aa_b2toa); i++) {
      if ((tmp=aascii[aa_b2toa[i]])<NA) aa_btof[i]=tmp;
      else if (aa_b2toa[i]=='*') aa_btof[i]=aascii['X'];
      else aa_b2toa[i]=0;
/*    else aa_btof[i]=aascii['X']; */
    }
  }
  else {	/* reading DNA library */
    newname(tname,name,NT_INDEX_EXT,(int)sizeof(tname));
    newname(hname,name,NT_HEADER_EXT,(int)sizeof(hname));
    newname(sname,name,NT_SEARCHSEQ_EXT,(int)sizeof(sname));

  }
	
  /* open the index file */
  if ((ifile = fopen(tname,RBSTR))==NULL) {
    fprintf(stderr," cannot open %s (%s) INDEX file",tname,name);
    perror("...");
    return 0;
  }
  src_uint4_read(ifile,(unsigned *)&dbformat); /* get format DB version number */
  src_uint4_read(ifile,(unsigned *)&dbtype);   /* get 1 for protein/0 DNA */

  if (dbformat != FORMATDBV) {
    fprintf(stderr,"error - %s wrong formatdb version (%d/%d)\n",
	    tname,dbformat,FORMATDBV);
    return NULL;
  }

  if ((!ldnaseq && dbtype != AAFORMAT) || (ldnaseq && dbtype!=NTFORMAT)) {
    fprintf(stderr,"error - %s wrong format (%d/%d)\n",
	    tname,dbtype,(ldnaseq ? NTFORMAT: AAFORMAT));
    return NULL;
  }

  /* the files are there - allocate lmf_str */

  if ((m_fptr=(struct lmf_str *)calloc(1,sizeof(struct lmf_str)))==NULL) {
    fprintf(stderr," cannot allocate lmf_str\n");
    return NULL;
  }

  /* open the header file */
  if ((m_fptr->hfile = fopen(hname,RBSTR))==NULL) {
    fprintf(stderr," cannot open %s header file\n",hname);
    goto error_r;
  }

  /* ncbl2_ranlib is used for all BLAST2.0 access */
  m_fptr->ranlib = ncbl2_ranlib;

  if (ldnaseq==1) {
    m_fptr->getlib = ncbl2_getlibn;
    m_fptr->sascii = nascii;
  }
  else {
    m_fptr->getlib = ncbl2_getliba;
    m_fptr->sascii = aascii;
  }
  strncpy(m_fptr->lb_name,sname,MAX_FN);

  /* open the sequence file */

#if defined (USE_MMAP) 
  m_fptr->mm_flg=((m_fptr->mmap_fd=open(sname,O_RDONLY))>=0);
  if (!m_fptr->mm_flg) {
    fprintf(stderr," cannot open %s",sname);
    perror("...");
  }
  else {
    if(fstat(m_fptr->mmap_fd, &statbuf) < 0) {
      fprintf(stderr," cannot fstat %s",sname);
      perror("...");
      m_fptr->mm_flg = 0;
    }
    else {
      m_fptr->st_size = statbuf.st_size;
      if((m_fptr->mmap_base = 
	  mmap(NULL, m_fptr->st_size, PROT_READ,
	       MAP_FILE | MAP_SHARED, m_fptr->mmap_fd, 0)) == (char *) -1) {
	fprintf(stderr," cannot mmap %s",sname);
	perror("...");
	m_fptr->mm_flg = 0;
      }  
      else {
	m_fptr->mmap_addr = m_fptr->mmap_base;
	m_fptr->mm_flg = 1;
      }
    }
    /* regardless, close the open()ed version */
    close(m_fptr->mmap_fd);
  }
#else
  m_fptr->mm_flg = 0;
#endif

  if  (!m_fptr->mm_flg) {
    if ((m_fptr->libf = fopen(sname,RBSTR))==NULL) {
      fprintf(stderr," cannot open %s sequence file",sname);
      perror("...");
      goto error_r;
    }
  }

/* all files should be open */

  src_uint4_read(ifile,(unsigned *)&title_len);

  if (title_len > 0) {
    if ((title_str = calloc((size_t)title_len+1,sizeof(char)))==NULL) {
      fprintf(stderr," cannot allocate title string (%d)\n",title_len);
      goto error_r;
    }
    fread(title_str,(size_t)1,(size_t)title_len,ifile);
  }
  
  src_uint4_read(ifile,(unsigned *)&date_len);

  if (date_len > 0) {
    if ((date_str = calloc((size_t)date_len+1,sizeof(char)))==NULL) {
      fprintf(stderr," cannot allocate date string (%d)\n",date_len);
      goto error_r;
    }
    fread(date_str,(size_t)1,(size_t)date_len,ifile);
  }
  
  m_fptr->lpos = 0;
  src_uint4_read(ifile,(unsigned *)&m_fptr->max_cnt);
  src_long4_read(ifile,&ltmp);
  m_fptr->tot_len = ltmp;
  src_long4_read(ifile,&ltmp);
  m_fptr->max_len = ltmp;

  /* currently we are not using this information, but perhaps later */
  if (title_str!=NULL) free(title_str);
  if (date_str!=NULL) free(date_str);

#ifdef DEBUG
    fprintf(stderr,"%s format: BL2 (%s)  max_cnt: %d, totlen: %ld, maxlen %ld\n",
	    name,m_fptr->mm_flg ? "mmap" : "fopen", 
	    m_fptr->max_cnt,m_fptr->tot_len,m_fptr->max_len);
#endif

  /* allocate and read hdr indexes */
  if ((f_pos_arr=(unsigned int *)calloc((size_t)m_fptr->max_cnt+1,sizeof(int)))==NULL) {
      fprintf(stderr," cannot allocate tmp header pointers\n");
      goto error_r;
    }

  if ((m_fptr->d_pos_arr=(MM_OFF *)calloc((size_t)m_fptr->max_cnt+1,sizeof(MM_OFF)))==NULL) {
      fprintf(stderr," cannot allocate header pointers\n");
      goto error_r;
    }

  /* allocate and read sequence offsets */
  if ((m_fptr->s_pos_arr=(MM_OFF *)calloc((size_t)m_fptr->max_cnt+1,sizeof(MM_OFF)))==NULL) {
      fprintf(stderr," cannot allocate sequence pointers\n");
      goto error_r;
    }

  /*
  for (i=0; i<=m_fptr->max_cnt; i++) src_uint4_read(ifile,&m_fptr->d_pos_arr[i]);
  for (i=0; i<=m_fptr->max_cnt; i++) src_uint4_read(ifile,&m_fptr->s_pos_arr[i]);
  */
  if (fread(f_pos_arr,(size_t)4,m_fptr->max_cnt+1,ifile)!=m_fptr->max_cnt+1) {
    fprintf(stderr," error reading hdr offsets: %s\n",tname);
    goto error_r;
  }

  for (i=0; i<=m_fptr->max_cnt; i++)
#ifdef IS_BIG_ENDIAN
    m_fptr->d_pos_arr[i] = f_pos_arr[i];
#else
    m_fptr->d_pos_arr[i] = bl2_uint4_cvt(f_pos_arr[i]);
#endif

  if (fread(f_pos_arr,(size_t)4,m_fptr->max_cnt+1,ifile)!=m_fptr->max_cnt+1) {
    fprintf(stderr," error reading seq offsets: %s\n",tname);
    goto error_r;
  }
  for (i=0; i<=m_fptr->max_cnt; i++) {
#ifdef IS_BIG_ENDIAN
    m_fptr->s_pos_arr[i] = f_pos_arr[i];
#else
    m_fptr->s_pos_arr[i] = bl2_uint4_cvt(f_pos_arr[i]);
#endif
  }

  if (dbtype == NTFORMAT) {
    /* allocate and ambiguity  offsets */
    if ((m_fptr->a_pos_arr=(MM_OFF *)calloc((size_t)m_fptr->max_cnt+1,sizeof(MM_OFF)))==NULL) {
      fprintf(stderr," cannot allocate sequence pointers\n");
      goto error_r;
    }

    /*
    for (i=0; i<=m_fptr->max_cnt; i++) src_uint4_read(ifile,&m_fptr->a_pos_arr[i]);
    */

    if (fread(f_pos_arr,(size_t)4,m_fptr->max_cnt+1,ifile)!=m_fptr->max_cnt+1) {
      fprintf(stderr," error reading seq offsets: %s\n",tname);
      goto error_r;
    }
    for (i=0; i<=m_fptr->max_cnt; i++) {
#ifdef IS_BIG_ENDIAN
      m_fptr->a_pos_arr[i] = f_pos_arr[i];
#else
      m_fptr->a_pos_arr[i] = bl2_uint4_cvt(f_pos_arr[i]);
#endif
    }
  }

  /*
  for (i=0; i < min(m_fptr->max_cnt,10); i++) {
    fprintf(stderr,"%d: %d %d %d\n",i,m_fptr->s_pos_arr[i],m_fptr->a_pos_arr[i],m_fptr->d_pos_arr[i]);
  }
  */

  /* all done with ifile, close it */
  fclose(ifile);

  free(f_pos_arr);

  if (!m_fptr->mm_flg) {
    tmp = fgetc(m_fptr->libf);
    if (tmp!=NULLB)
      fprintf(stderr," phase error: %d:%d found\n",0,tmp);
  }

  m_fptr->bl_lib_pos = 1;
  amb_cnt = 0;
  return m_fptr;

 error_r:
  /* here if failure after m_fptr allocated */
  free(m_fptr);
  return NULL;
}

void ncbl2_closelib(struct lmf_str *m_fptr)
{
  if (m_fptr->s_pos_arr !=NULL) {
    free(m_fptr->s_pos_arr);
    m_fptr->s_pos_arr = NULL;
  }
  if (m_fptr->a_pos_arr!=NULL) {
    free(m_fptr->a_pos_arr);
    m_fptr->a_pos_arr = NULL;
  }

  if (m_fptr->hfile !=NULL ) {
    fclose(m_fptr->hfile); m_fptr->hfile=NULL;
    free(m_fptr->d_pos_arr); m_fptr->d_pos_arr = NULL;
  }

#ifdef use_mmap
  if (m_fptr->mm_flg) {
    munmap(m_fptr->mmap_base,m_fptr->st_size);
    m_fptr->mmap_fd = -1;
  }
  else 
#endif
    if (m_fptr->libf !=NULL ) {fclose(m_fptr->libf); m_fptr->libf=NULL;}
}

int
ncbl2_getliba(char *seq,
	      int maxs,
	      char *libstr,
	      int n_libstr,
	      FSEEK_T *libpos,
	      int *lcont,
	      struct lmf_str *m_fd,
	      long *l_off)
{
  register char *sptr, *dptr;
  int s_chunk, d_len, lib_cnt;
  long seqcnt;
  long tmp;
  char ch, *bp;
  static long seq_len;
  
  *l_off = 1;

  lib_cnt = m_fd->lpos;
  *libpos = (FSEEK_T)m_fd->lpos;

  if (*lcont==0) {
    if (lib_cnt >= m_fd->max_cnt) return -1;	/* no more sequences */
    seq_len = m_fd->s_pos_arr[lib_cnt+1] - m_fd->s_pos_arr[lib_cnt]; /* value is +1 off to get the NULL */
    if (m_fd->mm_flg) m_fd->mmap_addr = m_fd->mmap_base+m_fd->s_pos_arr[lib_cnt];
#if !defined(DEBUG) && !defined(PCOMPLIB)
    libstr[0]='\0';
#else
    /* get the name from the header file */
    fseek(m_fd->hfile,m_fd->d_pos_arr[lib_cnt],0);

    d_len = min(n_libstr-1,m_fd->d_pos_arr[lib_cnt+1]-m_fd->d_pos_arr[lib_cnt]-1);
    fread(libstr,(size_t)1,(size_t)d_len,m_fd->hfile);
    libstr[d_len]='\0';
#endif
    }

  if (seq_len <= maxs) { /* sequence fits */
    seqcnt = seq_len;
    m_fd->lpos++;
    *lcont = 0;
  }
  else {		/* doesn't fit */
    seqcnt = maxs-1;
    (*lcont)++;
  } 

  if (m_fd->mm_flg) sptr = m_fd->mmap_addr;
  else {
    if ((tmp=fread(seq,(size_t)1,(size_t)seq_len,m_fd->libf))!=(size_t)seq_len) {
      fprintf(stderr," could not read sequence record: %ld %ld != %ld\n",
	      *libpos,tmp,seq_len);
      goto error; 
    }
    sptr = seq;
  }
  if (seq_len <= maxs) {seqcnt = --seq_len;}

  /* everything is ready, set up dst. pointer, seq_len */
  dptr = seq;

  if (aa_b2toa[sptr[seq_len-1]]=='*') seq_len--;
  s_chunk = seqcnt/16;
  while (s_chunk-- > 0) {
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
    *dptr++ = aa_btof[*sptr++];
  }
  while (dptr < seq+seqcnt) *dptr++ = aa_btof[*sptr++];

  if (m_fd->mm_flg) m_fd->mmap_addr = sptr;

  /* we didn't get it all, so reset for more */
  if (*lcont) seq_len -= seqcnt;

  seq[seqcnt]= EOSEQ;
  return (seqcnt);
  
error:	fprintf(stderr," error reading %ld at %ld\n",libstr,*libpos);
  fflush(stderr);
  return (-1);
}

char tmp_amb[4096];

int
ncbl2_getlibn(char *seq,
	      int maxs,
	      char *libstr,
	      int n_libstr,
	      FSEEK_T *libpos,
	      int *lcont,
	      struct lmf_str *m_fd,
	      long *l_off)
{
  register char *sptr, *tptr, stmp;
  long seqcnt;
  int s_chunk, lib_cnt;
  size_t tmp;
  char ch;
  static long seq_len;
  static int c_len,c_pad;
  int c_len_set, d_len;
  char *bp;

  *l_off = 1;

  lib_cnt = m_fd->lpos;
  *libpos = (FSEEK_T)lib_cnt;
  if (*lcont==0) {	/* not a continuation of previous */
    if (lib_cnt >= m_fd->max_cnt) return (-1);
    c_len = m_fd->a_pos_arr[lib_cnt]- m_fd->s_pos_arr[lib_cnt];
    if (!m_fd->mm_flg) {
      if (m_fd->bl_lib_pos != m_fd->s_pos_arr[lib_cnt]) { /* are we positioned to read? */
	amb_cnt++;
	if ((m_fd->bl_lib_pos - m_fd->s_pos_arr[lib_cnt]) < sizeof(tmp_amb)) {
	  /* jump over amb_ray */
	  fread(tmp_amb,(size_t)1,(size_t)(m_fd->s_pos_arr[lib_cnt]-m_fd->bl_lib_pos),m_fd->libf);
	}
	else {	/* fseek over amb_ray */
	  fseek(m_fd->libf,m_fd->s_pos_arr[lib_cnt],0);
	}
	m_fd->bl_lib_pos = m_fd->s_pos_arr[lib_cnt];
      }
    }
    else m_fd->mmap_addr = m_fd->mmap_base + m_fd->s_pos_arr[lib_cnt];
#if !defined(DEBUG) && !defined(PCOMPLIB)
    libstr[0]='\0';
#else
    /* get the name from the header file */
    fseek(m_fd->hfile,m_fd->d_pos_arr[lib_cnt],0);

    d_len = min(n_libstr-1,m_fd->d_pos_arr[lib_cnt+1]-m_fd->d_pos_arr[lib_cnt]-1);
    fread(libstr,(size_t)1,(size_t)d_len,m_fd->hfile);
    libstr[d_len]='\0';
#endif
  }			/* end of *lcont==0 */

  /* To avoid the situation where c_len <= 1; we must anticipate what
     c_len will be after this pass.  If it will be <= 64, back off this
     time so next time it will be > 64 */

  seq_len = c_len*4;

  if ((seq_len+4 > maxs) && (seq_len+4 - maxs  <= 256)) {
    /* we won't be done but we will have less than 256 to go */
    c_len -= 64; seq_len -= 256; c_len_set = 1; maxs -= 256;}
  else c_len_set = 0;

  /*
  fprintf(stderr," lib_cnt: %d %d %d %d\n",lib_cnt,c_len,seq_len,maxs);
  */

  /* does the rest of the sequence fit? */
  if (seq_len <= maxs-4 && !c_len_set) {
    seqcnt = c_len;
    if (!m_fd->mm_flg) {
      if ((tmp=fread(seq,(size_t)1,(size_t)seqcnt,m_fd->libf))!=(size_t)seqcnt) {
	fprintf(stderr,
		" could not read sequence record: %s %ld %ld != %ld: %d\n",
		libstr,*libpos,tmp,seqcnt,*seq);
	goto error; 
      }
      m_fd->bl_lib_pos += tmp;
      sptr = seq + seqcnt;
    }
    else sptr = m_fd->mmap_addr+seqcnt;

    *lcont = 0;		/* this is the last chunk */
    lib_cnt++;		/* increment to the next sequence */
    /* the last byte is either '0' (no remainder) or the last 1-3 chars and the remainder */
    c_pad = *(sptr-1);
    c_pad &= 0x3;	/* get the last (low) 2 bits */
    seq_len -= (4 - c_pad);	/* if the last 2 bits are 0, its a NULL byte */
  }
  else {	/* get the next chunk, but more to come */
    seqcnt = ((maxs+3)/4)-1;
    if (!m_fd->mm_flg) {
      if ((tmp=fread(seq,(size_t)1,(size_t)(seqcnt),m_fd->libf))!=(size_t)(seqcnt)) {
	fprintf(stderr," could not read sequence record: %ld %ld/%ld\n",
		*libpos,tmp,seqcnt);
	goto error;
      }
      m_fd->bl_lib_pos += tmp;
      sptr = seq + seqcnt;
    }
    else {
      sptr = m_fd->mmap_addr+seqcnt;
      m_fd->mmap_addr += seqcnt;
    }
    seq_len = 4*seqcnt;
    c_len -= seqcnt;
    if (c_len_set) {c_len += 64; maxs += 256;}
    (*lcont)++;
/*  hopefully we don't need this because of c_len -= 64. */
/*
    if (c_len == 1) {
#if !defined (USE_MMAP)
      c_pad = fgetc(m_fd->libf);
      *sptr=c_pad;
#else
      c_pad = *m_fd->mmap_addr++;
      sptr = m_fd->mmap_addr;
#endif
      c_pad &= 0x3;
      seq_len += c_pad;
      seqcnt++;
      lib_cnt++;
      *lcont = 0;
    }
*/
  }

  /* point to the last packed byte and to the end of the array
     seqcnt is the exact number of bytes read
     tptr points to the destination, use multiple of 4 to simplify math
     sptr points to the source, note that the last byte will be read 4 cycles
     before it is written
     */
  
  tptr = seq + 4*seqcnt;
  s_chunk = seqcnt/8;
  while (s_chunk-- > 0) {
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
  }
  while (tptr>seq) {
    stmp = *--sptr;
    *--tptr = (stmp&3) +1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
    *--tptr = ((stmp >>= 2)&3)+1;
  }
  /*
    for (sptr=seq; sptr < seq+seq_len; sptr++) {
    printf("%c",nt[*sptr]);
    if ((int)(sptr-seq) % 60 == 59) printf("\n");
    }
    printf("\n");
    */

  m_fd->lpos = lib_cnt;
  if (seqcnt*4 >= seq_len) {	/* there was enough room */
    seq[seq_len]= EOSEQ;
    /* printf("%d\n",seq_len); */
    return seq_len;
  }
  else {				/* not enough room */
    seq[seqcnt*4]=EOSEQ;
    seq_len -= 4*seqcnt;
    return (4*seqcnt);
  }
  
error:	fprintf(stderr," error reading %ld at %ld\n",libstr,*libpos);
  fflush(stderr);
  return (-1);
}

void
ncbl2_ranlib(char *str,
	     int cnt,
	     FSEEK_T libpos,
	     char *libstr,
	     struct lmf_str *m_fd)
{
  int llen, lib_cnt;
  char *bp;

  lib_cnt = (int)libpos;
  llen = m_fd->d_pos_arr[lib_cnt+1]-m_fd->d_pos_arr[lib_cnt];
  if (llen >= cnt) llen = cnt-1;
  fseek(m_fd->hfile,m_fd->d_pos_arr[libpos],0);

  fread(str,(size_t)1,(size_t)(llen),m_fd->hfile);
  str[llen]='\0';

  bp = str;
  while((bp=strchr(bp,'\001'))!=NULL) {*bp++=' ';}

  if (!m_fd->mm_flg) fseek(m_fd->libf,m_fd->s_pos_arr[libpos],0);


  m_fd->lpos = lib_cnt;
  m_fd->bl_lib_pos = m_fd->s_pos_arr[lib_cnt];
}

unsigned int bl2_uint4_cvt(unsigned int val)
{
  unsigned int res;
#ifdef IS_BIG_ENDIAN
  return val;
#else /* it better be LITTLE_ENDIAN */
  res = ((val&255)*256)+ ((val>>8)&255);
  res = (res<<16) + (((val>>16)&255)*256) + ((val>>24)&255);
  return res;
#endif
}  

unsigned int bl2_long4_cvt(long val)
{
  int val4;
  unsigned int res;
#ifdef IS_BIG_ENDIAN
  val4 = val;
  return val4;
#else /* it better be LITTLE_ENDIAN */
  res = ((val&255)*256)+ ((val>>8)&255);
  res = (res<<16) + (((val>>16)&255)*256) + ((val>>24)&255);
  return res;
#endif
}  

int64_t bl2_long8_cvt(int64_t val)
{
  int64_t res;
#ifdef IS_BIG_ENDIAN
  return val;
#else /* it better be LITTLE_ENDIAN */
  res = ((val&255)*256)+ ((val>>8)&255);
  res = (res<<16) + (((val>>16)&255)*256) + ((val>>24)&255);
#ifdef BIG_LIB64
  res = (res<<16) + (((val>>32)&255)*256) + ((val>>40)&255);
  res = (res<<16) + (((val>>48)&255)*256) + ((val>>56)&255);
#else
  fprintf(stderr,"Cannot use bl2_long8_cvt without 64-bit longs\n");
  exit(1);
#endif
  return res;
#endif
}  

void src_uint8_read(FILE *fd, unsigned long *val)
{
#ifdef IS_BIG_ENDIAN
  fread((char *)val,sizeof(long),(size_t)1,fd);
#else
  unsigned char b[8];

  fread((char *)&b[0],(size_t)1,(size_t)8,fd);
  *val = 0;
  *val = (long)((((((long)((long)(b[0]<<8)+(long)b[1]<<8)+(long)b[2]<<8)
		   +((long)b[3]<<8))+((long)b[4]<<8))+((long)b[5]<<8))
		   +((long)b[6]<<8))+(long)b[7];
#endif
}

void src_int4_read(FILE *fd,  int *val)
{
#ifdef IS_BIG_ENDIAN
  fread((char *)val,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *val = 0;
  *val = (int)((int)((int)(b[0]<<8)+(int)b[1]<<8)+(int)b[2]<<8)
	  +(int)b[3];
#endif
}

void src_long4_read(FILE *fd,  long *valp)
{
  int val4;
#ifdef IS_BIG_ENDIAN
  fread(&val4,(size_t)4,(size_t)1,fd);
  *valp = val4;
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  val4 = 0;
  val4 = (int)((int)((int)(b[0]<<8)+(int)b[1]<<8)+(int)b[2]<<8)
	  +(int)b[3];
  *valp = val4;
#endif
}

void src_uint4_read(FILE *fd,  unsigned int *valp)
{
#ifdef IS_BIG_ENDIAN
  fread(valp,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *valp = 0;
  *valp = (unsigned int)((int)((int)(b[0]<<8)+(int)b[1]<<8)+(int)b[2]<<8)
	  +(int)b[3];
#endif
}

void src_long8_read(FILE *fd,  long *val)
{
#ifdef IS_BIG_ENDIAN
  fread((void *)val,(size_t)8,(size_t)1,fd);
#else
  unsigned char b[8];

  fread((char *)&b[0],(size_t)1,(size_t)8,fd);
  *val = 0;
  *val = (long)((((((long)((long)(b[0]<<8)+(long)b[1]<<8)+(long)b[2]<<8)
		  +(long)b[3]<<8)+(long)b[4]<<8)+(long)b[5]<<8)
		+(long)b[6]<<8)+(long)b[7];
#endif
}

void src_char_read(FILE *fd, char *val)
{
  fread(val,(size_t)1,(size_t)1,fd);
}

void src_fstr_read(FILE *fd, char *val,  int slen)
{
  fread(val,(size_t)slen,(size_t)1,fd);
}

void
newname(char *nname, char *oname, char *suff, int maxn)
{
  strncpy(nname,oname,maxn-1);
  strncat(nname,".",1);
  strncat(nname,suff,maxn-strlen(nname));
}
