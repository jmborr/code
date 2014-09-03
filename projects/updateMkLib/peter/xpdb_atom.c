
/*
 *   XPDB -> Extract a set of PDB sequences
 *   Written by Piotr Rotkiewicz, 2004
 *   piotr@pirx.com
 */                                           

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>

#include "zlib.h"
 
#define DEBUG 

#define SEQID_CUTOFF 0.35
#define LENGTH_CUTOFF 0.75

#define MAX_SEQ_LENGTH 3000
#define MIN_SEQ_LENGTH 10

#define UPDATE_FRQ 100

char AA_NAMES[21][4] =  
  { "GLY", "ALA", "SER", "CYS", "VAL",
    "THR", "ILE", "PRO", "MET", "ASP",
    "ASN", "LEU", "LYS", "GLU", "GLN",
    "ARG", "HIS", "PHE", "TYR", "TRP",
    "UNK" };

char SHORT_AA_NAMES[22] = { "GASCVTIPMDNLKEQRHFYWX" };

int AA_NUMS[256];

typedef struct _cluster_type {
  char pdbid[6];
  struct _cluster_type *next;
} cluster_type;

typedef struct _database_type {
  struct _database_type *prev;
  cluster_type *cluster;
  int length;
  char *seq;
  char *name;
  char *expdta;
  char *resol;
  struct _database_type *next;
} database_type;

database_type *database;

float **SCORE;
float **VAL;
char **MOV;

float GAP_OPEN_PENALTY=-11.0;
float GAP_EXTN_PENALTY=-1.0;

FILE *outlog, *outfasta;

float align(int height, int width)
{
  register int i,j;
  float score;
  register float h,v,d;

    for (i=1; i<=height; i++) {
      VAL[i][0] = 0;
      MOV[i][0] = 0;
    }

    for (j=1; j<=width; j++) {
      VAL[0][j] = 0;
      MOV[0][j] = 0;
    }

    VAL[0][0] = 0;
    MOV[0][0] = 0;

    for (i=1;i<height+1;i++)
      for (j=1;j<width+1;j++) {
        h=VAL[i-1][j];
        if (MOV[i-1][j]!=1)
          h+=GAP_OPEN_PENALTY;  
        else
          h+=GAP_EXTN_PENALTY;  
        v=VAL[i][j-1];
        if (MOV[i][j-1]!=2)
          v+=GAP_OPEN_PENALTY;
        else
          v+=GAP_EXTN_PENALTY;  
        d=VAL[i-1][j-1]+SCORE[i-1][j-1];
        if (d>=h && d>=v) {
          VAL[i][j]=d;
          MOV[i][j]=0;
        } else
        if (h>=d && h>=v) {
          VAL[i][j]=h;
          MOV[i][j]=1;
        } else {
          VAL[i][j]=v;
          MOV[i][j]=2;
        }
      }

    score = VAL[height][width];

  return score;
}


float get_seqid(int height, int width, char *seq1, char *seq2)
{
  int i,j,nali,id,al;
  char c;

    i=height;
    j=width;
    id=al=0;
/*    
    fprintf(outlog,"\n");

    fprintf(outlog,"%d %d\n", height, width);
    for (i=0;i<height;i++)
      fprintf(outlog,"%c",seq1[i]);
    fprintf(outlog,"\n");
            
    for (i=0;i<width;i++)
      fprintf(outlog,"%c",seq2[i]);
    fprintf(outlog,"\n");
            
    i=height;
    j=width;
    id=al=0;
*/    
    do {
      switch (MOV[i][j]) {
        case 0:
//	  fprintf(log,"%3d %3d %c%c ", i, j, seq1[i-1], seq2[j-1]);
          if (seq1[i-1]==seq2[j-1]) {
	    id++;
//	    fprintf(log,"*\n");
	  } //else fprintf(log,"\n");    
	  al++;
          i--;
          j--;
        break;
        case 1:
//	  fprintf(log,"%3d %3d %c%c\n", i, j, seq1[i-1], '-');
          i--;
        break;
        case 2:
//	  fprintf(log,"%3d %3d %c%c\n", i, j, '-', seq2[j-1]);
          j--;
        break;
      }
    } while (i>0 && j>0);

    if (height<width) nali=height; else nali=width;
            
//    fprintf(log,"\nHEIGHT: %d WIDTH: %d SEQID: %d (%6.2f%%) NALI: %d\n", height, width, id, 100.*(float)id/(float)al, nali);
    
  return (float)id/(float)nali;
}


float MUT_MATRIX[21][21];

void read_mutation_matrix(char *name)
{
  FILE *input;
  int i, j;
  char *pc;
  int mmcodes[20];
  char line[200];
  float val;

    for (i=0; i<21; i++)
      for (j=0; j<21; j++)
        MUT_MATRIX[i][j] = 0.0;

    input = fopen(name, "rt");
    if (input) {
      do {
        fgets(line, 200, input);
        if (feof(input)) {
          exit(-1);
        }
      } while (line[0]!='I');
      pc = strchr(line, '=')+2;
      for (i=0; i<20; i++)
        mmcodes[i] = AA_NUMS[(int)pc[i]];
      for (i=0; i<20; i++)
        for (j=0; j<=i; j++) {
          fscanf(input, "%f", &val);
          MUT_MATRIX[mmcodes[i]][mmcodes[j]]=val;
          MUT_MATRIX[mmcodes[j]][mmcodes[i]]=val;
    }
      fclose(input);
    } else {
      exit(-1);
    }
}


int read_mutation_matrix_blast(char *name)
{
  FILE *input;
  int i, j;
  char *pc, dum[10];
  int mmcodes[21];
  char line[200];
  float val;

    for (i=0; i<21; i++)
      for (j=0; j<21; j++)
        MUT_MATRIX[i][j] = 0;

    input = fopen(name, "rt");
    if (input) {
      do {
        fgets(line, 200, input);
        if (feof(input)) {
          fclose(input);
          return 0;
        }
      } while (!strstr(line,"A  R  N  D"));
      pc = strchr(line, 'A');
      for (i=0; i<20; i++) {
        mmcodes[i] = AA_NUMS[(int)(*pc)];
        pc+=3;
      }  
      for (i=0; i<20; i++) {
        fscanf(input,"%s",dum);
        for (j=0; j<20; j++) {
          fscanf(input,"%f",&val);
          MUT_MATRIX[mmcodes[i]][mmcodes[j]]=(char)val;
          MUT_MATRIX[mmcodes[j]][mmcodes[i]]=(char)val;
        }  
        fgets(line,200,input);
      }
      fclose(input);
/*      
      for (i=0;i<20;i++) {
        for (j=0;j<20;j++)
          printf("%2d ", MUT_MATRIX[i][j]);
        printf("\n");
      }
*/          
      return 1;
    } else {
      return 0;
    }
}


database_type* check_seqid(char *seq, int len)
{
  database_type *tmpdata;
  int store, i, j;
  float ratio, seqid;
  
    if (!database) return NULL;
    
    tmpdata=database;

    while (tmpdata->next) 
      tmpdata=tmpdata->next;
      
    while (tmpdata) {
      ratio=(float)len/(float)tmpdata->length;
      if (ratio>1.) ratio=1./ratio;
      if (ratio>LENGTH_CUTOFF) { // check seq id
        for (i=0;i<len;i++)
          for (j=0;j<tmpdata->length;j++)
            SCORE[i][j]=MUT_MATRIX[AA_NUMS[seq[i]]][AA_NUMS[tmpdata->seq[j]]]; 
        align(len, tmpdata->length);
        seqid=get_seqid(len, tmpdata->length, seq, tmpdata->seq);        
        if (seqid>SEQID_CUTOFF) {
//          fprintf(log,"GOT: %5.2f, ITER: %d\n", seqid, iter);
          return tmpdata;                                  
        }  
      }  
      tmpdata=tmpdata->prev;
    } 

  return NULL; // new sequence - store!
}

int check_pdbid(char *pdbid)
{
  database_type *tmpdata;
  cluster_type *tmpcluster;
  
    tmpdata=database;
    while (tmpdata) {
      tmpcluster=tmpdata->cluster;
      while (tmpcluster) {     
        if (!strncmp(pdbid,tmpcluster->pdbid,4)) return 1;
        tmpcluster=tmpcluster->next;
      }  
      tmpdata=tmpdata->next;
    } 

  return 0; 
}

int check_pdbid_cid(char *pdbid)
{
  database_type *tmpdata;
  cluster_type *tmpcluster;
  
    tmpdata=database;
    while (tmpdata) {
      tmpcluster=tmpdata->cluster;
      while (tmpcluster) {     
        if (!strncmp(pdbid,tmpcluster->pdbid,5)) return 1;
        tmpcluster=tmpcluster->next;
      }  
      tmpdata=tmpdata->next;
    } 

  return 0; 
}

char setseq(char* aaname)
{
  int i;

    for (i=0; i<21; i++)
      if ((aaname[0]==AA_NAMES[i][0]) &&
          (aaname[1]==AA_NAMES[i][1]) &&
          (aaname[2]==AA_NAMES[i][2]))
         break;
      if (i==21) {
        if (!strncmp(aaname, "GLX", 3))
          return 'E';
        if (!strncmp(aaname, "HID", 3))
          return 'H';
        if (!strncmp(aaname, "CYX", 3))
          return 'C';
        if (!strncmp(aaname, "ASX", 3))
          return 'D';
        i--;
    }

  return SHORT_AA_NAMES[i];
}

int N_DATA=0;

void save_database(char *name)
{
  FILE *out;
  database_type *tmpdata;
  cluster_type *tmpcluster;
  
    out=fopen(name,"w");
    tmpdata=database;
    while (tmpdata) {
      tmpcluster=tmpdata->cluster;
      fprintf(out,">");
      while (tmpcluster) {
        fprintf(out," %5s", tmpcluster->pdbid);
        tmpcluster=tmpcluster->next;
      }             
      fprintf(out,"\n");
      fprintf(out,"%s\n", tmpdata->seq);
      fprintf(out,"\n");
      tmpdata=tmpdata->next;
    }    
    fclose(out);
}

void save_database_fasta(char *name)
{
  FILE *out;
  database_type *tmpdata;
  cluster_type *tmpcluster;
  
    out=fopen(name,"w");
    tmpdata=database;
    while (tmpdata) {
      tmpcluster=tmpdata->cluster;
      fprintf(out,">%5s %s\n", tmpcluster->pdbid, tmpdata->name);
      fprintf(out,"%s\n", tmpdata->seq);
      fprintf(out,"\n");
      tmpdata=tmpdata->next;
    }    
    fclose(out);
}

void save_clusters(char *name)
{
  FILE *out;
  database_type *tmpdata;
  cluster_type *tmpcluster;
  
    out=fopen(name,"w");
    tmpdata=database;
    while (tmpdata) {
      tmpcluster=tmpdata->cluster;
      fprintf(out,">");
      while (tmpcluster) {
        fprintf(out,"%5s", tmpcluster->pdbid);
        tmpcluster=tmpcluster->next;
        if (tmpcluster) fprintf(out," "); else fprintf(out,"\n\n");
      }             
      tmpdata=tmpdata->next;
    }    
    fclose(out);
}

void process(char *fullname, char *name)
{   
//  FILE *inp;
  char buf[1000];
  char pdbid[6];
  char expdta[5];
  char lname[2];
  char ca_label[2][20]= {" ", "CA_ONLY!" };
  char rname[100], prevrnum[100], rnum[100];
  char *cmp_name[100], *cmp_chain[100], *ptr;
  float res;
  int n_cmp_names, n_cmp_chains, unknown;
  char *tmp, *tmp2;
  database_type *newdata, *tmpdata;
  cluster_type *newcluster, *tmpcluster;
  int i, j, ok, length, reallength, model, seqres;
  char cid, rid, prevrid, prevcid, aa, *aaname;
  int resnum, prevresnum, nucleic, ca_only;
  char seq[10*MAX_SEQ_LENGTH];
  gzFile *inp;
  
    pdbid[0]=0; pdbid[4]=0; pdbid[5]=0;
    expdta[5]=0;
    res = -1.0;
    
    sprintf(expdta,"    ");
    
    fprintf(outlog,"-----------------------------------\nPROCESSING %s\n", name);
    
    inp=gzopen(fullname,"rb");
    if (!inp) {
      fprintf(outlog,"ERROR: File not found\n");
      return;
    }

    tmpdata=database;
    if (tmpdata)
      while (tmpdata->next) tmpdata=tmpdata->next;
    
// extract a pdb code

    while (!gzeof(inp)) {
      if (gzgets(inp,buf,1000)==buf) {
        if (!strncmp(buf,"HEADER",6) && strlen(buf)>65) {                                                              
          if (buf[62]>='0' && buf[62]<='9' && buf[63]!=' ' && buf[64]!=' ' && buf[65]!=' ') {
            strncpy(pdbid,&buf[62],5);
            break;
          }  
        } else {
          if (buf[72]>='0' && buf[72]<='9' && buf[73]!=' ' && buf[74]!=' ' && buf[75]!=' ') {
            strncpy(pdbid,&buf[72],5);
            break;
          } 
        }   
      }  
    }

    cmp_name[0]=NULL;
    n_cmp_names=n_cmp_chains=0;

    length=0;
    unknown=0;
    ca_only=0;
    nucleic=0;
        
    if (pdbid[0]) {

#ifdef DEBUG
fprintf(outlog,"Found PDBID : %s\n", pdbid);
#endif

      for (i=0;i<4;i++)
        if (pdbid[i]>='A' && pdbid[i]<='Z') pdbid[i]+='a'-'A';
    
// process the file
      gzrewind(inp);

      for (i=0;i<100;i++) {
        cmp_name[i]=NULL;
        cmp_chain[i]=NULL;
      }  
            
        strcpy(prevrnum,"aaaa");

//      if (!check_pdbid(pdbid)) 
//      if (1) 
        while (!gzeof(inp)) {
          if (gzgets(inp,buf,1000)==buf) {
            if (!strncmp(buf,"COMPND",6)) { // extract name of the protein
              if (!strstr(buf,"MOL_ID")) { // old type of compound record
                buf[71]=0;
                cmp_name[n_cmp_names]=(char*)malloc(63);
                strcpy(cmp_name[n_cmp_names],&buf[10]);
                do {
                  if (gzgets(inp,buf,1000)==buf && !strncmp(buf,"COMPND",6)) {
                    buf[71]=0;
                    tmp=(char*)malloc(63);
                    strcpy(tmp,&buf[10]);
                    tmp2=(char*)malloc(strlen(cmp_name[n_cmp_names])+strlen(tmp)+1);
                    strncpy(tmp2,cmp_name[n_cmp_names],strlen(cmp_name[n_cmp_names]));
                    i=strlen(cmp_name[n_cmp_names])-1;
                    while (cmp_name[n_cmp_names][i]==' ') cmp_name[n_cmp_names][i--]=0;
                    cmp_name[n_cmp_names][i]=' ';
                    strcpy(&tmp2[i+1],tmp);
                    free(tmp);
                    free(cmp_name[n_cmp_names]);
                    cmp_name[n_cmp_names]=tmp2;
                    i=strlen(cmp_name[n_cmp_names])-1;
                    while (cmp_name[n_cmp_names][i]==' ') cmp_name[n_cmp_names][i--]=0;
                  }          
                } while (!strncmp(buf,"COMPND",6) && !gzeof(inp));
                n_cmp_names=1;
              } else { // new type: MOL_ID, MOLECULE, CHAIN
                do {
                  if (gzgets(inp,buf,1000)==buf && !strncmp(buf,"COMPND",6)) {
                    tmp=strstr(buf,"MOLECULE:");
                    if (tmp) {
                      tmp+=10;
                      tmp2=strstr(tmp,";");
                      if (tmp2) tmp2[0]=0;
                      if (tmp[strlen(tmp)-1]=='\n') tmp[strlen(tmp)-1]=0;
                      cmp_name[n_cmp_names]=(char*)malloc(strlen(tmp)+1);
                      strcpy(cmp_name[n_cmp_names],tmp);
                      n_cmp_names++;
                    }
                    tmp=strstr(buf,"CHAIN:");
                    if (tmp) {
                      tmp+=7;
                      tmp2=strstr(tmp,";");
                      if (tmp2) tmp2[0]=0;
                      if (tmp[strlen(tmp)-1]=='\n') tmp[strlen(tmp)-1]=0;
                      cmp_chain[n_cmp_chains]=(char*)malloc(strlen(tmp)+1);
                      strcpy(cmp_chain[n_cmp_chains],tmp);
                      n_cmp_chains++;
                    }
                  }
                } while (!strncmp(buf,"COMPND",6) && !gzeof(inp));
              }
              
            } // COMPND
              
            if (!strncmp(buf,"EXPDTA",6)) {
              if (strstr(buf,"X-RAY"))
                sprintf(expdta,"XRAY");
              if (strstr(buf,"NMR"))
                sprintf(expdta," NMR");
              if (strstr(buf,"THEORETICAL"))
                sprintf(expdta,"THEO");
              if (strstr(buf,"NEUTRON"))
                sprintf(expdta,"NEUT");
              if (strstr(buf,"FLUORESCENCE"))
                sprintf(expdta,"FLUO");
              if (strstr(buf,"ELECTRON"))
                sprintf(expdta,"ELEC");
              if (strstr(buf,"FIBER"))
                sprintf(expdta,"FIBR");
            } // EXPDTA
            
    				if (!strncmp(buf,"REMARK",6)) {
    				  if (strstr(buf,"RESOLUTION") && strstr(buf,"ANGSTROMS")) {
    				    ptr = strstr(buf,"RESOLUTION");
    				    sscanf(ptr+11,"%f",&res);
    				  }  
    				  if (expdta[2]==' ' && (strstr(buf,"CRYSTALLOGRAPHIC") || strstr(buf,"X-RAY") || strstr(buf,"XPLOR")))
                sprintf(expdta,"XRAY");    				  
    				} // REMARK
/**    
            if (!strncmp(buf,"SEQRES",6)) {
              length=0;
              ok=0;

              do {
                sscanf(&buf[6],"%d",&resnum);
                cid=buf[11];
                if (!length) {
                  prevresnum=resnum;
                  prevcid=cid;
                }  
                sscanf(&buf[12],"%d",&reallength);
                aaname=&buf[19];
                for (i=0;i<13;i++) {
                  aa=setseq(aaname);
                  seq[length++]=aa;
                  if (aa!='X') ok++; 
                  aaname+=4;
                  if (length>=reallength) break;
                }  
                
                seq[length]=0;
                
                gzgets(inp,buf,1000);
                                
                sscanf(&buf[6],"%d",&resnum);
                cid=buf[11];
                
                if ((gzeof(inp) || strncmp(buf,"SEQRES",6) || 
                     cid!=prevcid || resnum<=prevresnum) && length) { // store the sequence
                  if (prevcid==' ') prevcid='_';
                  pdbid[4]=prevcid;
                  if (ok && length>=MIN_SEQ_LENGTH) {
                    fprintf(outlog,"ID: %s", pdbid);
                    fprintf(outlog,"   LENGTH: %4d", length);
                    fprintf(outlog,"   REAL LENGTH: %4d", reallength);
                    if (length!=reallength) fprintf(outlog, "\nWARNING: DIFFERENT LENGTHS!!!\n");

                    i=0; j=0;
                    if (n_cmp_names>1) {
                      for (i=0;i<n_cmp_chains;i++) {
                        if (strchr(cmp_chain[i],prevcid)) {
                          j=1;
                          break;
                        }  
                      }  
                    }
                    if (!j) i=0;
                    fprintf(outfasta,"> %s %s %5.2f %s\n", pdbid, expdta, res, cmp_name[i]);
                    fprintf(outfasta,"%s\n\n", seq);
                  }    
                  prevresnum=resnum;
                  prevcid=cid;
                  length=0;
                  ok=0;
                  fprintf(outlog,"\n");
                }      
              } while (!strncmp(buf,"SEQRES",6) && !gzeof(inp));
            } // SEQRES
**/

            if (!strncmp(buf,"ATOM",4)) {
              strncpy(rnum,&buf[22],5);
              rnum[5]=0;
              cid=buf[21];              
              if (cid==' ') pdbid[4]='_'; else pdbid[4]=cid;
              if (strcmp(rnum,prevrnum)) {
                if (lname[0]=='C' && lname[1]=='A' && buf[13]=='C' && buf[14]=='A') {
fprintf(outlog,"CA_ONLY: %s %s %s\n", rnum, rname, pdbid);              
                  ca_only=1;
                }  
                lname[0]=buf[13];
                lname[1]=buf[14];
                strcpy(prevrnum,rnum);
                strncpy(rname,&buf[17],3); 
                rname[3]=0;              
                if (!strncmp(rname,"  A",3)||!strncmp(rname,"  T",3)||!strncmp(rname,"  G",3)||!strncmp(rname,"  C",3)) {
                  nucleic=1;
fprintf(outlog,"NUCLEIC: %s %s %s (%d)\n", rnum, rname, pdbid, nucleic);              
                }
                seq[length] = setseq(rname);              
                if (seq[length]=='X') unknown++;
                length++;
              }  
            }

            if (!strncmp(buf,"END",3)||!strncmp(buf,"TER",3)) {
              if (length && !nucleic && length>unknown) {                
                fprintf(outlog,"ID: %s (N: %d)", pdbid, nucleic);
                fprintf(outlog,"   LENGTH: %4d", length);
                fprintf(outlog,"\n");
                i=0; j=0;
                if (n_cmp_names>1) {
                  for (i=0;i<n_cmp_chains;i++) {
                    if (strchr(cmp_chain[i],prevcid)) {
                      j=1;
                      break;
                    }  
                  }  
                }
                if (!j) i=0;                
                fprintf(outfasta,"> %s %s %5.2f %s %s\n", pdbid, expdta, res, ca_label[ca_only], cmp_name[i]);                
                seq[length]=0;
                fprintf(outfasta,"%s\n\n", seq);
              }  
              strcpy(prevrnum,"aaaa");
              unknown=0;
              nucleic=0;
              length = 0;
              ca_only=0;
            }

            if (!strncmp(buf,"ENDMDL",6)) { // just a single model
              break;
            }
          }  
        }
    }

    for (i=0;i<n_cmp_names;i++) {
      if (cmp_name[i]) free(cmp_name[i]);
      if (cmp_chain[i]) free(cmp_chain[i]);
    }  
    
    gzclose(inp);  
    
}

void read_database(char *dataname, char *fastaname)
{
  FILE *inp, *fasta;
  char buf[10000], buf2[10000];
  char *line, *pdbid, *name;
  database_type *newdata, *tmpdata;
  cluster_type *tmpcluster, *newcluster, *startcluster;
  int len, namelen, i;
  
    database=newdata=tmpdata=NULL;
    name=NULL;
    
    fasta=fopen(fastaname,"r");    
    inp=fopen(dataname,"r");       
    
    if (inp) {
      while (!feof(inp)) {
        if (fgets(buf,10000,inp)==buf) {
          if (buf[0]=='>') {
            if (fasta && !feof(fasta)) {

              do {
                if (fgets(buf2,10000,fasta)==buf2) {
                  if (buf2[0]=='>') { // get name from FASTA file
                    i=strlen(buf2)-1;
                    while (i>0 && buf2[i]=='\n') buf2[i--]=0;
                    namelen=strlen(&buf2[7]);
                    if (namelen>1) {
                      name=(char*)malloc(namelen+1);
                      if (name) strcpy(name,&buf2[7]);
                    }  
                    break;
                  }
                }
              } while (!name && !feof(fasta));                        

            }
            newdata=(database_type*)malloc(sizeof(database_type));
            if (name) {
              newdata->name=name;
              name=NULL;
            }  
            newcluster=startcluster=NULL;
            line=&buf[1];
            while (pdbid=(char*)strtok(line," ")) {
              if (pdbid[0]!='\n') {
                tmpcluster=(cluster_type*)malloc(sizeof(cluster_type));
                strncpy(tmpcluster->pdbid,pdbid,6);
                tmpcluster->pdbid[5]=0;
                if (newcluster) newcluster->next=tmpcluster;
                newcluster=tmpcluster;
                if (!startcluster) startcluster=tmpcluster;
              }  
              line=NULL;
            }
            newdata->cluster=startcluster;
          } else 
          if (strlen(buf)>2) {
            len=strlen(buf)-1;
            while (buf[len]<'A' || buf[len]>'Z') {
              len--;         
            }
            len++;
            buf[len]=0;  
            if (newdata) {
              newdata->seq=(char*)malloc(len+1);
              newdata->length=len;
              strncpy(newdata->seq,buf,len+1);
              newdata->seq[len]=0;
              if (!database) database=newdata;
              if (tmpdata) tmpdata->next=newdata;
              newdata->prev=tmpdata;
              tmpdata=newdata;
              newdata=NULL;
            }  
          } 
        }  
      }  
      fclose(inp);
      if (fasta) fclose(fasta);
    }  
}

int check_file(char *name)
{ 
  FILE *inp;
  DIR *dp;          
  
    dp=opendir(name);
    if (dp) {
      closedir(dp);
      return 0;
    }
    
    inp=fopen(name,"r");
    if (inp) {
       fclose(inp);
      return 1;
    }
  return 0;  
}

int scandir_pdb(char *dirname)
{ 
  DIR *dp;          
  struct dirent *ep;
  char buf[512], buf2[512];

    if (dirname[strlen(dirname)-1]!='/') 
      sprintf(buf,"%s/",dirname);
    else  
      sprintf(buf,"%s",dirname);
    dp=opendir(buf);
    if (!dp) return 0;
    while (ep=readdir(dp)) {
      sprintf(buf2,"%s%s",buf,ep->d_name);
      if (check_file(buf2) && strcmp(ep->d_name,".") && strcmp(ep->d_name,"..")) { // check file
        process(buf2, ep->d_name);
      } else { // check if directory
        if (strcmp(ep->d_name,".") && strcmp(ep->d_name,"..")) {
          sprintf(buf2,"%s%s/",buf,ep->d_name);
          scandir_pdb(buf2);  
        }  
      }
    }  
    closedir(dp);
    
  return 1;
}
            
int main(int argc, char **argv)
{          
  int i;
  
    setbuf(stdout,0);
    
    if (argc<2) {
      printf("Use: pdbref <database_path>\n");
      return 0;
    }
    
    outlog=fopen("XPDB_ATOM.LOG","w");
    setbuf(outlog,0);
    
    outfasta=fopen("XPDB_ATOM.FASTA","w");
    setbuf(outfasta,0);
    
    SCORE=(float**)malloc(sizeof(float*)*MAX_SEQ_LENGTH);
    for (i=0;i<MAX_SEQ_LENGTH;i++)
      SCORE[i]=(float*)malloc(sizeof(float)*MAX_SEQ_LENGTH);

    VAL=(float**)malloc(sizeof(float*)*(MAX_SEQ_LENGTH+1));
    for (i=0;i<MAX_SEQ_LENGTH+1;i++)
      VAL[i]=(float*)malloc(sizeof(float)*(MAX_SEQ_LENGTH+1));

    MOV=(char**)malloc(sizeof(char*)*(MAX_SEQ_LENGTH+1));
    for (i=0;i<MAX_SEQ_LENGTH+1;i++)
      MOV[i]=(char*)malloc(sizeof(char)*(MAX_SEQ_LENGTH+1));

    for (i=0; i<255; i++)
      AA_NUMS[i] = 20; // dummy aa code 
    for (i=0; i<20; i++)
      AA_NUMS[(int)SHORT_AA_NAMES[i]] = i;

    scandir_pdb(argv[1]);
 
    fprintf(outlog, "\nScanning done.\n");
    
    fclose(outlog);
      
  return 0;
}
