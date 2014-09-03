using namespace std;
#include<iostream>
#include<cstdio>
#include<fstream>
#include<algorithm>
#include<cstring>
#include "allocate.h"
#include "rmsd.h"

/******************* PSEUDOCODE *******************
1.create temporary file with byte offsets for the header of each snapshot
2.record for each snapshot the byte offset for the header and the RMSD score of the snapshot
3.go bin by bin in the RMSD score span
  3.1.collect all snapshots withinn the bin
  3.2.shuffle the collected snapshots
  3.3.filter the snapshots according to rmsd cut-off
  3.4.output filtered snapshots to the filtered file
  3.5.go next RMSD score bin
***************************************************/

/* global variables */
int maxnu=200000; /*maximum number of unfiltered snapshots*/
int maxnf=1000; /*maximum number of filtered snapshots per bin*/

typedef struct snapshot{
  int offset;   /*byte offset position in the trajectory file where the snapshot begins*/
  double rmsd;    /*rmsd score to native*/
  char *h;      /*header line*/
  double **xyz; /*coordinates*/
}snapshot;

/*===============================================================*/
/*go to the byte offset stated by sn.offset in the trajectory file and
  load header line onto char *line, and coordinates onto double **xyx */
void load_coords(snapshot *sn, FILE *pin, double **&xyz, char *&line,int L){
  /*printf("xyz=%d line=%d\n",xyz,line);*/
  if(!line){
    //printf("empty line\n");
    line=new char[264];   /*printf("sn->offset=%d\n",sn->offset);*/
  }
  if(!xyz){
    //printf("empty xyz\n");
    xyz=alloc2<double>(L,3);
  }
  fseek(pin,sn->offset,SEEK_SET); /*printf("ftell=%d\n",ftell(pin));exit(1);*/
  fgets(line,264,pin);  /*printf("%s",line); exit(1);*/
  for(int i=0;i<L;i++) fscanf(pin,"%lf %lf %lf",&xyz[i][0],&xyz[i][1],&xyz[i][2]);
}
/*===============================================================*/
void fill_coords(snapshot *sn, double **xyz, char *line,int L){
  if(!sn->h){
    //printf("empty sn->h\n");
    sn->h=new char[264];
  }
  if(!sn->xyz){
    //printf("empty sn->xyz\n");
    sn->xyz=alloc2<double>(L,3);
  }
  sprintf(sn->h,"%s",line);
  for(int i=0;i<3*L;i++) sn->xyz[0][i]=xyz[0][i];
}
/*===============================================================*/
/*dump header and coordinates of snapshot onto char *lines */
void snapshotToChar(snapshot *sn,char *lines,int L){
  char *line=new char[128];
  double **xyz=sn->xyz;
  sprintf(lines,"%s",sn->h); /*printf("%s",lines);exit(1);*/
  for(int i=0;i<L;i++){
    sprintf(line,"%10.3lf %10.3lf %10.3lf\n",xyz[i][0],xyz[i][1],xyz[i][2]);
    strcat(lines,line); /*append line to lines*/
  }
}
/*===============================================================*/
/*deallocate header and coordinates*/
void dealloc_info(snapshot **f, int nfbin){
  for(int s=0;s<nfbin;s++){
    delete [] (f[s]->h);
    dealloc2(f[s]->xyz); /*deallocate a double ** variable */
  }
}
/*===============================================================*/
/*welcoming message*/
int message(){
  system("clear");
  printf("Usage: ./filterByRmsd.x [options]\n");
  printf("  Optional:\n");
  printf("  -a rep1RgF.tra input file with snapshots  (default=\"rep1RgF.tra\")\n");
  printf("  -b rep1RMSD.tra output file with filtered snapshots (def=\"./rep1rmsdF.tra\")\n");
  printf("  -c minRMSD:  minimum RMSD (def=2.0)\n");
  printf("  -d maxRMSD:  maximum TM-score (def=8.0)\n");
  printf("  -e nbin:   number of bins (def=6)\n");
  printf("  -f nf:     maximum number of filtered snapshots per bin (def=1000)\n");
  printf("  -g rmsdco: rmsd cut-off when filtering (def=1.0Angstroms)\n");
  printf("  -i histf:  file name for histrogram of RMSD scores (def: rmsdHisto.dat)\n");
  printf("\n");
  return 1 ; /*failure!*/
}
/*======================================================================*/
/*error handling messages*/
typedef enum {_MNF_,_MUF_} errflag;
int error(errflag f){
  switch(f){
  case _MNF_:
    cout<<"ERROR: maximum number of filtered snaphots per bin is "<<maxnf<<endl;
    break;
  case _MUF_:
    cout<<"ERROR: maximum number of unfiltered snaphots is "<<maxnu<<endl;
    break;
  }
  return 1;
}
/*======================================================================*/
/*assume the header line of inpf is like:
L=100 E= -1112.6 n=     1 rep=11 cycle=    1 rmsd=0.716......
*/
int getSeqLength(char *inpf){
  int L;
  char *cmd=new char[264];
  sprintf(cmd,"head -1 %s|tr -s ' '|cut -d' ' -f 2",inpf);
  FILE *fp=popen(cmd, "r");
  fscanf(fp,"%d",&L);
  delete []cmd;
  return L;
}
/*======================================================================*/
/*execution starts here*/
int main(int argc, char **argv){
  
  /*bunch of variables*/
  int L;        /*sequence length*/
  int N;        /*number of snapshots in the input file with the unfiltered snapshots*/
  int nu;       /*number of unfiltered snapshots*/
  int nf=maxnf; /*maximum number of filtered snapshots per bin*/
  char *inpf=NULL;  /*input file with unfiltered snapshots*/
  char *outf=NULL;  /*output file with filtered snapshots*/
  char *histf=NULL; /*file name for histrogram of TM scores*/
  double minRMSD=2.0; /*default minimum TM score*/
  double maxRMSD=8.0; /*default maximum TM score*/
  int nbin=6;        /*default number of TM score bins*/
  double drmsd;        /*increase in RMSD score when shifting by one bin*/
  double rmsdco=2.5; /*rmsd cut-off when filtering*/

  {/*BLOCK TO PARSE INPUT*/
    int option_char,f=0;
    extern char *optarg;
    extern int optind, optopt;
    while ((option_char = getopt(argc, argv, ":a:b:c:d:e:f:g:i:h")) != EOF){
      switch (option_char){  
      case 'a': inpf=optarg; break;
      case 'b': outf=optarg;  break;
      case 'c': minRMSD=atof(optarg); break;
      case 'd': maxRMSD=atof(optarg); break;
      case 'e': nbin=atoi(optarg); break;
      case 'f': nf=atoi(optarg); break;
      case 'g': rmsdco=atof(optarg); break;
      case 'i': histf=optarg; break;
      case 'h': return message(); 
      default: return message();
      }
    }
    if(f<0){ return message(); }
    if(nf>maxnf) return error(_MNF_);
    
    /*initialize input if not passed as arguments*/
    if(!inpf){ inpf=new char[14];  sprintf(inpf,"rep1RgF.tra"); }
    if(!outf){ outf=new char[14];  sprintf(outf,"rep1rmsdF.tra"); }
    if(!histf){histf=new char[14]; sprintf(histf,"rmsdHisto.dat"); }
    /*printf("%s %s %lf %lf %d %d %lf\n",inpf,outf,minRMSD,maxRMSD,nbin,nf,rmsdco);exit(1);*/
  }/*end of block to parse input*/

  L=getSeqLength(inpf); /*sequence length*/ /*printf("L=%d\n",L);exit(1);*/

  /*write byte offsets for each header*/
  char *cmd=new char[264];
  sprintf(cmd,"grep -b rmsd %s > junk.filterByRmsd",inpf); system(cmd); /*exit(1);*/
  
  /*total number of unfiltered snapshots*/
  {
    FILE *fp=popen( "wc -l junk.filterByRmsd", "r");
    fscanf(fp,"%d",&nu);
    if(nu>maxnu) return error(_MUF_);
    fclose(fp);    /*printf("nu=%d\n",nu);exit(1);*/
  }

  /*record byte offset and TM score of each snapshot*/
  snapshot *sns=new snapshot[nu];
  {
    for(int i=0;i<nu;i++){ sns[i].h=NULL; sns[i].xyz=NULL; } /*initialize*/
    sprintf(cmd,"cat junk.filterByRmsd | cut -d':' -f 1"); /*printf("%s\n",cmd);exit(1);*/
    FILE *fp=popen(cmd,"r");    
    for(int i=0;i<nu;i++) fscanf(fp,"%d",&(sns[i].offset)); /*byte offsets*/
    fclose(fp);
    sprintf(cmd,"cat junk.filterByRmsd|cut -d '=' -f 2");
    /*printf("%s\n",cmd);exit(1);*/
    fp=popen(cmd,"r");    
    for(int i=0;i<nu;i++) fscanf(fp,"%lf",&(sns[i].rmsd));/*RMSD scores*/
    fclose(fp);
    /*for(int i=0;i<nu;i++) printf("of=%d rmsd=%lf\n",sns[i].offset,sns[i].rmsd); exit(1);*/
  }

  FILE *pout=fopen(outf,"w");
  FILE *pin=fopen(inpf,"r");
  FILE *ph=fopen(histf,"w");
  fprintf(ph," RMSD N(RMSD)\n");

  snapshot **unf=new snapshot*[nu]; /*pointers to unfiltered snapshots within the RMSD bin*/
  snapshot **f=new snapshot*[nu];   /*pointers to filtered snapshots within the RMSD bin*/
  int nubin;                        /*actual number of unfiltered snapshots in the RMSD bin*/
  int nfbin;                        /*actual number of filtered snapshots in the RMSD bin*/
  drmsd=(maxRMSD-minRMSD)/nbin;     /*printf("drmsd=%lf\n",drmsd);exit(1);*/
  double rmsdbegin=minRMSD;
  double rmsdend=rmsdbegin+drmsd;
  int bin;
  double rmsd_sn;
  double **xyz=alloc2<double>(L,3); /*stores coordinates for one snapshot*/
  bool nonanalogous;
  char *line=new char[256];
  char *lines=new char[128*L];
  snapshot *sn;
  double rmsd,minrmsd;
  for(bin=0; bin<nbin; bin++){
    rmsdco=rmsdbegin/2.0;
    /*printf("bin=%d rmsdbegin=%lf rmsdend=%lf\n",bin,rmsdbegin,rmsdend);*/
    /*collect all snapshots with rmsd within the RMSD bin*/
    sn=sns;  /*pointer to current snapshot info*/
    nubin=0; /*number of unfiltered bins with TM score in [tmbegin,tmend)*/
    for(int s=0;s<nu;s++){
      rmsd_sn=sn->rmsd; //printf("rmsd=%lf rmsdbegin=%lf rmsdend=%lf\n",rmsd_sn,rmsdbegin,rmsdend);
      if( rmsd_sn>=rmsdbegin and rmsd_sn<rmsdend) unf[nubin++]=sn; /*point to the snapshot info*/
      sn++; /*go to next snapshot info*/
    }
    printf("bin=%d rmsdbegin=%lf rmsdend=%lf nubin=%d\n",bin,rmsdbegin,rmsdend,nubin);
    if(!nubin){ /*no snapshots in this RMSD score bin. Thus, go to next bin*/
      fprintf(ph, "%5.2lf %5d\n",(rmsdbegin+rmsdend)/2,0); /*output histogram of rmsd scores*/
      rmsdbegin=rmsdend ; rmsdend+=drmsd ;
      continue;
    }
    /*randomize the snapshots*/
    /*for(int i=0;i<nubin;i++) printf("of=%d\n",unf[i]->offset); printf("***********\n");*/
    random_shuffle(unf, unf+nubin); /*using random shuffling of "algorithm" library*/
    /*for(int i=0;i<nubin;i++) printf("of=%d\n",unf[i]->offset);exit(1);*/

    /*obtain nf (or less) filtered structures from the unfiltered structures in the RMSD bin*/
    f[0]=unf[0];                             /*initialize f with one snapshot*/
    load_coords(f[0],pin,f[0]->xyz,f[0]->h,L); /*read coords and header line from .tra file*/
    snapshotToChar(f[0],lines,L); /*dump header and coordinates onto variable lines*/
    fprintf(pout,"%s",lines);

    /*for(int i=0;i<nubin;i++) printf("of=%d\n",unf[i]->offset);exit(1);*/
    /*snapshotToChar(f[0],lines,L);printf("%s",lines);exit(1);*/
    nfbin=1;                                 /*current number of filtered structures*/
    for(int s=1;s<nubin;s++){ /*go through all remaining unfiltered snapshots in the RMSD bin*/
      nonanalogous=true;
      load_coords(unf[s],pin,xyz,line,L); /*load coords of unfiltered snapshot onto temp. xyz*/
      minrmsd=1000.0;
      /*compare to filtered structures. Add as filtered if not structurally similar*/
      for(int i=0;i<nfbin;i++){
	rmsd=getrmsd(xyz,f[i]->xyz,L); /*printf("rmsd=%lf\n",rmsd);exit(1)*/
	if(rmsd<minrmsd) minrmsd=rmsd;
	if( rmsd<rmsdco ){
	  nonanalogous=false; /*current snapshot is analogous to one filtered snapshot*/
	  /*printf("analogous rmsd=%lf\n",rmsd); exit(1);*/
	  break;
	}
      }

      if(nonanalogous){
	f[nfbin]=unf[s]; /*adding a filtered snapshot*/
	fill_coords(f[nfbin],xyz,line,L);
	/*output filtered structures*/
	snapshotToChar(f[nfbin-1],lines,L); /*dump header and coordinates onto variable lines*/
	fprintf(pout,"%s",lines);
	nfbin++;
      }
      
      if(nfbin==nf) break; /*don't collect more than nf snapshots per bin*/
      /*printf("bin=%2d s=%5d nfbin=%5d minrmsd=%5.2lf\n",bin,s,nfbin,minrmsd);*/
    }
    fprintf(ph, "%5.2lf %5d\n",(rmsdbegin+rmsdend)/2,nfbin); /*output histogram of rmsd scores*/
    /*reclaim memory space by deallocating header and coordinates of filtered snapshots*/
    dealloc_info(f,nfbin);   /*exit(1);*/
    rmsdbegin=rmsdend ; rmsdend+=drmsd ; /*shift to next bin by update of bin boundaries*/
  }

  fclose(pout);
  fclose(pin);
  fclose(ph);
  system("/bin/rm junk.filterByRmsd"); /*some clean up*/
  return 0;
}/*end of main(..)*/