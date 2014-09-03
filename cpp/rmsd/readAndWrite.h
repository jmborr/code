#include <string.h>
/*======================================================================*/
int readline(FILE *in, char *buf){
  int n=0;
  char c;
  buf[0]='\0';
  c=getc(in);
  if( c==EOF || c=='\n') return 0; /*empty line or end of file*/
  else buf[n++]=c;
  do{ 
    c=getc(in);
    buf[n++]=c; 
  }while(c!='\n');
  if(c=='\n') buf[--n]='\0'; /*overwrite new-line character*/
  return n;
}
/*=======================================================================*/
int read_pdb(char *name, double **r, char *type){
  FILE *in=fopen(name,"r");
  char buf[100],prev[30],next[50];
  int n=0;
  double x,y,z;
  if(!in){
    printf("ERROR: can't open file %s\n",name);
    return 0; /*failure*/
  }

  while(readline(in,buf)){
    /*printf("%s\n",buf);*/
    if((char*)strstr(buf,"TER")==buf) return n;
    if((char*)strstr(buf,"ATOM")==buf){
      if(!type or (char*)strstr(buf,type)!=NULL){
	sscanf(buf+30,"%lf %lf %lf",&r[n][0],&r[n][1],&r[n][2]);
	/*printf("%8.3lf%8.3lf%8.3lf\n",r[n][0],r[n][1],r[n][2]);*/
	n++;
      }
    }
  }
  fclose(in);
  return n; /*length of sequence*/
}
/*=======================================================================*/
int write_pdb(char *outf, char *nat, double **r){
  int i,n=0,l;
  char buf[124];
  FILE *out=fopen(outf,"w");
  FILE *in=fopen(nat,"r");
  if(!in){
    printf("ERROR: can't open file %s\n",nat);
    return 0; /*unsuccessful opening of file*/
  }

  while(readline(in,buf)){
    if((char*)strstr(buf,"TER")==buf){
      fprintf(out,"TER\n");
      fclose(in);
      fclose(out);
      return n;
    }
    if((char*)strstr(buf,"ATOM")==buf && (char*)strstr(buf," CA ")!=NULL){
      l=strlen(buf);
      for(i=1;i<=30;i++) fprintf(out,"%c",buf[i]);
      fprintf(out,"%8.3f%8.3f%8.3f",r[n][0],r[n][1],r[n][2]);
      for(i=55;i<=l;i++) fprintf(out,"%c",buf[i]);
      fprintf(out,"\n");
      n++;
    }
  }

  fclose(in);
  fclose(out);
  return n;
}
/*=======================================================================*/

