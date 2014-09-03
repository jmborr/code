#include <stdio.h> /* standart input,output */
#include <ctype.h> /* for charachter recognition */
#include <stdlib.h> /* for conversion from char to dec */
#include <strings.h>
#include <time.h>
#include <math.h>
#include "bcp.h"
#include "controls.h"
#include "fs_output.h"
#include "make_system.h"

enum fs_keys 
{MIN_KEY=0,TEXT_KEY=MIN_KEY,dtext_KEY,STORY_KEY,rate_KEY,
store_KEY,KK_KEY,JJ_KEY,
NN_KEY,MAXKK_KEY,DATA_KEY,TIME_KEY,SAVE_KEY,END_KEY,MAX_KEY};


static int save_pos=0;
static time_t next_t;
static char next_ct_name[60];
static time_t curr_t; 
static double next_ct=0;
static int KK=0;
static int JJ=0;
static int NN=1;
static int MAXKK=1000;
static double rate=1;
static double store=2;
static double dtext=600;
static char TEXT_NAME[80]="bcp_restart";
static char STORY_NAME[80]="bcp_story";
static char DATA_NAME[80]="bcp_data";
static char SAVE_NAME[80]="bcp_data";
static double init_time=0;
static int good=0;
static char ** keywords;
static double maxtime;
static int col_num=0;
static int n_atom;
static int n_dim;
static double param[4];
static int n_mes=0;
static int KK_length;
static int NN_length;

extern int fs_ok(void)
{
  return good;
}


void fs_add(double dat0,double dat1,double dat2,double dat3)
{
  param[0]+=dat0;
  param[1]+=dat1;
  param[2]+=dat2;
  param[3]+=dat3;
  n_mes++;
  return;
}

void init_fs_param(void)
{
  param[0]=0;
  param[1]=0;
  param[2]=0;
  param[3]=0;
  n_mes=0;
  return;
}

double next_update(void)
{
int jj=JJ+1;
int kk=KK;
  if (jj>NN){kk++;jj=0;}
return 	rate*(pow(store,(double)NN)*kk+pow(store,(double)jj));
}

void write_story1()
{
  FILE * path;
  char new_name[100];
  int fErr=0;
  int i,k;
  double timeb=get_timeb();
  double newtime=get_time();
  double timec=get_timec();
  dimensions * bound=get_bounds();
  moved_iatom * a=(moved_iatom *)get_atom();
  advance_timeb(); 
  moveatoms();
  set_timeb(timeb);
  sprintf(new_name,"%s-%0*d-%0*d",STORY_NAME,KK_length,KK,NN_length,JJ);
  printf("%s\n",new_name);
  path=fopen(new_name,"w");
  if(!path)return;
  if(n_mes)
    for(i=0;i<4;i++)
      param[i]/=n_mes;
  else
    {
      param[0]=get_temperature();
      param[1]=-get_avePot();
      param[2]=param[1]+param[0]*n_atom*n_dim*0.5;
      param[3]=get_pressure();
    }
  
  fprintf(path,"%lf %lf\n",newtime,timec);
  fprintf(path,"%d\n",n_atom);
  fprintf(path,"%lf %lf %lf\n",bound[0].length,bound[1].length,
	  bound[2].length);
  fprintf(path,"%lf %lf %lf\n",param[0]*n_atom*n_dim*0.5,param[1],param[2]);
  fprintf(path,"%le %lf %lf %lf\n",param[3],param[0],get_coeff(),
	  get_temp_limit());
  init_fs_param();
  for(i=0;i<n_atom;i++)
    {
      for(k=0;k<n_dim;k++)
	fprintf(path,"%lf ",a[i].r[k]);
	fprintf(path,"\n");
    }
  corr_vel();
  for(i=0;i<n_atom;i++)
    {
      for(k=0;k<n_dim;k++)
	fprintf(path,"%lf ",a[i].u[k]);
	fprintf(path,"\n");
    }
  for(i=0;i<n_atom;i++)
    fprintf(path,"%d\n",a[i].c);
 
  fflush(path);
  fclose(path);  
  sprintf(new_name,"gzip %s-%0*d-%0*d",STORY_NAME,KK_length,KK,NN_length,JJ);
  system(new_name);
  return;		

}


void write_story2()
{
  FILE * path;
  char new_name[100];
  int fErr=0;
  int i,k;
  double timeb=get_timeb();
  double newtime=get_time();
  double timec=get_timec();
  dimensions * bound=get_bounds();
  moved_iatom * a=(moved_iatom *)get_atom();
  advance_timeb(); 
  moveatoms();
  set_timeb(timeb);
  sprintf(new_name,"%s-%s",STORY_NAME,next_ct_name);
  printf("%s\n",new_name);
  path=fopen(new_name,"w");
  if(!path)return;
  if(n_mes)
    for(i=0;i<4;i++)
      param[i]/=n_mes;
  else
    {
      param[0]=get_temperature();
      param[1]=-get_avePot();
      param[2]=param[1]+param[0]*n_atom*n_dim*0.5;
      param[3]=get_pressure();
    }
  
  fprintf(path,"%lf %lf\n",newtime,timec);
  fprintf(path,"%d\n",n_atom);
  fprintf(path,"%lf %lf %lf\n",bound[0].length,bound[1].length,
	  bound[2].length);
  fprintf(path,"%lf %lf %lf\n",param[0]*n_atom*n_dim*0.5,param[1],param[2]);
  fprintf(path,"%le %lf %lf %lf\n",param[3],param[0],get_coeff(),
	  get_temp_limit());
  init_fs_param();
  for(i=0;i<n_atom;i++)
    {
      for(k=0;k<n_dim;k++)
	fprintf(path,"%lf ",a[i].r[k]);
	fprintf(path,"\n");
    }
  corr_vel();
  for(i=0;i<n_atom;i++)
    {
      for(k=0;k<n_dim;k++)
	fprintf(path,"%lf ",a[i].u[k]);
	fprintf(path,"\n");
    }
  for(i=0;i<n_atom;i++)
    fprintf(path,"%d\n",a[i].c);
 
  fflush(path);
  fclose(path);  
  sprintf(new_name,"gzip %s-%s",STORY_NAME,next_ct_name);
  system(new_name);
  return;		

}


void save_data(void)
{
  double current_time=get_time();
  FILE * path =fopen(DATA_NAME,"w");
  if(path)
    {
      fprintf(path,"%-6s %s\n",keywords[TEXT_KEY],TEXT_NAME);
      fprintf(path,"%-6s %lf\n",keywords[dtext_KEY],dtext);
      fprintf(path,"%-6s %s\n",keywords[STORY_KEY],STORY_NAME);
      if(good==1)
	{      
	  fprintf(path,"%-6s %lf\n",keywords[rate_KEY],rate);
	  fprintf(path,"%-6s %lf\n",keywords[store_KEY],store);
	  fprintf(path,"%-6s %d\n",keywords[KK_KEY],KK);
	  fprintf(path,"%-6s %d\n",keywords[JJ_KEY],JJ);
	  fprintf(path,"%-6s %d\n",keywords[NN_KEY],NN);
	  fprintf(path,"%-6s %d\n",keywords[MAXKK_KEY],MAXKK);
	}
      if(good==2)
	fprintf(path,"%-6s %s\n",keywords[SAVE_KEY],SAVE_NAME);
      fprintf(path,"%-6s %s\n",keywords[DATA_KEY],DATA_NAME);
      fprintf(path,"%-6s %lf\n",keywords[TIME_KEY],current_time);
      fprintf(path,"%-6s\n",keywords[END_KEY]);
      fflush(path);
      fclose(path);
    }
}

extern void fs_close(void)
{
fs_output();
  good=0;
}


/*all these keywords will show up in the DATA file*/
int init_fs_keywords()
{ 
  int i;   
  keywords=(char **)malloc(MAX_KEY*sizeof(unsigned char *));
  if(!keywords) return 0;
  keywords[TEXT_KEY]="TEXT";
  keywords[dtext_KEY]="dtext";
  keywords[STORY_KEY]="STORY";
  keywords[rate_KEY]="rate";
  keywords[store_KEY]="store";
  keywords[KK_KEY]="KK";
  keywords[JJ_KEY]="JJ";
  keywords[NN_KEY]="NN";
  keywords[MAXKK_KEY]="MAXKK";
  keywords[DATA_KEY]="DATA";
  keywords[TIME_KEY]="TIME";
  keywords[SAVE_KEY]="SAVE";
  keywords[END_KEY]="END";
  return 1;
}

extern int fs_init(void)
{
  FILE * infile;
  FILE * save=NULL;
  char name[100];
  char value[60];
  double dummy;   

  int i;
  good=0;
  
  printf("we do not record story\n");
  if(yes())return good;

  n_atom= get_atom_number(); /*number of atoms*/
  n_dim=get_dimension();
  
  printf("What is DATA file name?\n");
  scanf("%s",name);
  infile=fopen(name,"r");
  if(!infile)return good;
  
  if(!init_fs_keywords())return good;
  init_fs_param();
  while(!feof(infile))
    {
      char keyword[60];
      enum fs_keys key;
      FILE *try;
      fscanf(infile,"%s%s",keyword,value);
      for(key=MIN_KEY;key<MAX_KEY;key++)
	if(!strcmp(keyword,keywords[key]))
	  {
	    switch(key){
	    case TEXT_KEY:{
	      try=fopen(value,"a");
	      if(try){
		fclose(try);
		strcpy(TEXT_NAME,value);}
	      break;
	    }
            case STORY_KEY:{
	      try=fopen(value,"a");
	      if(try){
		fclose(try);
		strcpy(STORY_NAME,value);}
	      break;
	    }
            case DATA_KEY:{
	      try=fopen(value,"a");
	      if(try){
		fclose(try);
		strcpy(DATA_NAME,value);}
	      break;
	    }
            case SAVE_KEY:{
	      save=fopen(value,"r");
	      if(save)strcpy(SAVE_NAME,value);
	      break;
	    }
  
	    case TIME_KEY:{dummy=atof(value);init_time=dummy;break;}
	    case dtext_KEY:{dummy=atof(value);if(dummy>=1)dtext=dummy;break;}
            case rate_KEY:{dummy=atof(value);if(dummy>0)rate=dummy;break;}
	    case store_KEY:{dummy=atof(value);if(dummy>1)store=dummy;break;}
            case KK_KEY:{dummy=atof(value);if(dummy>=-1)KK=dummy;break;} 
            case JJ_KEY:{dummy=atof(value);if(dummy>=0)JJ=dummy;break;} 
            case NN_KEY:{dummy=atof(value);if(dummy>=0)NN=dummy;break;} 
            case MAXKK_KEY:{dummy=atof(value);if(dummy>=1)MAXKK=dummy;break;}
            case END_KEY: goto finish;   
	    }
	    break;
	  }
    }
 finish:
  if(save)
    {
      while(!feof(save))
      {
	fscanf(save,"%s",value);
	dummy=atof(value);
	if(dummy>init_time)break;
      }
      if(dummy>init_time)
	{ 
	  next_ct=dummy;
	  strcpy(next_ct_name,value); 
	  set_time(init_time);
	  save_pos=ftell(save);
	  maxtime=next_ct;
	  while(!feof(save))
	    {
	      fscanf(save,"%s",value);
	      dummy=atof(value);
	      if(maxtime<dummy)maxtime=dummy;
	    }
	  printf("Maxtime=%lf\n",maxtime);
	  col_num=n_atom-1;
	  time(&next_t);
	  n_mes=0;
	  good=2;
	  set_max_time(maxtime);
	}
      fclose(save); 
    }
  if(good!=2)
    {
      if(JJ>NN){JJ=NN;printf("JJ>NN\n");}
      if(KK>MAXKK){MAXKK=KK;printf("KK>MAXKK\n");}
      next_ct=next_update();
      if(next_ct<init_time){
	printf("TT does not correspond to KK and JJ\n");
	init_time=next_ct;
      } 
      set_time(init_time);
      maxtime=rate*(pow(store,(double)NN)*MAXKK);
      printf("Maxtime=%lf\n",maxtime);
      n_mes=0;
      KK_length=0;
      i=MAXKK;
      while(i)
	{
	  i/=10;
	  KK_length++;
	}
      i=NN;
      while(i)
	{
	  i/=10;
	  NN_length++;
	}
      good=1;
      set_max_time(maxtime);
    }

  return good; 
}

double read_save(void)
{
  char value[60];
  double dummy=next_ct;
  double curr_ct;
  double newmaxtime;
  FILE * save=fopen(SAVE_NAME,"r");
  if(save)
    {
      fseek(save,save_pos,SEEK_SET);
      while(!feof(save))
	{
	  fscanf(save,"%s",value);
	  dummy=atof(value);
          if(dummy>next_ct)break;
	}
      if(dummy>next_ct)
	{     	
          curr_ct=dummy;
          newmaxtime=dummy;
          strcpy(next_ct_name,value);
	  save_pos=ftell(save);
	  while(!feof(save))
	    {
	      fscanf(save,"%s",value);
	      dummy=atof(value);
	      if(newmaxtime<dummy)newmaxtime=dummy;
	    }
	  if(newmaxtime !=maxtime)
	    {
	    maxtime=newmaxtime;
	    set_max_time(maxtime);
	  fclose(save);
	  return curr_ct;
	} 
    }
  maxtime=0;
  set_max_time(maxtime);
  fclose(save);
  good=0;
  return 0;
    }
}



void fs_output(void)
{
if(good)
  {
    double new_time=get_time();
    int save_restart=0; 
    col_num++;
    if(col_num==n_atom)
      {
	time(&curr_t);
	col_num=0;    
        if(curr_t>=next_t){
	  save_restart=1;
	  next_t=curr_t+dtext;
	}

      }
    if(new_time>=next_ct)
      {
	if(good==1)
	  {
	    JJ++;
	    if(JJ>NN)
	      {
		KK++;
		JJ=0;
	      }  
	    printf("%lf %lf %d %d\n",new_time,next_ct,JJ,KK);
	    save_restart=1;
	    write_story1();
	    next_ct=next_update();
	  }
	if(good==2)
	  {
	    printf("%lf %lf\n",new_time,next_ct);
	    save_restart=1;
	    write_story2();
	    next_ct=read_save();
	  }
      }
    if(save_restart)
      {
	FILE * path=fopen(TEXT_NAME,"wb");
	if(path)
	  {
            int res;
	    double timeb=get_timeb();
            time(&curr_t); 
	    next_t=curr_t+dtext;
	    advance_timeb();
	    res=write_key_coord(path);
            fflush(path);
            fclose(path); 
	    set_timeb(timeb);
	    if(res==noErr)save_data();
	  }
      }

  }      
}
