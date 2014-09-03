#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

static short n_frames;
static int frame_length;
static int version=0;
static short nat;
static short n_atom;
static short n_param=0;
static char ** param_names;
static double ** r, **r1;
static double * param_values;
static double box[3];
static double scale_factor[3];
static double dt;
static short * atom_colors;
static short n_dim;
static double * atom_radii;
static short * atom_types;
static char *s;
static double * dummy;
static short s_dummy;
static int bo=0;

long header_length;
int read_header(FILE * path);
int read_frame(FILE * path,double ** r);
char * read_java(char * s,void * input,int nbyte);

int byte_order(void)
{ 
  unsigned short i=1;
  unsigned char * s;
  s=(char*)&i;
  //printf("size of short: %ld\n", sizeof(short));
  if(s[1]>s[0])return 1;
  return 0;
}

/*encapsulated fread*/
void fread_bo(void* inpt, int size, int nitem, FILE* path, int bo)
{
  int iread,i;
  static unsigned char buf[100];
  void* top = inpt;
  for(iread=0; iread<nitem; iread++){
    fread(buf, size, 1, path);
    if(bo){//big_eiden
      for(i=0; i<size; i++)
	*((char*)top+i) = buf[i];
    }
    else{//small_eiden
      for(i=0; i<size; i++)
	*((char*)top+i) = buf[size-1-i];
    }
    top = (char*)top + size;
  }
}

int read_header(FILE * movie_file)
{
  char * password="BCPkino";
  char * password1="BCPkina";
  int l_pass=strlen(password);
  int l_param;
  int i,j;
  char * s;
  s=(char *)malloc(1000);
  fread(s,sizeof(char),l_pass,movie_file);
  s[l_pass]=(char)0;
  version=0;
  if(strcmp(password,s)){
    version=1;
    if(strcmp(password1,s))
      {
	printf("wrong movie format\n");
	return 0;
      }
  }
  
  fread_bo(&n_frames,sizeof(n_frames),1,movie_file,bo);
  printf("n_frames=%ld\n",n_frames);
  /* system size components */
  fread_bo(box,sizeof(double),3,movie_file,bo);
  printf("X=%lf Y=%lf Z=%lf\n",box[0],box[1],box[2]);
  /* delta T */
  fread_bo(&dt,sizeof(double),1,movie_file,bo);
  printf("dt=%lf\n",dt);
  /* number of diff. types */
  fread_bo(&nat,sizeof(nat),1,movie_file, bo);
  printf("nat=%d\n",nat);
  
  atom_colors=(short *)malloc(nat*sizeof(short));
  atom_radii=(double *)malloc(nat*sizeof(double));
   /* types of atoms */
  fread_bo(atom_colors,sizeof(short),nat,movie_file,bo);
  /* radius */
  fread_bo(atom_radii,sizeof(double),nat,movie_file,bo);
  /*number of atoms*/
  fread_bo(&n_atom,sizeof(n_atom),1,movie_file,bo);
  printf("n_atom=%d\n",n_atom);
  
  r=(double **)malloc(sizeof(double *)*n_atom);
  r[0]=(double *)malloc(sizeof(double)*n_atom*3);
  for(i=1;i<n_atom;i++)
    r[i]=r[i-1]+3;
  
  /*list of atoms types */
  atom_types=(short *)malloc(n_atom*sizeof(short));
  fread_bo(atom_types,sizeof(short),n_atom,movie_file,bo);  
  /* list of bonds , skip it!*/
  /*---------------FORMAT----------------
    0, ..., 0 (... means the index of friends that are beger than i)
    1, ..., 1
    ...
    ...
    ---------------END-------------------*/
  for(i=0;i<n_atom;i++)
    {
      fread_bo(&s_dummy,sizeof(s_dummy),1,movie_file,bo);
      do{
	fread_bo(&s_dummy,sizeof(s_dummy),1,movie_file,bo);
      }while(s_dummy!=i);
    }
  /* number of params */ 
  fread_bo(&n_param,sizeof(n_param),1,movie_file,bo);
  printf("n_param=%d\n",n_param);
  param_values=(double*)malloc(n_param*sizeof(double));
  /*list of Parameters*/
  for(i=0; i<n_param; i++){
    fread_bo(&s_dummy, sizeof(s_dummy),1,movie_file, bo);
    fread(s, s_dummy, 1, movie_file);
    s[s_dummy]='\0';
    printf("%s\n",s);
  }
  free(s);
  if(box[2]!=0)n_dim=3;
  else n_dim=2;
  for(i=0;i<n_dim;i++)
    scale_factor[i]=box[i]/(double)65536;
  frame_length=(2+n_atom*n_dim)*sizeof(short);
  return ftell(movie_file);
}

int read_frame(FILE * movie_file,double ** r){
  if(movie_file)
    {
      int i,j,k,n_bytes;
      short caption_length=0;
      unsigned short coord;
      char c_dummy;
      short s_dummy;
      /*new list of atoms types ?*/
      fread_bo(&c_dummy, sizeof(c_dummy),1, movie_file,bo);
      if(c_dummy){/*read the new type list*/
	fread_bo(atom_types,sizeof(short),n_atom,movie_file,bo);  
      }
      /*new list of bonds ?*/
      fread_bo(&c_dummy, sizeof(c_dummy),1, movie_file,bo);
      if(c_dummy){
	for(i=0;i<n_atom;i++)
	  {
	    fread_bo(&s_dummy,sizeof(s_dummy),1,movie_file,bo);
	    do{
	      fread_bo(&s_dummy,sizeof(s_dummy),1,movie_file,bo);
	    }while(s_dummy!=i);
	  }
      }
      /* coords */
      if(r){
	if(version)
	  for(i=0;i<n_atom;i++)
	    for(j=0;j<n_dim;j++){
	      fread_bo(&coord, sizeof(coord),1, movie_file, bo);
	      r[i][j]=(coord)*scale_factor[j];
	    }
	else
	  for(i=0;i<n_atom;i++)
	    for(j=0;j<n_dim;j++)
	      {
		fread_bo(&coord, sizeof(coord),1, movie_file, bo);
		r[i][j]=coord;
	      }
      }
      else{
	int nbyte = n_atom*n_dim*sizeof(short);
	fread(s, 1, nbyte, movie_file);
      }
      
      /*params */
      fread_bo(param_values,sizeof(double), n_param,movie_file, bo);
      for(i=0; i<n_param; i++){
	printf("param: %lf\n", param_values[i]);
      }
      
      /*END of FRAME*/
      do{
	fread_bo(&s_dummy,sizeof(s_dummy), 1 ,movie_file, bo);
      }while(s_dummy!=0);
      
      return 1;
    }
  return 0;
}

int main()
{ 
  int i,j,at,n_atom1;
  int delta_t,nfr1;
  long pos;
  FILE *path,*path1,*fp;
  char fname[80];
  bo = byte_order();
  printf("what is movie file name\n");
  scanf("%s",fname);
  path=fopen(fname,"rb");
  if(!path)
    return 0;
  if(!(header_length=read_header(path))){fclose (path);return;}
  
  printf("Version: %d\n", version);
  s=(char*)malloc(header_length>frame_length?header_length:frame_length);	    
  printf("what is out file name ?\n");
  scanf("%s",fname);
  
  printf("number of frames skip ?\n");
  scanf("%ld",&delta_t); 
  printf("number of frames read ?\n");
  scanf("%ld",&nfr1);
  if(!n_frames)
    {
      pos=ftell(path);
      while(read_frame(path,0))
	n_frames++;
      fseek(path,pos,SEEK_SET);
    }
  if(nfr1>n_frames-delta_t)nfr1=n_frames-delta_t;
  if(nfr1<=0)return;
  /* if succesfful, skipping delta_t frames in both files */	
  for(i=0;i<delta_t;i++)
    {
      if(!read_frame(path,r))return;
    }
  
  
  
  
  {
    double dd,maxdd=0;
    double box2[3];
    int nfr=0;
    int yes;
    int res;
    fp=fopen(fname, "w");    
/*    fprintf(fp,"%hd %hd %lf %lf %lf\n",n_atom,nat,box[0],box[1],box[2]);

      for(i=0;i<n_atom;i++)
	{
	  fprintf(fp,"%hd %lf %hd\n",atom_types[i],atom_radii[atom_types[i]-1],i);
	}
*/
    for(i=0;i<3;i++)
      box2[i]=box[i]/2;

    for(j=0;j<nfr1;j++)
      {
	if(yes=read_frame(path,r))
	  {
	/*    fprintf(fp,"\n"); */ 
	    for(i=0;i<n_atom;i++)	
	      fprintf(fp,"%lf %lf %lf\n",r[i][0], r[i][1],r[i][2]);
	  }
      }
  }
  fclose(fp);      
}
