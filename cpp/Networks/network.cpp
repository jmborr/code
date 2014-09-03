#include "network.h"

/*================================================================*/
void   dump_clusters(int **c, int *cs, int nc){
  for(int i=0;i<nc;i++){
    cout<<"cs["<<i<<"]="<<cs[i]<<endl;
    for(int j=0;j<cs[i];j++){
      cout<<c[i][j]<<" ";
    }
    cout<<endl;
  }
}
/*================================================================*/
/* cm:contact map, d:dimension, c:clusters, cs:clusters sizes, nc: number of clusters, f:forest, ff: forest focus, bc: burning cluster, isb: is burning?, bf: burning focus, bl: burning last, n:number of clusters, tid: tree ID*/
void burn_forest_clusters( int **cm, int d, int**&c, int*&cs, int &nc){
  int *f, ff, fid, *bc, *isb, bf, bl, n, tid;
  nc=0;
  f =new int [d+1]; /*indexes of "f" refer to the tree ID, hence "d+1"*/
  for(int i=1;i<=d;i++){f[i]=1;} /*"1" means a green tree (not burnt)*/
  bc=new int [d];
  isb=new int [d+1];
  c =new int*[d];
  cs=new int [d]; for(int i=0;i<d;i++){cs[i]=0;} /*At most, "d" clusters*/
  ff=1; /*we always start the fire by igniting the tree with ID==1*/
  /*cout<<"d="<<d<<endl;*/
  do{/*while there are green trees*/
    for(int i=0;i<d;i++){bc[i]=0;}/*initialize the burning cluster*/
    for(int i=1;i<=d;i++){isb[i]=0;}/*no burning trees yet*/
    tid=ff; /*tree ID set to the tree we are igniting*/
    bc[0]=tid; /*first tree to burn within the cluster*/
    isb[tid]=1;
    bf=0;/*tree that is burning*/
    bl=1;
    while(bf<bl){/*fill the cluster as we go on burning neighbors*/
      /*cout<<"tid="<<tid<<" bf="<<bf<<" bl="<<bl<<endl;*/
      for(int i=1; i<=d; i++){ 
	if(cm[tid][i]!=0 && isb[i]==0){/*"i" is neighbor tree and not burning*/
	  bc[bl]=i;bl++; /*add neighboring green trees to the burning cluster*/
	  isb[i]=1;
	}
      }
      bf++; tid=bc[bf];/*go on to burn next neighbor*/
    }
    nc++;
    cs[nc-1]=bl; /*size of the cluster number "nc".Note shift in array index*/
    c[nc-1]=new int[bl]; /*Note shift in array index*/
    for(int i=0;i<bl;i++){
      /*cout<<"bc["<<i<<"]="<<bc[i]<<endl;*/
      c[nc-1][i]=bc[i];/*fill the cluster with the tree ID;s*/
      f[bc[i]]=0; /*in the forest list, set the burnt tree as burnt*/
    }
    while(f[ff]==0 && ff<=d){ff++;} /*go to next green tree we can find*/
    /*cout<<"ff="<<ff<<endl;*/
  }while(ff<=d);
  /*cout<<"nc="<<nc<<endl;
    dump_clusters(c,cs,nc);*/
}

/*int main(){
  int **cm, d=5, **c, *cs, nc;

  cm=new int*[6]; 
  cm[0]=new int[6]; for(int i=0;i<36;i++){cm[0][i]=0;}
  for(int i=1;i<6;i++){cm[i]=cm[i-1]+6;}
  
  cm[1][3]=cm[3][1]=1; cm[2][5]=cm[5][2]=1;   cm[2][4]=cm[4][2]=1;
  cm[3][5]=cm[5][3]=1; cm[4][5]=cm[5][4]=1;

  for(int i=1;i<6;i++){
    for(int j=1;j<6;j++){cout<<cm[i][j]<<" ";}cout<<endl;
  }

  burn_forest_clusters(cm, d, c, cs, &nc);
  return 0;
}
*/
