#include "pseudoPep.h"
#include "movie.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc < 6){
    cout << "usage: movie2rms.linux movie(bin) inputMethod inputFile origPDB rmsd_f" << endl;
    cout << " inputMethod : 0/1" << endl;
    cout << "    0 --- input seq file only" << endl;
    cout << "    1 --- input PDB file but using SEQ only" << endl;
    exit(1);
  }
  pseudoPep p(argv[3], (input_t)atoi(argv[2]));
  int natoms = p.makeIndex();
  
  
  movie m(argv[1]);
  if(m.getErrorStatus()){
    cout << m.getMsg() << endl;
  }
  double** r = m.getCoords();
  double boundary = m.getDimensions()[0];
  
  int natoms_m = m.getNAtoms();
  if((natoms_m%natoms)){
    cout << "mismach of input pdb and movie" << endl;
  }
  int nProts = natoms_m / natoms;

  char tmp_fname[100];
  sprintf(tmp_fname, ".tmp_pdb_%s", argv[5]);
  cout << tmp_fname << " " << m.getNFrames() << endl;
  char buf[1024];
  ofstream tmp_pdb(tmp_fname);
  if(!tmp_pdb){
    cout << "can not open " << tmp_fname  << endl;
    exit(1);
  }
  for(int iFrame=0; iFrame<m.getNFrames(); iFrame++){
    tmp_pdb.seekp(0, ios_base::beg);
    m.nextFrame();
    {
      int& np = nProts;
      int& nres = natoms;
      int& npart = natoms_m;

      double cm[np][3];
      for(int ip=0; ip<np; ip++){
	cm[ip][0]=r[ip*nres][0];
	cm[ip][1]=r[ip*nres][1];
	cm[ip][2]=r[ip*nres][2];
	for(int i=1+ip*nres; i<nres*(ip+1); i++){
	  if(r[i][0]-r[i-1][0] > boundary/2.0f ) r[i][0]-=boundary;
	  if(r[i][0]-r[i-1][0] < -boundary/2.0f) r[i][0]+=boundary;
	  if(r[i][1]-r[i-1][1] > boundary/2.0f ) r[i][1]-=boundary;
	  if(r[i][1]-r[i-1][1] < -boundary/2.0f) r[i][1]+=boundary;
	  if(r[i][2]-r[i-1][2] > boundary/2.0f ) r[i][2]-=boundary;
	  if(r[i][2]-r[i-1][2] < -boundary/2.0f) r[i][2]+=boundary;
	  cm[ip][0]+=r[i][0];
	  cm[ip][1]+=r[i][1];
	  cm[ip][2]+=r[i][2];
	} 
	for(int i=0; i<3; i++)cm[ip][i]/=(double)nres;
      }
      
      int alloc[np];
      for(int ip=0; ip<np; ip++)alloc[ip]=1;
      int next=0;
      alloc[next]=0;
      double cmCurr[3];
      
      for(int ncount=1; ncount<np; ncount++){
	for(int j=0; j<3; j++) cmCurr[j]=0;
	int ialloc=0;
	for(int i=0; i<np; i++)if(!alloc[i]){
	  ialloc++;
	  for(int j=0;j<3;j++)cmCurr[j]+=cm[i][j];
	}
	for(int j=0; j<3; j++) cmCurr[j]/=(double)ialloc;
	
	double dist=0, min=3000*pow(boundary,2);
	int imin;
	for(int i=0; i<np; i++)if(alloc[i]){
	  for(int j=0; j<3;j++){
	    double temp = fabs(cm[i][j]-cmCurr[j]);
	    if(temp>boundary/2.0f)temp-=boundary;
	    dist+=temp*temp;
	  }
	  if(min>dist){
	    min = dist;
	    imin = i;
	  }
	}
	
	if(cm[imin][0]-cmCurr[0]>boundary/2.0f){
	  cm[imin][0]-=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][0]-=boundary;
	}
	if(cm[imin][0]-cmCurr[0]<-boundary/2.0f){
	  cm[imin][0]+=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][0]+=boundary;
	}
	
	if(cm[imin][1]-cmCurr[1]>boundary/2.0f){
	  cm[imin][1]-=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][1]-=boundary;
	}
	if(cm[imin][1]-cmCurr[1]<-boundary/2.0f){
	  cm[imin][1]+=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][1]+=boundary;
	}
	
	if(cm[imin][2]-cmCurr[2]>boundary/2.0f){
	  cm[imin][2]-=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][2]-=boundary;
	}
	if(cm[imin][2]-cmCurr[2]<-boundary/2.0f){
	  cm[imin][2]+=boundary;
	  for(int i= imin*nres; i<(imin+1)*nres; i++) r[i][2]+=boundary;
	}
	alloc[imin]=0;
      }
      
      for(int i=0; i<3; i++)cmCurr[i]=cm[0][i];
      for(int ip=1; ip<np;ip++)
	for(int i=0;i<3;i++)cmCurr[i]+=cm[ip][i];
      for(int i=0; i<3; i++)cmCurr[i]/=(double)np;
      for(int i=0; i<npart; i++){
	r[i][0]-=(cmCurr[0]-boundary/2.0);
	r[i][1]-=(cmCurr[1]-boundary/2.0);
	r[i][2]-=(cmCurr[2]-boundary/2.0);
      }
    }

    int index=0;
    for(int iProts=0; iProts<nProts; iProts++){
      for(int iAA=0; iAA<p.getLength(); iAA++){
	pseudoAA* theAA = p.getResidue(iAA);
	*theAA->getN()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	index++;
	*theAA->getC()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	index++;
	*theAA->getO()->getR() =  vec(r[index][0], r[index][1], r[index][2]);
	index++;
	*theAA->getCA()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	index++;
	if(theAA->getCB()){
	  *theAA->getCB()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	  index++;
	}
	if(theAA->getG()){
	  *theAA->getG()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	  index++;
	}
	if(theAA->getG2()){
	  *theAA->getG2()->getR() = vec(r[index][0], r[index][1], r[index][2]);
	  index++;
	}
      }
      p.constructHydrogen();
    }
    p.printPDB(tmp_pdb);
    sprintf(buf,"krms.linux %s %s >> %s", tmp_fname, argv[4], argv[5]);
    system(buf);
  }  
  tmp_pdb.close();
  sprintf(buf, "rm -f %s", tmp_fname);
  system(buf);
}
