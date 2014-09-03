#include "second.h"
#include "PDBLib.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc<2){
    cout << "usage: pdb2seconds.linux PDB" << endl;
    exit(1);
  }
  protein prot(argv[1]);
  Dihedrals* fy=new Dihedrals[prot.getLength()];
  double phi, psi;
  for(int i=0; i<prot.getLength(); i++){
    prot.getPhiPsi(i, phi, psi);
    fy[i].init(phi, psi);
  }

  
  unsigned char* ss = new unsigned char[prot.getLength()];
  for(int i=0; i<prot.getLength(); i++)
    ss[i]=getSeconds(*(fy+i), *(fy+i-1), *(fy+i+1));
  /*check HELIX*/
  int pt=0;
  int start=0;
  while(pt<prot.getLength()){
    if(ss[pt]==HELIX){
      if(!start){//begin
	start=1;
	int j=1;
	while(static_cast<seconds>(ss[pt-j])==TURN1 && pt>j){
	  ss[pt-j]=HELIX;
	  j++;
	}
      }
    }
    else{
      if(start){//end of helix
	start=0;
	int j=0;
	while(static_cast<seconds>(ss[pt+j])==TURN1 && pt+j<prot.getLength()){
	  ss[pt+j]=HELIX;
	  j++;
	}
	pt+=j;
      }
    }
    pt++;
  }
  /*check STRAND*/
  pt=0;
  start=0;
  while(pt<prot.getLength()){
    if(ss[pt]==STRAND){
      if(!start){//begin
	start=1;
	if(pt>0 && ss[pt-1]==COIL)ss[pt-1]=STRAND;
      }
    }
    else{
      if(start){//end of strand
	start=0;
	if(ss[pt]==COIL)ss[pt]=STRAND;
      }
    }
    pt++;
  }
  
  char* buf = new char[prot.getLength()];
  for(int i=0; i<prot.getLength(); i++)
    buf[i]=second2Char(static_cast<seconds>(ss[i]));
  
  char tmp[80];
  int  gpt;
  int end=0;
  for(int lines=0; lines< static_cast<int>(prot.getLength()/50)+1; lines++){
    pt=0;
    for(int i=0; i<10; i++){
      if(lines*50+i*5+1<=prot.getLength())
	pt+=sprintf(tmp+pt,"%-6ld", lines*50+i*5+1);
      else break;
    }
    tmp[60] = '\0';
    cout << tmp << endl;
    pt=0;
    gpt=lines*50;
    for(int i=0; i<10; i++){
      for(int j=0; j<5; j++){
	tmp[pt++]=buf[gpt++];
	if(gpt==prot.getLength()){
	  end=1;
	  break;
	}
      }
      if(end)break;
      tmp[pt++]=' ';
    }
    tmp[pt] = '\0';
    cout << tmp << endl;
  }
}
