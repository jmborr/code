#include "second.h"
#include "PDBLib.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc<3){
    cout << "usage: pdb2fy.linux input_PDB out_fy" << endl;
    exit(1);
  }
  protein prot(argv[1]);
  ofstream out(argv[2]);
  double phi, psi;
  char line[100];
  for(int i=1; i<prot.getLength()-1; i++){
    prot.getPhiPsi(i, phi, psi);
    sprintf(line, "%4d %7.4f %7.4f\n",i+1,phi,psi);
    out << line;
  }

  return 0;
}
