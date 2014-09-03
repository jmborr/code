#include "pseudoAA.h"
#include "pseudoPep.h"
#include "PDBLib.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc<4){
    cout << "usage: generatePPDB.linux inputMethod inputFile outPDB" << endl;
    cout << " inputMethod : 0/1/2/3" << endl;
    cout << "    0 --- input seq file only" << endl;
    cout << "    1 --- input PDB file but using SEQ only" << endl;
    cout << "       these two case, the generated peptide is streched" << endl;
    cout << "    2 --- input PDB file and fix BACKBONE and use approximate SideChain INFOR" << endl;
    cout << "       the sidechain geometry is satisfied by the DESIGN" << endl;
    cout << "    3 --- input PDB file and fix BACKBONE and SIDECHAIN from PDB" << endl;
    cout << "       the sidechain geometry might be away from DESIGN, need relax" << endl;
    exit(1);
  }
  
  pseudoPep p(argv[2], (input_t)atoi(argv[1]));
  p.printPDB(argv[3]);
}
