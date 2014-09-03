#include "second.h"
#include "PDBLib.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc<3){
    cout << "usage: dmd9_pdbs2fy.x input_PDB out_fy" << endl;
    exit(1);
  }
  ifstream in(argv[1]);
  string line;
  char command[300];
  while(getline(in, line)){
    ofstream out("junk_pdbs2fy");
    while(line.find("END") == string::npos){
      out<<line<<endl;
      getline(in, line);
    }
    out<<"END\n";
    out.close();
    sprintf(command, "dmd9_pdb2fy.x junk_pdbs2fy junk2_pdbs2fy");
    system(command);
    sprintf(command, "cat junk2_pdbs2fy >> %s",argv[2]);
    system(command);
  }
  sprintf(command, "rm -f junk_pdbs2fy junk2_pdbs2fy");
  system(command);
  return 0;
}
