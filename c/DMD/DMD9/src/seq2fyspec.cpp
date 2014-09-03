#include "pseudoPep.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc<3){
    cout << "usage: sep2fyspec seq_f fy_spec_f" << endl;
    exit(1);
  }
  pseudoPep p(argv[1], static_cast<input_t>(0));
  p.makeIndex();
  ofstream out(argv[2]);
  out << p.getLength() << endl;
  for(int i=0; i<p.getLength(); i++){
    out << p.getResidue(i)->getN()->getIndex()  << " "
	<< p.getResidue(i)->getCA()->getIndex() << " " 
	<< p.getResidue(i)->getC()->getIndex() << endl;
  }
  out.close();
}
