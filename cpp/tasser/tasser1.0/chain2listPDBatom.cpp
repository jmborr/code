#include "chain2listPDBatom.h"

/*======================================================================*/
bool read_chain_line(ifstream &in, string &out){
  char buff[81],name[]=" CA ";
  int resSeq=-1;
  double r[3]={_INFc2l_,_INFc2l_,_INFc2l_};

  in>>resSeq>>r[0]>>r[1]>>r[2];
  if(resSeq<0) return false;
  for(int i=0;i<3;i++) if(r[i]==_INFc2l_) return false;

  /*cout<<"resSeq="<<resSeq<<endl; 
    cout<<r[0]<<endl<<r[1]<<endl<<r[2]<<endl; exit(0);*/

  sprintf(buff,
	  "ATOM  %5d %s      %4d    %8.3f%8.3f%8.3f                          ",
	  resSeq,name,resSeq,r[0],r[1],r[2]);
  out.assign(buff);

  return true;
}
/*======================================================================*/
int read_chain_file(const char *name, listsPDBatom &ltps){
  ifstream chain(name);
  int nt,lt;
  string out;
  PDBatom at; /*atom entry*/
  listPDBatom tpl; /*template entry*/

  chain>>nt>>out; /*number of templates*/
  for(int nti=0;nti<nt;nti++){
    chain>>lt; /*length of the current template*/
    tpl.empty();
    for(int ilt=0;ilt<lt;ilt++)
      if(read_chain_line(chain,out)){
	out>>at;  /*cout<<at<<endl;*/
	tpl.insertAtBack(at);
      }
    ltps.insertAtBack(tpl);
  }
}
/*======================================================================*/
