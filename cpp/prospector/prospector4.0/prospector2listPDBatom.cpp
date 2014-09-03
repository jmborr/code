#include "prospector2listPDBatom.h"

/*======================================================================*/
bool read_prosp_line(ifstream &in, string &out){
  getline(in,out); /*cout<<out<<endl;exit(0);*/
  int l=out.length();
  for(int i=l;i<80;i++) out.append(" "); /*line only misses proper length*/
  /*cout<<out.length()<<endl;*/ /*cout<<out<<endl;*/
}
/*======================================================================*/
int read_prosp_file(const char *name, listsPDBatom &ltps){
  char buff[256];
  string id,entry;
  int lt;
  double z;
  ifstream prosp(name);
  PDBatom at;
  listPDBatom tpl;

  while(!prosp.eof()){
    tpl.empty(); lt=0;
    prosp>>id>>lt; /*cout<<lt<<endl;*/ /*exit(0);*/ /*read header line*/
    prosp.ignore(256,'\n'); /*read rest of line*/
    for(int i=0;i<lt;i++){ /*read the ATOM entries*/
      read_prosp_line(prosp,entry);
      entry>>at; /*cout<<at<<endl;*/
      tpl.insertAtBack(at);
    }
    getline(prosp,entry); /*read TER line*/
    if(lt) ltps.insertAtBack(tpl,id); /*insert template and its 5-letter id*/
  }
}
/*========================================================================*/
