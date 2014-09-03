#include "pdbClasses2.h"

/*=====================================================================*/
int storeSingleChain(char *pdbf, string chainID, char *outf){
  ofstream out(outf);
  char *command=new char[512];
  PDBchain chain;
  if(!chain.importChain(pdbf,chainID)){
    cerr<<"ERROR from store_single_chain: unsuccessful importChain\n";
    return 0;
  }
  chain.simplify();
  /*1 means do not output atom info beyond coordinates field*/
  if(chain.exportChain(outf,1,"TER\n")==0){
    cerr<<"ERROR from store_single_chain: unsuccessful exportChain\n";
    return 0;
  }
  if( chain.numberOfAtoms()/chain.length() < 4){/*missing residues and backbone*/
    sprintf(command,"./peter/pulchra -vpc %s && /bin/mv rebuilt_%s %s",outf,outf,outf);
    if(system(command)){
      cerr<<"ERROR from store_single_chain: unsuccessful pulchra refinement\n";
      return 0;
    }
  }
  return 1;
}
/*=====================================================================*/
int extractSegment(char *pdbf, string chainID, char *outf,
		    int startResSeq, int endResSeq,
		    string startiCode, string endiCode){
  ofstream out(outf);
  char *command=new char[512];
  PDBchain chain,segment;
  if(!chain.importChain(pdbf,chainID)){
    cerr<<"ERROR from extract_segment: unsuccessful importChain\n";
    return 0;
  }
  if(!chain.extract_segment(segment,startResSeq,endResSeq,startiCode,endiCode)){
    cerr<<"ERROR from extract_segment: unsuccessful extract_segment\n";
    return 0;
  }
  segment.simplify();
  /*1 means do not output atom info beyond coordinates field*/
  if(segment.exportChain(outf,1,"TER\n")==0){
    cerr<<"ERROR from extract_segment: unsuccessful exportChain\n";
    return 0;
  }
  if( segment.numberOfAtoms()/segment.length() < 4){/*missing residues and backbone*/
    sprintf(command,"./peter/pulchra -vpc %s && /bin/mv rebuilt_%s %s",outf,outf,outf);
    if(system(command)){
      cerr<<"ERROR from extract_segment: unsuccessful pulchra refinement\n";
      return 0;
    }
  }
  return 1;
}
/*=====================================================================*/
string OneLetterSeq(char *pdbf){
  PDBchain chain(pdbf);
  return chain.output_one_letter_sequence();
}
/*======================================================================*/
int main(){
  char *pdbf=new char[256],*outf=new char[256],*outf2=new char[256];
  /*  sprintf(pdbf,"./pdb12e8.ent");*/
  sprintf(pdbf,"./pdb2cvx.ent");
  sprintf(outf,"/gpfs1/active/jose/code/projects/updateMkLib/junk");
  sprintf(outf2,"/gpfs1/active/jose/code/projects/updateMkLib/junk2");
  string id("A");
  int start=81,end=102;
  storeSingleChain(pdbf,id,outf);
  extractSegment(pdbf,id,outf2,start,end," "," ");
  cout<<OneLetterSeq(outf)<<endl;
  cout<<OneLetterSeq(outf2)<<endl;
  return 0;
}
/*=====================================================================*/
