#include "pdbClasses2.h"

/*=====================================================================*/
int storeSingleChain(char *pdbf, string chainID, char *outf, string endline){
  ofstream out(outf);
  char *command=new char[512];
  PDBchain chain;
  if(!chain.importChain(pdbf,chainID)){
    cerr<<"ERROR from store_single_chain: unsuccessful importChain for"
	<<pdbf<<"\n";
    return 0;
  }
  chain.simplify();
  /*1 means do not output atom info beyond coordinates field*/
  if(chain.exportChain(outf,1,endline)==0){
    cerr<<"ERROR from store_single_chain: unsuccessful exportChain for"
	<<pdbf<<"\n";
    return 0;
  }
  /*missing residues and backbone*/
  if( chain.numberOfAtoms()/chain.length() < 4){
    sprintf(command,"./peter/pulchra -vpc %s && echo %s.rebuilt && /bin/mv %s.rebuilt %s",outf,outf,outf,outf);
    if(system(command)){
      cerr<<"ERROR from store_single_chain: "
	  <<"unsuccessful pulchra refinement for"<<pdbf<<"\n";
      return 0;
    }
  }
  return 1;
}
/*=====================================================================*/
int extractSegment(char *pdbf, string chainID, char *outf,
		    int startResSeq, int endResSeq,
		    string startiCode, string endiCode, string endline){
  ofstream out(outf);
  char *command=new char[512];
  PDBchain chain,segment;
  if(!chain.importChain(pdbf,chainID)){
    cerr<<"ERROR from extract_segment: unsuccessful importChain for"<<pdbf<<"\n";
    return 0;
  }
  if(!chain.extract_segment(segment,startResSeq,endResSeq,startiCode,endiCode)){
    cerr<<"ERROR from extract_segment: unsuccessful extract_segment for"<<pdbf<<"\n";
    return 0;
  }
  segment.simplify();
  /*1 means do not output atom info beyond coordinates field*/
  if(segment.exportChain(outf,1,endline)==0){
    cerr<<"ERROR from extract_segment: unsuccessful exportChain for"<<pdbf<<"\n";
    return 0;
  }
  /*missing residues and backbone*/
  if( segment.numberOfAtoms()/segment.length() < 4){
    sprintf(command,"./peter/pulchra -vpc %s && /bin/mv %s.rebuilt %s",
	    outf,outf,outf);
    if(system(command)){
      cerr<<"ERROR from extract_segment: "
	  <<"unsuccessful pulchra refinement for"<<pdbf<<"\n";
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
/*=====================================================================*/
int check_CAs(char *pdbf){
  PDBchain chain(pdbf);
  int rms=0,i=1,l=chain.length();
  /*printf("l=%d\n",l); cout<<chain;*/
  while(i<=l){
    /*printf("%d ",i);*/
    if( !chain.isAtomNameAtAminoAcidIndex(" CA ",i) ){
      rms+=chain.removeResAtIndex(i);
      l--; /*the chain has now one less residue*/
      i--; /*a new residue is occupying position "i" now*/
    }
    i++;
  }
  /*If we removed residues, then rearrange resSeq and atom number && rewrite pdb file*/
  if(rms){
    chain.renumberFully(1);
    chain.renumberFullyAtomSerialNumber(1);
    ofstream out(pdbf);  out<<chain;  out.close();
  }
  return rms;
}
/*======================================================================
int main(){
  char *pdbf=new char[256],*outf=new char[256],*outf2=new char[256];
  //  sprintf(pdbf,"./pdb12e8.ent");
  sprintf(pdbf,"./pdb1amt.ent");
  sprintf(outf,"/gpfs1/active/jose/code/projects/updateMkLib/junk");
  sprintf(outf2,"/gpfs1/active/jose/code/projects/updateMkLib/junk2");
  string id("A");
  int start=81,end=102;
  storeSingleChain(pdbf,id,outf,"TER\n");
  extractSegment(pdbf,id,outf2,start,end," "," ","TER\n");
  cout<<OneLetterSeq(outf)<<endl;
  cout<<OneLetterSeq(outf2)<<endl;
  return 0;
}
=====================================================================*/
