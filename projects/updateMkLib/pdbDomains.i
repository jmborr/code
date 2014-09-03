/*compilation on my local computer:
 swig -c++ -python -o pdbDomains_wrap.cpp pdbDomains.i &&  g++ -Wno-deprecated -fPIC -c pdbDomains.cpp pdbClasses2.cpp && g++ -fPIC -c pdbDomains_wrap.cpp -o pdbDomains_wrap.o -I/usr/include/python2.4 -I/usr/lib/python2.4/config/ && g++ -shared pdbDomains_wrap.o pdbDomains.o pdbClasses2.o -o _pdbDomains.so

  compilation on the cluster:
 swig -c++ -python -o pdbDomains_wrap.cpp pdbDomains.i &&  g++ -Wno-deprecated -fPIC -c pdbDomains.cpp pdbClasses2.cpp && g++ -fPIC -c pdbDomains_wrap.cpp -o pdbDomains_wrap.o -I/usr/include/python2.3 -I/usr/lib64/python2.3/config/ && g++ -shared pdbDomains_wrap.o pdbDomains.o pdbClasses2.o -o _pdbDomains.so
*/

%module pdbDomains
%include "std_string.i"
%{
using namespace std;
#include<string>
#include "pdbClasses2.h"
%}

extern int storeSingleChain(char *pdbf,std::string chaindID,char *outf,std::string endline);

extern int extractSegment(char *pdbf,std::string chaindID,char *outf,int startResSeqm, int endResSeq, std::string startiCode,std::string endiCode,std::string endline);

extern std::string OneLetterSeq(char *pdbf);

extern int check_CAs(char *pdbf);
