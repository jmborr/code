#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

class movie{
 private:
  /*------------------private entries--------------------------*/
  FILE* mov_path;
  short version;
  short nFrames;               /*# of frames*/
  int i_nFrames;
  short nAtomTypes;            /*# of atom type*/
  short nAtoms;                /*# of atoms*/
  short nParameters;           /*# of parameter*/
  char  ** parameterNames;     /*parameter names*/
  double * parameterValues;    /*parameter values*/
  double ** r;                 /*buffer to read and store the corrd*/
  double box[3];              /*box size array*/
  double scaleFactors[3];      /*scale factor*/
  double dt;                   /*dt between consecutive frames*/
  short nDimensions;           /*# of dimension*/
  int bo;                      /*byte order*/
  short  * atomColors;          
  double * atomRadii;          /*list of atom radii*/
  short  * atomTypes;          /*list of atom type*/
  short** bonds;               /*bond lists*/
  vector<long> frameOffset;    /*frame offsets*/
  vector<char> new_type;
  vector<char> new_bond;
  short s_dummy; 
  char  c_dummy;
  int DEBUG;
  int l_pass;                 /*length of password*/
  int errStatus;
  string errMsg;
  int readed;
  int current;

  /*-----------------private methods--------------------------*/
  /* determine the byte order of storage of current machine
     if it is BIG-ENDIAN, return 1;otherwise, return 0.  */
  int byte_order(void){
    unsigned short i=1;
    unsigned char * s;
    s=(unsigned char*)&i;
    if(s[1]>s[0])return 1;
    return 0;
  }
  
  /* because the binary is stored in big endian byte order,
     we need to read the data in according to the byte order*/
  int fread_bo(void* inpt, int size, int nitem, FILE* path){
    int iread,i;
    static unsigned char buf[64];
    void* top = inpt;
    for(iread=0; iread<nitem; iread++){
      if(fread(buf, size, 1, path)<1) return iread;
      if(bo){/*big endian*/
	for(i=0; i<size; i++)
	  *((char*)top+i) = buf[i];
      }
      else{/*little endian*/
	for(i=0; i<size; i++)
	  *((char*)top+i) = buf[size-1-i];
      }
      top = (char*)top + size;
    }
    return nitem;
  }
  
  int read_bond(){
    for(int i=0; i<nAtoms; i++){
      if(fread_bo(&s_dummy, sizeof(s_dummy), 1, mov_path) < 1)
	return 0;
      if(s_dummy!=i){
	errMsg = "wrong format at bonding.";
	return 0;
      }
      int j=0;
      while(1){
	if(fread_bo(&s_dummy, sizeof(s_dummy), 1, mov_path) < 1) 
	  return 0;
	if(s_dummy!=i){
	  bonds[i][j] = s_dummy;
	  j++;
	}
	else {
	  bonds[i][j] = -1;
	  break;
	}
      }
    }
    return 1;
  }
  
  int read_header(void);
  
  void make_offset(void);
  
  int readChanges(int type_frame, int bond_frame);
  
  int read_frame(char& typeChanges, char& bondChanges);
  
 public:
  /*----------PUBLIC METHODS-------------*/
  movie(char* movie_name, int DEBUG_OPT=0);
  
  ~movie();
  
  int nextFrame(void);
  
  int prevFrame(void);
  
  int jumpto(int target);
  
  int getNFrames(){return i_nFrames;}
  
  int getNAtomTypes(){return nAtomTypes;}
  
  double* getRadii(){return atomRadii;}
  
  int getNAtoms(){return nAtoms;}
  
  short* getAtomTypes(){return atomTypes;}
  
  double** getCoords(){return r;}
  
  short** getBonds(){return bonds;}
  
  int getNParameters(){return nParameters;}

  char** getParameterNames(){return parameterNames;}

  double* getParameterValues(){return parameterValues;}

  int getNDimensions(){return nDimensions;}

  double* getDimensions(){return box;}
  
  double getDeltaT(){return dt;}
  
  int getErrorStatus(){return errStatus;}

  string& getMsg(){return errMsg;}

  int getCurrentFrame(){return current;}

  /*
  int getReaded(){
    return readed;
  }

  int getSize(){
    return frameOffset.size();
  }
  */
};
