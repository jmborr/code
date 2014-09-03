#include "movie.h"
using namespace std;
/*-------------------------------------------
  This program read the movie file using movie class and 
  gives the simple informations in the movie file.

  Author: Feng Ding
  fding@bu.edu
 -------*/

int main(int argc, char* argv[]){
  if(argc<2){
    cout << "usage: readMovHead.linux movie(bin)" << endl;
    cout << "---- To get the infor of the movie files" << endl;
    cout << "[OUTPUT]: # of Frames, Detla T, Box Dimension, # of Atoms" << endl;
    exit(1);
  }
  
  movie m(argv[1]);
  if(m.getErrorStatus()){
    cout << m.getMsg() << endl;
    exit(1);
  }
  
  char buf[1024];
  sprintf(buf, "%ld %lf %lf %ld", m.getNFrames(), m.getDeltaT(),
	  m.getDimensions()[0], m.getNAtoms());
  cout << buf << endl;
}
