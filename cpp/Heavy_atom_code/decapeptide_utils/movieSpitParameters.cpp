#include "pdbClasses2.h"
#include "movie.h"
using namespace std;

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( (number_of_arguments != 4) && (number_of_arguments != 5) ){
    system("clear");
    cout << "Usage: ./movieSpitParameters.x movie(bin) startFrame nFrame outFile [nskip]\n" ;
    cout << "movie(bin): output movie from a DMD simulation\n" ;
    cout << "model.pdb : any configuration of the system in PDB format\n" ;
    cout << "startFrame: first frame to translate to PDB format\n" ;
    cout << "nFrame    : number of frames to translate\n" ;
    cout << "outFile : output file to store the parameter values (temp,energy,...)\n" ;
    cout << "[nskip]   : skip nskip frames in between two consecutive translated frames (default is 0 )\n" ;
    cout << "\n\n" ;
    return false ;
  }
  else return true ;
}
/*=====================================================*/


int main(int argc, char* argv[]){

  /*check arguments*/
  if( test_input(argc, argv) == false ) { return 1 ; }
  
  /*store arguments in memory*/
  movie m( argv[ 1 ] ) ;
  if( m.getErrorStatus( ) ){ cout << m.getMsg( ) << endl ; }
  int startFrame = atoi(argv[2]); /*cout<<"startFrame="<<startFrame<<endl;*/
  int nFrame = atoi( argv[ 3 ] ) ; /*cout<<"nFrame="<<nFrame<<endl;*/
  ofstream OUT( argv[ 4 ] ) ;
  int nskip=0 ;
  int narg = argc -1 ;
  if( narg == 5 ){ nskip = atoi( argv[ 5 ] ) ; }

  int nParameters=m.getNParameters();
  double *values= m.getParameterValues();
  char **names=m.getParameterNames();
  m.jumpto( startFrame ) ;
  m.nextFrame( ) ;
  m.prevFrame( ) ;
  char buff[1024],buff2[64];;

  for(int i=0;i<nParameters;i++){ OUT<<names[i]<<"  "; }  OUT<<endl;
  for( int iFrame = 0 ; iFrame<nFrame ; iFrame++ ){
    /*cout<<"iFrame="<<iFrame<<endl;*/
    m.nextFrame( ) ;

    sprintf(buff,"");
    for(int i=0;i<nParameters;i++){ 
      sprintf(buff2,"%12.7f ",values[i]);
      strncat(buff,buff2,strlen(buff2));
    }
    OUT<<buff<<endl;

    for( int i=0; i<nskip; i++ ){ m.nextFrame( ) ; }
  }/*matches for(int iFrame=0; iFrame<nFrame; iFrame++)*/

  OUT.close( ) ;
  return 0 ;
}
