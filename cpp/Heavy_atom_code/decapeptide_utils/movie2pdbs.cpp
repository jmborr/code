#include "pdbClasses2.h"
#include "movie.h"
using namespace std;

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( (number_of_arguments != 5) && (number_of_arguments != 6) ){
    system("clear");
    cout << "Usage: ./movie2pdbs.x movie(bin) model.pdb startFrame nFrame movie.pdb [nskip]\n" ;
    cout << "movie(bin): output movie from a DMD simulation\n" ;
    cout << "model.pdb : any configuration of the system in PDB format\n" ;
    cout << "startFrame: first frame to translate to PDB format\n" ;
    cout << "nFrame    : number of frames to translate\n" ;
    cout << "movie.pdb : output file to store the frames in PDB format\n" ;
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
  ifstream MODEL( argv[ 2 ] ) ;  
  PDBchain chain( MODEL ) ; 
  int startFrame = atoi( argv[ 3 ] ) ; /*cout<<"startFrame="<<startFrame<<endl;*/
  int nFrame = atoi( argv[ 4 ] ) ; /*cout<<"nFrame="<<nFrame<<endl;*/
  ofstream PDBMOV( argv[ 5 ] ) ;
  int nskip=0 ;
  int narg = argc -1 ;
  if( narg == 6 ){ nskip = atoi( argv[ 6 ] ) ; }
  int natoms = chain.numberOfAtoms( ) ; /*cout<<"natoms="<<natoms<<endl;*/
  chain.renumberFully( 1 ) ;
  chain.renumberFullyAtomSerialNumber( 1 ) ;

  double** r = m.getCoords( ) ;
  double boundary = m.getDimensions( )[ 0 ] ;

  int n=0 ;
  int natoms_m = m.getNAtoms( ) ;/*atoms per movie frame*/
  /*cout<<"natoms_m="<<natoms_m<<endl;*/
  if( ( natoms_m % natoms ) ){
    cout << "mismach of input pdb and movie" << endl ; 
    return 1 ;
  }
  /*int total_number_frames = m.getNFrames() ;
  if( total_number_frames < startFrame + nFrame*(nskip+1) ){
    cout << "Movie "<<argv[ 1 ]<<" has only "<<total_number_frames<<" frames!\n" ; 
    return 1 ;
    }*/
  int nProts = natoms_m / natoms; /*if there is more than one
				    peptide in the system*/
  /*cout<<"nProts="<<nProts<<endl;*/
  m.jumpto( startFrame ) ;
  m.nextFrame( ) ;
  m.prevFrame( ) ;
 
  for( int iFrame = 0 ; iFrame<nFrame ; iFrame++ ){
    /*cout<<"iFrame="<<iFrame<<endl;*/
    double cm[ nProts ][ 3 ] ;/*center of mass of each peptide*/
    for( int ip=0; ip<nProts; ip++ ){/*calculate the cm's*/
      cm[ ip ][ 0 ] = r[ ip * natoms ][ 0 ] ;
      cm[ ip ][ 1 ] = r[ ip * natoms ][ 1 ] ;
      cm[ ip ][ 2 ] = r[ ip * natoms ][ 2 ] ;
      for( int i=1+ip*natoms; i<natoms*(ip+1); i++ ){
	for( int k=0; k<3; k++ ){
	  if( r[i][k] - r[i-1][k] >  boundary/2.0f ){ r[i][k] -= boundary ; }
	  if( r[i][k] - r[i-1][k] < -boundary/2.0f ){ r[i][k] += boundary ; }
	  cm[ip][k] += r[i][k] ;
	}
      } 
      for( int i=0; i<3; i++ ){ cm[ip][i] /= (double)natoms ; }
    }

    int alloc[ nProts ] ;
    for( int ip=0; ip<nProts; ip++ ){ alloc[ ip ] = 1 ; }
    int next = 0 ;
    alloc[ next ] = 0 ;
    double cmCurr[ 3 ] ;

    for( int ncount=1; ncount<nProts; ncount++ ){
      for( int j=0; j<3; j++ ){ cmCurr[ j ] = 0 ; }
      int ialloc = 0 ;
      for( int i=0; i<nProts; i++ ){
	if( !alloc[i] ){
	  ialloc++;
	  for( int j=0; j<3; j++ ){ cmCurr[j] += cm[i][j] ; }
	}
      }
      for( int j=0; j<3; j++){ cmCurr[j] /= (double)ialloc ; }

      double dist = 0, min = 3000 * pow( boundary, 2 ) ;
      int imin ;
      for( int i=0; i<nProts; i++ ){
	if(alloc[i]){
	  for( int j=0; j<3; j++ ){
	    double temp = fabs( cm[i][j] - cmCurr[j] ) ;
	    if( temp>boundary/2.0f ){ temp -= boundary ; }
	    dist += temp * temp ;
	  }
	  if( min > dist ){
	    min = dist ;
	    imin = i ;
	  }
	}
      }

      for( int k=0; k<3; k++ ){
	if( cm[imin][k] - cmCurr[k] >  boundary / 2.0f ){
	  cm[imin][k] -= boundary ;
	  for( int i=imin*natoms ; i<(imin+1)*natoms; i++ ){ r[i][k]-=boundary; }
	}
	if( cm[imin][k] - cmCurr[k] < -boundary / 2.0f ){
	  cm[imin][k] += boundary ;
	  for( int i= imin*natoms; i<(imin+1)*natoms; i++ ){ r[i][k]+=boundary; }
	}
      }
      alloc[imin]=0;
    }/*matches for( int ncount=1; ncount<nProts; ncount++ )*/

    /*calculate center of mass of whole system*/
    for( int i=0; i<3; i++ ){ cmCurr[i] = cm[0][i] ; }
    for( int ip=1; ip<nProts; ip++ ){
      for( int i=0; i<3; i++){ cmCurr[i] += cm[ip][i] ; }
    }
    for( int i=0; i<3; i++){ cmCurr[i] /= (double)nProts ; }
    for( int i=0; i<natoms_m; i++ ){ 
      for( int k=0; k<3; k++){
	r[i][k] -= ( cmCurr[k] - boundary/2.0 ) ;
      }
    }

    int index=0;
    /*double *junk = new double[3] ;*/
    for(int iProts=0; iProts<nProts; iProts++){
      for( int i=1; i<=natoms; i++ ){
	/*junk[0]=r[index][0]/2.3f;junk[1]=r[index][1]/2.3f;junk[2]=r[index][2]/2.3f;
	  chain.changeCoordOfAtomWithIndex( i, junk ) ;*/
	chain.changeCoordOfAtomWithIndex( i, r[index] ) ;
	index++ ;
      }
      PDBMOV << chain ;
      if( iProts<nProts-1 ){ PDBMOV <<"TER\n"; }
    }
    PDBMOV << "END" << endl ;
    m.nextFrame( ) ;
    for( int i=0; i<nskip; i++ ){ m.nextFrame( ) ; }
  }/*matches for(int iFrame=0; iFrame<nFrame; iFrame++)*/

  MODEL.close( ) ;
  return 0 ;
}
