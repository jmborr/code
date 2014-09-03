#include "pdbClasses2.h"
#include "movie.h"
using namespace std;

/*=================================================================*/
void shift_to_center( double **r, const int &natoms, const double &boundary ){
  double bH = boundary/2.0 ;
  double dr[3], old_r[3]={r[0][0], r[0][1], r[0][2] } ;
  for( int i=1; i<natoms; i++ ){/*remove periodic boundary conditions*/
    for( int j=0; j<3; j++ ){
      dr[j] = r[i][j]-old_r[j] ;
      if( dr[j] > bH ){ dr[j] -= boundary ; }
      else if( dr[j] < -bH ){ dr[j] += boundary ; }
      old_r[j] = r[i][j] ;
      r[i][j] = r[i][j] + dr[j] ;
    }
  }
  double gc[3]={ 0.0, 0.0, 0.0 } ;/*geometric center*/
  for(int i=0;i<natoms;i++){for(int j=0;j<3;j++){gc[j]+=r[i][j];}}
  for(int j=0;j<3;j++){ gc[j] /= natoms ; }
  for(int i=0;i<natoms;i++){for(int j=0;j<3;j++){r[i][j]+=bH-gc[j];}}
}
/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( (number_of_arguments !=7) && (number_of_arguments != 8) ){
    system("clear");
    cout << "Usage: ./rmsd_subset.x movie(bin) ts nat.pdb subset startFrame nFrame outFile [nskip]\n" ;
    cout << "movie(bin): output movie from a DMD simulation\n" ;
    cout << "ts        : time step between consecutive frames\n" ;
    cout << "nat.pdb   : native configuration of the system in PDB format\n" ;
    cout << "subset    : file with sets of atoms to find rmsd\n" ;
    cout << "startFrame: first frame to translate to PDB format\n" ;
    cout << "nFrame    : number of frames to translate\n" ;
    cout << "outFile   : output file to store the rmsd values\n" ;
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
  double ts = atof( argv[ 2 ] ) ;
  ifstream NATIVE( argv[ 3 ] ) ;  PDBchain native( NATIVE ) ; NATIVE.close( ) ;
  int natoms = native.numberOfAtoms( ) ;
  native.renumberFully( 1 ) ;
  native.renumberFullyAtomSerialNumber( 1 ) ;
  int startFrame = atoi( argv[ 5 ] ) ; 
  int nFrame = atoi( argv[ 6 ] ) ; 
  ofstream OUT( argv[ 7 ] ) ;
  int nskip=0 ;
  int narg = argc -1 ;
  if( narg == 8 ){ nskip = atoi( argv[ 8 ] ) ; }

  /*store substets into list of PDBatoms*/
  ifstream SUBSET( argv[ 4 ] ) ;
  listPDBatom subsets[100] ;
  int n_subsets = 0 ;
  string line ;
  do{ SUBSET >> subsets[ n_subsets++ ] ; }while( !SUBSET.eof( ) ) ;
  SUBSET.close( ) ; n_subsets-- ;

  /*allocate coordinates of subsets to arrays*/
  int **indexes = new int *[ n_subsets ] ;
  double ***r_ss = new double **[ n_subsets ] ;
  double ***t_ss = new double **[ n_subsets ] ;
  /*double ***junk = new double **[ n_subsets ] ;*/
  int *l = new int[ n_subsets ] ;
  for( int i=0; i<n_subsets; i++ ){
    l[i] = subsets[i].length( ) ;
    indexes[i] = new int[ l[i] ] ;
    r_ss[ i ] = new double*[ l[i] ] ;
    t_ss[ i ] = new double*[ l[i] ] ;
    /*junk[ i ] = new double*[ l[i] ] ;*/
    r_ss[ i ][ 0 ] = new double[ l[i] * 3 ] ;
    t_ss[ i ][ 0 ] = new double[ l[i] * 3 ] ;
    /*junk[ i ][ 0 ] = new double[ l[i] * 3 ] ;*/
    for( int j=1; j<l[i]; j++ ){
      r_ss[ i ][ j ] = r_ss[ i ][ j-1 ] + 3 ;
      t_ss[ i ][ j ] = t_ss[ i ][ j-1 ] + 3 ;
      /*junk[ i ][ j ] = t_ss[ i ][ j-1 ] + 3 ;*/
    }
  }


  /*begin reading movie and check number of atoms*/
  double t=0.0, **r = m.getCoords( ) ;
  double boundary = m.getDimensions( )[ 0 ] ;
  int n=0 ;
  int natoms_m = m.getNAtoms( ) ;/*atoms per movie frame*/
  if( ( natoms_m % natoms ) ){
    cout << "mismach of input pdb and movie" << endl ; 
    return 1 ;
  }

  /*fill arrays r_ss and indexes. Put GM of each list in the center of box*/
  PDBatom *atom ;
  for( int i=0; i<n_subsets; i++ ){
    for( int j=0; j<l[i]; j++ ){
      atom = subsets[i].getAtomAt( j+1 ) ;
      indexes[i][j] = atom->getSerialNumber( ) ; 
      indexes[i][j]-- ;
      atom->getCoordToDoubleArrayII( r_ss[i][j] ) ;
      /*atom->getCoordToDoubleArrayII( junk[i][j] ) ;*/
      delete atom ;
    }
    shift_to_center( r_ss[i], l[i], boundary ) ;
    /*shift_to_center( junk[i], l[i], boundary ) ;*/
  }

  m.jumpto( startFrame ) ;
  m.nextFrame( ) ; t += ts ;
  m.prevFrame( ) ;
  double rmsd ;
  char cbuf[512], cbuf2[64] ;

  /*output header*/
  sprintf( cbuf, "#    time   " ) ;
  for( int i=0; i<n_subsets; i++ ){
    sprintf( cbuf2, " subset%d ", i+1 ) ;
    strcat( cbuf, cbuf2 ) ;
  }
  sprintf( cbuf2, "\n" ) ; strcat( cbuf, cbuf2 ) ;
  OUT<<cbuf;

  for( int iFrame = 0 ; iFrame<nFrame ; iFrame++ ){

    /*fill array t_ss with fresh coordinates*/
    for( int i=0; i<n_subsets; i++ ){
      for( int j=0; j<l[i]; j++ ){
	/*cout<<"natoms="<<natoms<<" indexes["<<i<<"]["<<j<<"]="<<indexes[i][j]<<endl;*/
	for( int k=0; k<3; k++ ){
	  t_ss[i][j][k] = r[ indexes[i][j] ][ k ] ;
	}
      }
      /*we nedd GM of r_ss[i] and t_ss[i] coincide for get_rms to work*/
      shift_to_center( t_ss[i], l[i], boundary ) ;
    }

    /*calculate the rmsd's and output*/
    sprintf( cbuf, "%10.1f  ", t ) ; 
    for( int i=0; i<n_subsets; i++ ){
      rmsd = get_rms( r_ss[i], t_ss[i], l[i] ) ;
      /*rmsd = get_rms( junk[i], t_ss[i], l[i] ) ;*/
      sprintf( cbuf2, "%7.2f  ", rmsd ) ;
      strcat( cbuf, cbuf2 ) ;
    }
    sprintf( cbuf2, "\n" ) ; strcat( cbuf, cbuf2 ) ;
    OUT<<cbuf;

    /*now junk becomes t_ss
    for( int i=0; i<n_subsets; i++ ){
      for( int j=0; j<l[i]; j++ ){
	for( int k=0; k<3; k++ ){ junk[i][j][k] = t_ss[i][j][k] ; }
      }
      }*/

    m.nextFrame( ) ;  t += ts ;
    for( int i=0; i<nskip; i++ ){ m.nextFrame( ) ;  t += ts ; }
  }/*matches for(int iFrame=0; iFrame<nFrame; iFrame++)*/

  OUT.close( ) ;
  return 0 ;
}
