#include "movie.h"
#include "go.h"
using namespace std;

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( (number_of_arguments !=10) && (number_of_arguments != 11) ){
    system("clear");
    cout << "Usage: ./fraction_nat.x movie(bin) ts nat.pdb frc_fld f_go min_sep subset startFrame nFrame outFile [nskip]\n" ;
    cout << "movie(bin): output movie from a DMD simulation\n" ;
    cout << "ts        : time step between consecutive frames\n" ;
    cout << "nat.pdb   : native configuration of the system in PDB format\n" ;
    cout << "frc_fld   : force field file with info on atoms (mass, radius, bonded\n" ;
    cout << "f_go      : (1+go_co_f)*(VdW_i+VdW_j) will be Go range for i,j\n" ;
    cout << "min_sep   : atoms belonging to amino acids that differ in their residue\n" ;
    cout << "numbers by min_sep or less, will not have native contacts\n" ;
    cout << "subset    : file with sets of atoms for fraction of native contacts\n" ;
    cout << "startFrame: first frame to translate to PDB format\n" ;
    cout << "nFrame    : number of frames to translate\n" ;
    cout << "outFile   : output file to store the fraction values\n" ;
    cout << "[nskip]   : skip nskip frames in between two consecutive\n" ;
    cout << "            translated frames (default is 0 )\n" ;
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
  double ts = atof( argv[ 2 ] ) ; /*cout<<"ts="<<ts<<endl;*/
  ifstream NATIVE( argv[ 3 ] ) ;  
  ifstream DAT( argv[ 4 ] ) ;
  static const double f_go = atof( argv[ 5 ] ) ; /*cout<<"f_go="<<f_go<<endl;*/
  static const int min_sep=atoi(argv[ 6 ]);/*cout<<"min_sep="<<min_sep<<endl;*/
  ifstream SUBSET( argv[ 7 ] ) ;
  int startFrame = atoi(argv[ 8 ]) ; /*cout<<"startFrame="<<startFrame<<endl;*/
  int nFrame = atoi( argv[ 9 ] ) ;  /*cout<<"nFrame="<<nFrame<<endl;*/
  ofstream OUT( argv[ 10 ] ) ;
  int nskip=0 ;
  int narg = argc -1 ;
  if( narg == 11 ){ nskip = atoi( argv[ 11 ] ) ; }

  /*renumber chain*/
  PDBchain native( NATIVE ) ; NATIVE.close( ) ;
  int natoms = native.numberOfAtoms( ) ;
  native.renumberFully( 1 ) ;
  native.renumberFullyAtomSerialNumber( 1 ) ; /*cout<<native;*/

  /*store subsets into list of PDBatoms*/
  listPDBatom ss1[100], ss2[100] ;
  int n_ss = 0 ;
  string line ;
  do{ 
    SUBSET >> ss1[ n_ss ] ; 
    getline( SUBSET, line ) ;/*cout<<"line="<<line<<endl;*/
    if(line.find("END")!=string::npos){ ss2[n_ss] = ss1[n_ss] ; }
    else if(line.find("TER")!=string::npos){ SUBSET >> ss2[ n_ss ] ; }
    else{ cout<<"\nERROR:line neither END nor TER !\n" ; exit(1) ; }
    n_ss++ ;
  }while( !SUBSET.eof( ) ) ;
  SUBSET.close( ) ; n_ss-- ;  /*cout<<"n_ss="<<n_ss<<endl;*/
  /*for(int i=0;i<n_ss;i++){cout<<ss1[i]<<"TER\n"<<ss2[i]<<"END\n";}*/

  /*allocate and fill lists of native contacts for each subset*/
  int *n_nc_ss = new int[ n_ss ] ; /*number of native contacts per subset*/
  int **p_ss = new int*[ n_ss ],  **q_ss = new int*[ n_ss ] ;/*nat cont*/
  double **go_ranges = new double*[ n_ss ] ;/*range for each nat cont*/
  for( int i=0; i<n_ss; i++ ){
    p_ss[i]=NULL; q_ss[i]=NULL; go_ranges[i]=NULL, n_nc_ss[i]=0 ;
    output_nat_cont( native, ss1[i], ss2[i], min_sep, f_go, DAT,
		     p_ss[i], q_ss[i], go_ranges[i], n_nc_ss[i] ) ;
    /*for( int j=0; j<n_nc_ss[i]; j++ ){
      cout<<"( "<<p_ss[i][j]<<", "<<q_ss[i][j]<<" ), "<<go_ranges[i][j]<<endl;
      }*/
  }

  /*output header*/
  char cbuf[256], cbuf2[32] ;
  sprintf( cbuf, "#    time " ) ;
  for( int i=0; i<n_ss; i++ ){
    sprintf( cbuf2, " s%d-%2d", i+1,n_nc_ss[i] ) ;
    strcat( cbuf, cbuf2 ) ;
  }
  OUT<<cbuf<<endl;

  double f ; /*fraction of native contacts*/
  double t= startFrame*ts ; /*frame time*/
  double** r = m.getCoords( ) ;
  double boundary = m.getDimensions( )[ 0 ] ;
  int natoms_m = m.getNAtoms( ) ;/*atoms per movie frame*/
  if( ( natoms_m % natoms ) ){
    cout << "mismach of input pdb and movie" << endl ; return 1 ;
  }

  m.jumpto( startFrame ) ;
  m.nextFrame( ) ;
  m.prevFrame( ) ;
  
  for( int iFrame=0 ; iFrame<nFrame ; iFrame++ ){
    shift_to_center( r, natoms, boundary ) ;
    sprintf( cbuf, "%10.1f  ", t ) ; 
    /*if(iFrame==nFrame-1){
      PDBchain junkChain; 
      junkChain=native;
      for(int i=1;i<=natoms;i++){
	junkChain.changeCoordOfAtomWithIndex(i,r[i-1]);
      }
      cout<<junkChain<<"END\n";
      for(int i=1;i<=natoms;i++){
	cout<<i<<" "<<r[i-1][0]<<" "<<r[i-1][1]<<" "<<r[i-1][2]<<endl;
      }
    }*/
    for( int i=0; i<n_ss; i++ ){
      f = output_fraction( native, r, p_ss[i], q_ss[i], go_ranges[i], n_nc_ss[i] ) ;
      sprintf( cbuf2, "%4.2f  ", f ) ; strcat( cbuf, cbuf2 ) ;
    }
    OUT<<cbuf<<endl;
    m.nextFrame( ) ; t += ts ;
    for( int i=0; i<nskip; i++ ){ m.nextFrame( ) ; t += ts ; }
  }
  return 0 ;
}
