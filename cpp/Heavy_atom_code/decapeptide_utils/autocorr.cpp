using namespace std;
#include<iostream>
#include<fstream>
#include<fstream>
#include<cstring>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./autocorr.x inFile outFile\n";
    cout << "inFile : one-column file with signal values\n" ;
    cout << "outFile: autocorrelation function of the signal\n";
    cout << "\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char* argv[]){

  /*check arguments*/
  if( test_input(argc, argv) == false ) { return 1 ; }

  /*retrieve input*/
  ifstream IN(argv[1]);
  ofstream OUT(argv[2]);
  /*count number of measurements and calculate average*/
  double v=0,av=0;
  string line;
  int N=-1;
  do{ N++; av+=v; }while(IN>>v); 
  av/=N; /*cout<<"N="<<N<<endl<<"av="<<av<<endl;*/
  IN.clear( ios::goodbit ) ;  IN.seekg( 0 ) ; /*rewind to origing of IN*/

  /*store file in memory*/
  double *x=new double[N];
  for(int i=0; i<N; i++){IN>>x[i]; x[i]-=av; /*cout <<x[i]<<endl;*/ }
  /*av=0;for(int i=0; i<N; i++){av+=x[i];} cout<<"av="<<av<<endl;*/
  IN.close();
  
  
  double c,c0=0.0;

  for(int n=0; n<N; n++){ c0+=x[n]*x[n]; } c0/=N;

  for(int m=0; m<N-1; m++){
    c=0;
    for(int n=0; n<N-m; n++){ c+=x[n+m]*x[n]; }
    c/=(c0*(N-m));
    OUT<<c<<endl;
  }

  OUT.close();
  return 0;
}
