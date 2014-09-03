#include "pdbClasses2.h"
using namespace std;
#include<fstream>
#include<string>

/*=================================================================*/
bool test_input(int argc, char ** argv ){
  int number_of_arguments = argc -1 ;
  if( number_of_arguments != 2 ){
    system("clear");
    cout << "Usage: ./dihedral_plot.x four_atom_file out\n";
    cout << "four_atom_filed : four atoms in PDB format\n" ;
    cout << "out : file were output is written. Output is two columns of distance and rotamer angle, in increase of one degree from 0 to 360\n";
    cout << "Four atoms define two planes, thus a rotamer/dihedral-angle\n";
    cout << "\n";
    return false ;
  }
  else return true ;
}
/*=====================================================*/
int main(int argc, char* argv[]){

  /*check arguments*/
  if( test_input(argc, argv) == false ) { return 1 ; }

  ifstream ATOMS(argv[1]);
  ofstream OUT(argv[2]);
  PDBatom ats[4];
  PDBvector v[4];
  for(int i=0;i<4;i++){
    ATOMS>>ats[i];
    v[i]=ats[i].getCoord();
  }
  PDBvector a,b,c;
  a=v[1]-v[0]; b=v[2]-v[1]; c=v[3]-v[2];
  double da,db,dc;
  da=sqrt(a.norm2());   db=sqrt(b.norm2());   dc=sqrt(c.norm2());
  /*cout<<da<<" "<<db<<" "<<dc<<endl;exit(1);*/
  double cos_a,sin_a,cos_c,sin_c;
  cos_a=(a*b)/(da*db); sin_a=sqrt(1-cos_a*cos_a);
  cos_c=(c*b)/(dc*db); sin_c=sqrt(1-cos_c*cos_c);
  /*cout<<cos_a<<" "<<sin_a<<" "<<cos_c<<" "<<sin_c<<endl;exit(1);*/
  double c1=da*sin_a, c2=dc*sin_c;
  double h=0,hi=PI/180.0,f1,f2,f3=da*cos_a+db+dc*cos_c,d;
  char buf[100];
  for(int chi=0;chi<=360;chi++){
    f1=c1+c2*cos(h);
    f2=c2*sin(h);   
    d=sqrt(f1*f1+f2*f2+f3*f3);
    h+=hi;
    sprintf(buf,"%3d %lf\n",chi,d);
    OUT<<buf;
  }
  return 0 ;
}
