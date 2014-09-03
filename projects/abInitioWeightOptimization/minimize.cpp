#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "Minuit/FCNBase.h"
#include "Minuit/FunctionMinimum.h"
#include "Minuit/MnPrint.h"
#include "Minuit/MnMigrad.h"

#include "stringUtils.h"
using namespace std;

#define nw 20 /*number of weights*/

/*
       How to compile in the cluster
g++ -static -c minimize.cpp -I/gpfs1/active/jose/code/cpp/Minuit && g++ -static -o minimize.x minimize.o /gpfs1/active/jose/code/cpp/Minuit/src/*.o

       How to run
./minimize.x -a /gpfs1/scratch/jose/abInitioWeightOptimization/100to199/all_energy_averages.dat

       Energy terms and corresponding weights in TASSER
      EHB
     $     =eh1*EHB1            !+1/r of Ca-SC
     $     +eh1a*EHB1a          !+1/r for non-parallel of Ca-Ca
     $     +eh1b*EHB1b          !excluded volumn of SC-SC
     $     +eh1c*EHB1c          !pair-wise potential of SC-SC
     $     +eh2*EHB2            !quarsi3 for SC-SC
     $     +eh3*EHB3            !enhance good piece
     $     +eh4*EHB4            !-1/r for parallel contact of Ca-Ca
     $     +eh5a*EHB5a          !H-bond energy (alpha)
     $     +eh5b*EHB5b          !H-bond energy (beta)
      ESHORT=
     $     +es2*ESHORT2       !bury potential for SG
     $     +er1*ESHORT3       !distance restrain from both threading (for C_a)
     $     +er3*ESHORT4       !contact restrain from threading (for SG)
     $     +er4*ESHORT4a      !deviation of contact restrain
     $     +es3*ESHORT5       !bias2,3: v(i)-v(i+4) anti/parallel; c(i)-c(i+2) anit/paralel
     $     +es3a*ESHORT5a     !crumpling 
     $     +es3b*ESHORT5b     !bias4 to predicted alpha/beta structure.
     $     +es3c*ESHORT5c     !bias1 to possible alpha/beta structure. 
     $     +es4*ESHORT6       !correlation of E13 of Ca
     $     +es5*ESHORT7       !correlation of E14, from both common and 2th specific data
     $     +es6*ESHORT8       !correlation of E15, from both common and 2th specific data
     $     +er5*ESHORT9       !contact restraints of CA
     $     +er6*ESHORT10      !Long-range distance restraints of CA
     $     +er7*ESHORT11      !RMSD deviation
      energy_tot=energy_tot
     $     +en1*eprofo        !environment potential
     $     +en2*E_cord        !deviation from predicted contact order
     $     +en3*E_cnum        !deviation from predicted contact number

*/

typedef vector<double> vecd; /*a shorthand*/

/*======================================================================*/
/*attributes and properties related to each query, specially averages over the decoys*/
class query{

private:
  string header;   /*pdb id of the query*/
  double R;        /*average of R=1-TMscore over all decoys, ie, <R>_d*/
  double R2;       /*<R*R>_d*/
  double e[nw];    /* average of each energy term over all decoys, <E_i>_d*/
  double en[nw];   /* native energy terms En_i*/
  double ee[nw][nw];/* average of correlations between energy terms, <E_i*E_j>_d*/
  double Re[nw];   /* <R*E_i>_d*/

public:

  /*read all atributes from a line pointer. Line format assumes that
    all info for a query resideds in the line in the following order:
    header R R2 Re[nw] en[nw] ee[nw,nw] e[nw]
  */
  query(const string &line){
    int n=0;
    vector<string> tokens;
    splitString(line," ",tokens,false);
    header=tokens[n++];
    R=atof(tokens[n++].c_str());
    R2=atof(tokens[n++].c_str());
    for(int i=0;i<nw;i++) Re[i]=atof(tokens[n++].c_str());
    for(int i=0;i<nw;i++) e[i]=atof(tokens[n++].c_str());
    for(int i=0;i<nw;i++) for(int j=0;j<nw;j++) ee[i][j]=atof(tokens[n++].c_str());
    for(int i=0;i<nw;i++) en[i]=atof(tokens[n++].c_str());
  }

  /*default destructor*/
  ~query(){};

  string getHeader() const { return header; }

  void printQuery() const{
    cout<<"##amt       EHB1       EHB1a       EHB1b        EHB1c       EHB2        EHB3         EHB4       EHB5a        EHB5b    ESHORT2     ESHORT5      ESHORT5a    ESHORT5b    ESHORT5c   ESHORT6      ESHORT7    ESHORT8      eprofo      E_cord      E_cnum\n";

    cout<<"#<mt> <mt^2>\n";
    printf("%8.6lf %8.6lf",R,R2);

    cout<<"\n#<mt*Ee> (1..M)\n";
    for(int i=0;i<nw;i++) printf(" %12.6lf",Re[i]);

    cout<<"\n#<Ee> (1..M)\n";
    for(int i=0;i<nw;i++) printf(" %12.6lf",e[i]);
	
    cout<<"\n#<Ee*Ee'> (1,1)..(1,M),..(M,M)\n";
    for(int i=0;i<nw;i++)
      for(int j=0;j<nw;j++) printf(" %15.6lf",ee[i][j]);

    cout<<"\n#Ene (1..M)\n";
    for(int i=0;i<nw;i++) printf(" %12.6lf",en[i]);

    cout<<endl;
  }

/*average of energy over all decoys < Sum_i(w_i*E_i) >_d == Sum_i(w_i*<E_i>_d)*/
  double E(const vecd &w) const {
    double z=0;
    for(int i=0;i<nw;i++) z+=w[i]*e[i];
    return z;
  }

/*energy of native state Sum_i(w_i*En_i) */
  double En(const vecd &w) const {
    double z=0;
    for(int i=0;i<nw;i++) z+=w[i]*en[i];
    return z;
  }

  /*average of the energy squared of all decoys < Sum_i(w_i*E_i) * Sum_i(w_i*E_i) >_d */
  double E2(const vecd &w) const {
    double x,y,z=0;
    for(int i=0;i<nw;i++) z+=w[i]*w[i]*ee[i][i]; /*self terms */
    for(int i=0;i<nw-1;i++){                      /*cross terms*/
      x=2*w[i];
      y=0;
      for(int j=i+1;j<nw;j++)y+=w[j]*ee[i][j];
      z+=x*y;
    }
    return z;    
  }

  /*<R*E>_d =< R * Sum_i(w[i]*E_i) >_d == Sum_i( w[i] * <R*E_i>_d )*/
  double RE(const vecd &w) const {
    double z=0;
    for(int i=0;i<nw;i++) z+=w[i]*Re[i];
    return z;
  }

  /*correlation coefficient between R and E*/
  double r(const vecd &w) const {
    double avE=this->E(w); /*cout<<"avE="<<avE<<endl; cout<<"E2="<<this->E2(w)<<endl;*/
    double dE2=this->E2(w)-avE*avE; /*cout<<"dE="<<sqrt(dE2)<<endl;*/
    double dR2=R2-R*R; /*cout<<"dR="<<sqrt(dR2)<<endl;*/
    double ther=( this->RE(w)- R*avE )/sqrt( dE2 * dR2 ); /*cout<<"r="<<ther<<"\n\n";*/
    /*if(dE2<0){
      this->printQuery();
      exit(1);
      }*/
    return ther;
  }

  /*Z-score. We take it as positive*/
  double Z(const vecd &w) const {
    double avE=this->E(w);
    double dE2=this->E2(w)-avE*avE;    
    return ( avE - this->En(w) )/sqrt(dE2);
  }

};/*Matches class query*/

typedef vector<query> qvec; /*a shorthand*/

/*======================================================================*/
/*class specifiying the set and functions to minimize*/
class set2minimize{
private:
  qvec queries; /*a vector container of all queries*/
  double r0;
  double Z0;
public:

  set2minimize(const qvec &theQueries, const double &r, const double &Z): queries(theQueries), r0(r), Z0(Z) {}
  
 /*native energy averaged over all native estates*/
  double avEn(const vecd &w){
    double en=0.0;
    static const int N=queries.size();
    for(int i=0;i<N;i++) en+=queries[i].En(w);
    return en;
  }

  double avr(const vecd &w){
    double av_r=0;
    static const int N=queries.size(); /*initialize only the first time G1 is invoked*/
    for(int i=0;i<N;i++) av_r+=queries[i].r(w);
    return av_r/N;
  }

  double avZ(const vecd &w){
    double av_Z=0;
    static const int N=queries.size(); /*initialize only the first time G1 is invoked*/
    for(int i=0;i<N;i++) av_Z+=queries[i].Z(w);
    return av_Z/N;
  }

  double G1(const vecd &w) const {
    double av_r=0;
    double g1;
    static const int N=queries.size(); /*initialize only the first time G1 is invoked*/
    for(int i=0;i<N;i++) av_r+=queries[i].r(w);    /*cout<<"av_r= "<<av_r/N;exit(1);*/
    g1=1/( 1+av_r/(r0*N) );
    return g1;
  }

  double G3(const vecd &w) const {
    double av_Z=0;
    double g3;
    static const int N=queries.size(); /*initialize only the first time G2 is invoked*/
    for(int i=0;i<N;i++) av_Z+=queries[i].Z(w);
    av_Z/=N; //printf("av_Z=%lf ",av_Z);
    g3=1/( 1+av_Z/Z0 );
    return g3;
  }

  double G1G3(const vecd &w) const {
    double g1=this->G1(w);
    double g3=this->G3(w);
    return g1*g3;
  }

  double minimize(const vecd &w) const{ /*selected function to minimize*/
    double g1=this->G1(w);
    double g3=this->G3(w);
    double g1g3=g1*g3;
    //printf("%5.3lf %5.3lf %5.3lf %6.3lf \n",g1,1.0/g1-1,g3,1.0/g3-1);
    double x=g1g3;
    return x;
  }

};
/*======================================================================*/
/*An instance of the following class is a neccesary container to pass to Minuit*/
class FNCDerived : public FCNBase{
private:
  set2minimize setx;
  double theErrorDef;

public:
  FNCDerived(const qvec &theQueries, const double &r, const double &Z): setx(theQueries,r,Z), theErrorDef(1.0) {}
  ~FNCDerived(){}

  virtual double up() const {return theErrorDef;}

  /*  virtual double operator()(const vecd &w ) const{*/
  virtual double operator()(const vecd &w ) const {
    assert(w.size()==nw);
    return  setx.minimize(w);
  }

};/*Matches class FNCDerived : public FCNBase*/

/*======================================================================*/
/*welcoming message*/
int wellcome(){
  system("clear");
  printf("Usage: ./minimze.x [options]\n");
  printf("  Required:\n");
  printf("  -a input file containing info on queries\n");
  printf("  -b input file containing the weights\n");
  printf("  -c input r0 (def:r0=1.0)\n");
  printf("  -d input Z0 (def:Z0=1.0)\n");
  printf("\n");
  return 1 ; /*failure!*/
}


/*======================================================================*/
/*read all atributes from a file pointer. File format assumes that all
  info for a query resideds in a single line in the following order:
  header R R2 Re[nw] en[nw] e[nw] ee[nw,nw]
  Comment lines begin with '#'
*/
void initQueries(qvec &theQueries, char *inpf){
  ifstream in(inpf);
  string line;
  int n=0;
  getline(in,line);
  while(line.length()>0){
    if(line[0]=='#'){ getline(in,line);  continue; }/*comment line*/
    theQueries.push_back( query(line) ); /*insert a query object*/
    getline(in,line);
  }
}/*Matches void initQueries(..)*/

/*======================================================================*/
void initW(vecd &w,char *inpw){ /*pass the weights from file to memory*/
  ifstream in(inpw);
  vector<string> tokens;
  string line;
  getline(in,line);
  while(line.length()>0){
    if(line[0]=='#'){ getline(in,line); continue; } /*comment line*/
    splitString(line," ",tokens,false);
    getline(in,line);
  }
  for(int i=0; i<tokens.size(); i++) w.push_back( atof(tokens[i].c_str()) );
}
/*======================================================================*/
void initErr(const vecd &w, vecd &err){ /*estimate uncertainities in the weights*/
  double x=0.005;
  for(int i=0; i<nw; i++) err.push_back(x*w[i]); /*estimate a x*100% error*/
}
/*======================================================================*/
/*EXECUTION STARTS HERE*/
int main(int argc, char **argv){

  char *inpf=NULL; /*pointer to file name of energy averages*/
  char *inpw=NULL; /*pointer to file name of weigths*/
  char *name= new char[16];
  char *buf= new char[512];
  double r0=1.0;
  double Z0=1.0;
  qvec theQueries; /*store all info for the queries*/
  vecd w0,w;           /*store the weights*/
  vecd err;         /*store uncertainities for the weights*/

  {/*BLOCK TO PARSE INPUT*/
    int option_char,required=2,nreq=0;
    extern char *optarg;
    while ((option_char = getopt(argc, argv, ":a:b:c:d:h")) != EOF){
      switch (option_char){ 
      case 'a': inpf=optarg; nreq++; break;
      case 'b': inpw=optarg; nreq++; break;
      case 'c': r0=atof(optarg); break;
      case 'd': Z0=atof(optarg); break;
      case 'h': return wellcome(); 
      default: return wellcome();
      }
    }
    if(nreq<required) return wellcome(); /*user did not supply all required arguments*/
  }

  initW(w0,inpw);               /*pass the weights from file to memory*/
  initQueries(theQueries,inpf); /*pass all info in file to memory*/
  initErr(w0,err);              /*estimate uncertainities in the weights*/

  //set2minimize setx(theQueries,r0,Z0);
  //cout<<setx.G1(w)<<"\n"<<setx.G3(w);

  FNCDerived f(theQueries,r0,Z0); /*create appropriate object that can be minimized*/


  /*Minuit minimization*/
  MnUserParameters mw;
  for(int i=0;i<nw;i++){
    sprintf(name,"%02d",i);
    mw.add(name,w0[i],err[i]);
    //mw.setLimits(i,0.0,10.0);
  }
  MnMigrad theMinimizer(f,mw); /*create minimizer (default constructor)*/
  FunctionMinimum min=theMinimizer(); /*run minimization procedure*/
  

  /*output results*/
  double lambda,r,Z;
  if( min.isValid() ){
    for(int i=0;i<nw;i++) w.push_back( theMinimizer.value(i) );
    lambda=fabs(setx.avEn(w0)/setx.avEn(w));
    /*cout<<setx.avr(w0)<<" "<<setx.avZ(w0)<<endl;*/
    r=setx.avr(w);
    Z=setx.avZ(w);
    sprintf(buf,"results= %5.3lf %5.3lf %5.3lf %5.3lf %6.2lf",r0,Z0,r,Z,lambda);
    for(int i=0;i<nw;i++){
      w[i]*=lambda; //normalize so that average native energies coincide
      sprintf(buf,"%s %5.3lf ",buf,w[i]);
    }
    printf("%s\n",buf);
  }
  

  return 0;
}/*Matches main(..)*/
