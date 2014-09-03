/*
 * diffSphere.h
 *
 *  Created on: Jan 28, 2012
 *      Author: Jose Borreguero
 */

#ifndef DIFFSPHERE_H_
#define DIFFSPHERE_H_

#include <string>
#include <vector>

using namespace std;

//structure to hold info on Volino's coefficients
struct xnlc{
  double x;
  size_t l;
  size_t n;
};

// simple structure to hold a linear interpolation of factor J around its numerical divergence point
struct linearJ{
  double slope;
  double intercept;
};

class diffSphere{
  public:
  /// Constructor
  diffSphere(const double &radius, const double &diffusion, const double &mom_transf);
  /// Destructor
  ~diffSphere() {}

  private:
    double R;
    double D;
    double Q;
    std::vector<xnlc> xnl;      //xnl coefficients
    std::vector<double> alpha; //certain coefficients invariant during fitting
    size_t lmax;         //maximum value of l in xnlist
    size_t ncoeff;      //number of coefficients
    double divZone;    //linear interpolation zone around the numerical divergence of factor J
    std::vector<linearJ> linearJlist;
    void initXnlCoeff();
    void initAlphaCoeff();
    void initLinJlist();

  public:
    std::string name()const{return "DiffSphere";}
    double getR(){return R;}
    double getD(){return D;}
    double getQ(){return Q;}
    size_t getlmax(){return lmax;}
    double getdivZone(){return divZone;}
    std::vector<xnlc> const getxnl(){ return xnl; }
    std::vector<double> const getalpha(){ return alpha; };
    std::vector<linearJ> const getlinearJlist(){ return linearJlist; };
    std::vector<double> LorentzianCoefficients();
};

#endif /* DIFFSPHERE_H_ */
