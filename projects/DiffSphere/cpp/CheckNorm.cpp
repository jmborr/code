/*
 * getCurve.cpp
 *
 *  Created on: Jan 28, 2012
 *      Author: Jose Borreguero
 */

#include <fstream>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include "diffSphere.h"
#include <boost/math/special_functions/bessel.hpp>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void CheckNorm(diffSphere &dS){

  double R=dS.getR();
  double D=dS.getD();
  double w=D/(R*R); //active parameters, "scaled" radius and diffusion
  std::vector<xnlc> xnl=dS.getxnl();

  std::vector<double>YJ = dS.LorentzianCoefficients();

  //check normalization of the Lorentzian Coefficients();
  double norm=0.0;
  for(std::vector<double>::const_iterator itYJ=YJ.begin(); itYJ!=YJ.end(); ++itYJ){
    norm += *itYJ;
  }
  double a=dS.getQ()*R;
  double z= 3*boost::math::sph_bessel(1,a)/a;
  norm = z*z + M_PI * norm;
  std::cout << "Norm of Lorentzian Coefficients = " << norm <<"\n";

} // end of CheckNorm

int main(int argc, char* argv[]){

  //parse command line arguments
  try{
    TCLAP::CmdLine cmd("CheckNorm adds up all the Lorentzian coefficients plus the elastic term.\nExample: ./CheckNorm -q 0.6 -d 3.0 -r 2.5", ' ', "0.0");
    TCLAP::ValueArg<double> radius("r","radius","sphere radius,",true,1.0,"double");
    cmd.add( radius );
    TCLAP::ValueArg<double> diffusion("d","diffusion","diffusion constant,",true,1.0,"double");
    cmd.add( diffusion );
    TCLAP::ValueArg<double> momTransf("q","momTransf","momentum transfer,",true,1.0,"double");
    cmd.add( momTransf );
    cmd.parse( argc, argv );

    double R = radius.getValue();
    double D = diffusion.getValue();
    double Q = momTransf.getValue();

    diffSphere dS(R,D,Q);
    CheckNorm(dS);

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return 0;
}

