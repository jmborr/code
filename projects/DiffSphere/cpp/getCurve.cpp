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


void getCurve(diffSphere &dS, const double &wi, const double &wf, const double &dw, std::vector<double> &dataX, std::vector<double> &dataY){

  //populate the frequencies
  double x = wi;
  while(x < wf){
    dataX.push_back(x);
    x += dw;
  }

  //calculate the structure factor dataY for the frequencies of dataX
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
  //std::cout << "\nNorm of Lorentzian Coefficients=" << norm <<"\n";

  for(std::vector<double>::const_iterator itx=dataX.begin(); itx!=dataX.end(); ++itx){
    double x = *itx;
    double y = 0.0;  //structure factor, sum of all the Lorentzians
    size_t i=0; //coefficient counter
    std::vector<double>::const_iterator itYJ=YJ.begin();
    //loop over all coefficients
    for(std::vector<xnlc>::const_iterator it=xnl.begin(); it!=xnl.end(); ++it){
      double cw = (it->x * it->x) * w;
      double L = cw/(cw*cw+x*x); //Lorentzian. L is dependent on parameter w and data x
      y +=  (*itYJ)*L ;  //update the structure factor value
      //printf("%2d %5.2f %1d %5.2f %7.2f %14.12f %7.5f %10.8f\n",i,it->x,l,fabs(a - it->x),Y,J,L,Y*J*L);
      //std::cout<<"i="<<i<<" x="<<it->x<<" l="<<l<<" fabs(a - it->x)="<< fabs(a - it->x)<<" Y="<<Y<<" J="<<J<<" L="<<L<<" y="<<Y*J*L<<"\n";
      ++itYJ;
      ++i;
    } // end of for(std::vector<xnlc>::const_iterator it=xnl.begin()
    dataY.push_back( y );
  } // end of for(std::vector<xnlc>::const_iterator itx=dataX.begin(); itx!=dataX.end(); ++itx)

} // end of getCurve

int main(int argc, char* argv[]){

  //parse command line arguments
  try{
    TCLAP::CmdLine cmd("getCurve produces S(Q,w) for diffusion within a sphere without the elastic term.\nExample: ./getCurve -o outfile.dat -i 0.0 -f 10.0 -s 0.1 -q 0.6 -d 3.0 -r 2.5", ' ', "0.0");
    TCLAP::ValueArg<double> radius("r","radius","sphere radius,",true,1.0,"double");
    cmd.add( radius );
    TCLAP::ValueArg<double> diffusion("d","diffusion","diffusion constant,",true,1.0,"double");
    cmd.add( diffusion );
    TCLAP::ValueArg<double> momTransf("q","momTransf","momentum transfer,",true,1.0,"double");
    cmd.add( momTransf );
    TCLAP::ValueArg<double> wi("i","initFreq","initial frequency,",true,0.0,"double");
    cmd.add( wi );
    TCLAP::ValueArg<double> wf("f","finalFreq","final frequency,",true,1.0,"double");
    cmd.add( wf );
    TCLAP::ValueArg<double> ws("s","spFreq","frequency step,",true,0.01,"double");
    cmd.add( ws );
    TCLAP::ValueArg<std::string> outfile("o","outFile","file to write output,",true,"junk.dat","string");
    cmd.add( outfile );
    cmd.parse( argc, argv );

    double R = radius.getValue();
    double D = diffusion.getValue();
    double Q = momTransf.getValue();
    double iw = wi.getValue();
    double fw = wf.getValue();
    double sw = ws.getValue();

    diffSphere dS(R,D,Q);
    std::vector<double> X,Y;
    getCurve(dS,iw,fw,sw,X,Y);

    std::string fileName=outfile.getValue();
    std::ofstream fout(fileName.c_str());
    fout << "#R="<<R<<"\n#D="<<D<<"\n#Q="<<Q<<"\n#frequency  S(Q,w)\n";
    std::vector<double>::const_iterator itY=Y.begin();
    for(std::vector<double>::const_iterator itX=X.begin(); itX!=X.end(); ++itX){
      fout << *itX <<" "<< *itY << "\n";
      ++itY;
    }

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return 0;
}

