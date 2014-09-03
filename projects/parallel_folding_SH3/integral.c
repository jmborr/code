/*obtain the integral from -inifinity to x of the standar gaussian distribution*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv){
  double x;
  x=atof(argv[1]);
  if(x<0) printf( "%e\n", 0.5*(1-erf(-x)));
  else printf( "%e\n", 0.5*(1+erf(x)) );
  return 0;
}
