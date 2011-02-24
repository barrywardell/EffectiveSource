#include <stdio.h>
#include <math.h>
#include "phis.h"

int main(int argc, char* argv[])
{
  double val;
  
  const double a = 0.5;
  const double M = 1.0;

  const double r1     = 9.;
  const double theta1 = M_PI_2;
  const double phi1   = 0.;

  phis_init(M, a, r1);
  
  double theta = M_PI_2;

//   for(double r=4.0; r<=14.0; r+=0.1)
//   {
//     for(double phi=-M_PI; phi<=M_PI; phi+=0.1)
//     {
//       val = phis(r, r1, theta, theta1, phi, phi1);
//       
//       printf("%g\t%g\t%g\t%g\n", r, theta, phi, val);
//     }
//   }

  double r = 10.0;
  double phi = 0.0;
  double vall = 0.;
  for(int i=0; i<10000000; i++)
  {
    val = phis(r, r1, theta, theta1, phi, phi1);
    vall += val;
  }
  
  printf("%g\n", vall);

  return 0;
}