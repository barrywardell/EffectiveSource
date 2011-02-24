#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phis.h"
#include "phisr9.h"
#include "src.h"
#include "srcr9.h"

int main(int argc, char* argv[])
{
  double val;
  const double r_p = 9.;
  const double a = 0.5;
  const double M = 1.0;
  int type = atoi(argv[1]);
  int num_iters = atoi(argv[2]);

  double r = 10.0;
  double phi = 0.0;
  double phip = 0.0;
  double vall = 0.;
  double theta = M_PI_2;
  double thetap = M_PI_2;
  Initialize_src();
  for(int i=0; i<10000000; i++)
  {
    src(r_p, a, M, r-r_p, theta - thetap, phi-phip, &val);
    vall += val;
  }
  Uninitialize_src();
  printf("%g\n", vall);
//   switch(type) {
//   case 1:
//     Initialize_phis();
//     for(int i=0; i<num_iters; i++)
//     {
//       double dr = ((double)rand())/RAND_MAX;
//       double dth = ((double)rand())/RAND_MAX;
//       double dph = ((double)rand())/RAND_MAX;
//   
//       phis(r_p, a, M, dr, dth, dph, &val);
//     }
//     Uninitialize_phis();
//     break;
//     
//   case 2:
//     Initialize_phisr9();
//     for(int i=0; i<num_iters; i++)
//     {
//       double dr = ((double)rand())/RAND_MAX;
//       double dth = ((double)rand())/RAND_MAX;
//       double dph = ((double)rand())/RAND_MAX;
//   
//       phisr9(dr, dth, dph, &val);
//     }
//     Uninitialize_phisr9();
//     break;
// 
//   case 3:
//     Initialize_src();
//     for(int i=0; i<num_iters; i++)
//     {
//       double dr = ((double)rand())/RAND_MAX;
//       double dth = ((double)rand())/RAND_MAX;
//       double dph = ((double)rand())/RAND_MAX;
//   
//       src(r_p, a, M, dr, dth, dph, &val);
//     }
//     Uninitialize_src();
//     break;
// 
//   case 4:
//     Initialize_srcr9();
//     for(int i=0; i<num_iters; i++)
//     {
//       double dr = ((double)rand())/RAND_MAX;
//       double dth = ((double)rand())/RAND_MAX;
//       double dph = ((double)rand())/RAND_MAX;
//   
//       srcr9(dr, dth, dph, &val);
//     }
//     Uninitialize_srcr9();
//     break;
//   
//     case 5:
//     Initialize_srcr9();
//     for(double dph=-2*M_PI; dph<2*M_PI; dph+=0.01)
//     {
//       double dr = 0;
//       double dth = 0;
//   
//       srcr9(dr, dth, dph, &val);
//       printf("%g\n", val);
//     }
//     Uninitialize_srcr9();
//     break;
// 
//     case 6:
//     Initialize_src();
//     for(double dph=-2*M_PI; dph<2*M_PI; dph+=0.01)
//     {
//       double dr = 0;
//       double dth = 0;
//   
//       src(r_p, a, M, dr, dth, dph, &val);
//       printf("%g\n", val);
//     }
//     Uninitialize_src();
//     break;
//   }

  return 0;
}