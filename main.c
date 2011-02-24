#include <stdio.h>
#include <math.h>
#include "effsource.h"

int main(int argc, char* argv[])
{
  struct coordinate x, xp;
  
  const double a = 0.5;
  const double M = 1.0;

  /* The particle's coordinate location */
  xp.r     = 9.;
  xp.theta = M_PI_2;
  xp.phi   = 0.0;
  xp.t     = 0.0;

  /* The field point */
  x.theta   = M_PI_2;
  x.t       = 0.;

  /* Initialize the coefficients */
  effsource_init(M, a, &xp);

  for(double r=4.0; r<=14.0; r+=0.1)
  {
    x.r = r;
    for(double phi=-M_PI; phi<=M_PI; phi+=0.1)
    {
      double phis, dphis_dr, dphis_dtheta, dphis_dphi, dphis_dt, src;
      x.phi     = phi;
      effsource_calc(M, a, &x, &xp, &phis, &dphis_dr, &dphis_dtheta, &dphis_dphi, &dphis_dt, &src);
      
      printf("%g\t%g\t%g\t%g\t%g\n", x.r, x.theta, x.phi, phis, src);
    }
  }

  return 0;
}