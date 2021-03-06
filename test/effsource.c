#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "effsource.h"
#include "decompose.h"

/* Return the value of the singular field at the point x */
double PhiS_calc(struct coordinate * x)
{
  double PhiS;
  effsource_PhiS(x, &PhiS);
  return PhiS;
}

/* Return the value of the effective source at the point x */
double src_calc(struct coordinate * x)
{
  double PhiS, dPhiS_dx[4], d2PhiS_dx2[10], src;
  effsource_calc(x, &PhiS, dPhiS_dx, d2PhiS_dx2, &src);
  return src;
}

int main(int argc, char* argv[])
{
  int orbit;

  if(argc != 2)
  {
    printf( "usage: %s orbit", argv[0] );
    return(0);
  } else {
    orbit = atoi(argv[1]);
  }

  /* Mass and spin of the central black hole */
  const double a = 0.5;
  const double M = 1.0;
  effsource_init(M, a);

  /* Orbital parameters */
  double r_p  = 10.0;
  double e, l, ur;
  struct coordinate xp;

  if(atoi(argv[1]) == 0)
  {
    /* Energy an angular momentum for an elliptic orbit between r=9 and r=11 (M=1) */
    e = sqrt((-434070 + 2471*a*a + 6*sqrt(110)*a*sqrt(6237 + 162*a*a + a*a*a*a))/(-474721 + 3960*a*a));
    l = ((261*a + a*a*a - 3*sqrt(110)*sqrt(6237 + 162*a*a + a*a*a*a))*
        sqrt((-434070 + 2471*a*a + 6*sqrt(110)*a*sqrt(6237 + 162*a*a + a*a*a*a))/(-474721 + 3960*a*a)))/(-630 + a*a);
    ur = -sqrt(-1 + e*e - 2*(-(l-a*e)*(l-a*e)*M/(r_p*r_p*r_p) + (l*l-a*a*(e*e-1))/(2.*r_p*r_p) - M/r_p));
  } else {
    /* Circular orbit of radius 10M */
    e = ((r_p-2.0*M)*sqrt(M*r_p)+a*M)/(sqrt(M*r_p)*sqrt(r_p*r_p-3.0*M*r_p+2.0*a*sqrt(M*r_p)));
    l = (M*(a*a+r_p*r_p-2.0*a*sqrt(M*r_p)))/(sqrt(M*r_p)*sqrt(r_p*r_p-3.0*M*r_p+2.0*a*sqrt(M*r_p)));
    ur = 0;
  }

  xp.t = 0;
  xp.r = r_p;
  xp.theta = M_PI_2;
  xp.phi = 0;

  effsource_set_particle(&xp, e, l, ur);

  /* The point where we measure the singular field/effective source */
  struct coordinate x = {0.0, r_p, M_PI_2, 0.0};

  /* Disable the GSL error handler so that it doesn't abort due to roundoff errors */
  gsl_set_error_handler_off();

  /* Output the singular field for the m=2 mode in the r-theta plane */
  int m = 2;
  for(double r=9.9; r<=10.1; r+=0.01)
  {
    x.r = r;
    for(double theta=M_PI_2-0.1; theta<=M_PI_2+0.1; theta+=0.011)
    {
      double PhiS[2], src_num[2], dPhiS[8], ddPhiS[20], src[2];
      x.theta     = theta;
      effsource_calc_m(m, &x, PhiS, dPhiS, ddPhiS, src);
      m_decompose(m, x, src_calc, src_num);

      printf("%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
        x.r-xp.r, x.theta-xp.theta, x.phi-xp.phi,
        src[0], src[1], src_num[0], src_num[1]);
    }
  }

  return 0;
}