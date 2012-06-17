#include <stdio.h>
#include <math.h>
#include "effsource.h"

double phis_calc(struct coordinate * x)
{
  double res;
  effsource_phis(x, &res);
  return res;
}

double src_calc(struct coordinate * x)
{
  double res1, res2, res3, res4, res5, res6;
  effsource_calc(x, &res1, &res2, &res3, &res4, &res5, &res6);
  return res6;
}

int main(int argc, char* argv[])
{
  struct coordinate x, xp, up;

  /* Mass and spin of the central black hole */
  const double a = 0.5;
  const double a2  = a*a;
  const double a3  = a2*a;
  const double a4  = a2*a2;
  const double M = 1.0;

  /* Energy an angular momentum for an elliptic orbit between r=9 and r=11 (M=1) */
  double e  = sqrt((-434070 + 2471*a2 + 6*sqrt(110)*a*sqrt(6237 + 162*a2 + a4))/(-474721 + 3960*a2));
  double l  = ((261*a + a3 - 3*sqrt(110)*sqrt(6237 + 162*a2 + a4))*
              sqrt((-434070 + 2471*a2 + 6*sqrt(110)*a*sqrt(6237 + 162*a2 + a4))/(-474721 + 3960*a2)))/(-630 + a2);

  /* 4-velocity for the elliptic orbit - choose ingoing motion */
  double rp  = 10.0;
  double rp2 = rp*rp;
  double rp3 = rp2*rp;
  double e2  = e*e;
  double ur2 = -1 + e2 - 2*(-(l-a*e)*(l-a*e)*M/rp3 + (l*l-a2*(e2-1))/(2.*rp2) - M/rp);

  up.t   = (-2*a*l*M + e*rp3 + a2*e*(2*M + rp))/(rp*a2 + rp*(-2*M + rp));
  up.phi = (2*a*e*M - 2*l*M + l*rp)/(a2*rp - 2*M*rp2 + rp3);
  up.r   = ur2 > 0.0 ? -sqrt(ur2) : 0.0;

  /* Circular orbit */
//   up.t   = (rp*sqrt(M*rp)+a*M)/(sqrt(M*rp)*sqrt(rp*rp-3.0*M*rp+2.0*a*sqrt(M*rp)));
//   up.phi = sqrt(M*rp)/(rp*sqrt(rp*rp-3.0*M*rp+2.0*a*sqrt(M*rp)));
//   up.r   = 0;
//   e = ((rp-2.0*M)*sqrt(M*rp)+a*M)/(sqrt(M*rp)*sqrt(rp*rp-3.0*M*rp+2.0*a*sqrt(M*rp)));
//   l = (M*(a2+rp2-2.0*a*sqrt(M*rp)))/(sqrt(M*rp)*sqrt(rp*rp-3.0*M*rp+2.0*a*sqrt(M*rp)));

  xp.t = 0;
  xp.r = rp;
  xp.theta = M_PI_2;
  xp.phi = 0;

  /* The point where we measure the effective source */
  x.t     = 0.0;
  x.r     = rp;
  x.theta = M_PI_2;
  x.phi   = 0.0;

  /* Initialize the background parameters */
  effsource_init(M, a);

  int use_el = 1;
  if( use_el )
  {
    effsource_set_particle_el(&xp, e, l, up.r);
  } else {
    effsource_set_particle(&xp, &up);
  }

  int m = 2;
  for(double r=9.9; r<=10.1; r+=0.01)
  {
    x.r = r;
    for(double theta=M_PI_2-0.1; theta<=M_PI_2+0.1; theta+=0.011)
    {
      double phis, dphis_dr, dphis_dtheta, dphis_dphi, dphis_dt, src, phis_num_re, phis_num_im;
      x.theta     = theta;
      effsource_calc_m(m, &x, &phis, &dphis_dr, &dphis_dtheta, &dphis_dphi, &dphis_dt, &src);
      m_decompose(m, x, phis_calc, &phis_num_re, &phis_num_im);
      
      printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
        x.r-xp.r, x.theta-xp.theta, x.phi-xp.phi,
        phis, dphis_dr, dphis_dtheta, dphis_dphi, dphis_dt, src, phis_num_re, phis_num_im);
    }
  }

  return 0;
}