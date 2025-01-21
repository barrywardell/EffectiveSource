/*******************************************************************************
 * Copyright (C) 2011 Barry Wardell
 ******************************************************************************/
#include <stdio.h>

#include "effsource.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>

int counter_phi, counter_theta;
enum component{Re, Im};

/* Parameters for the integration */
struct intpars
{
  struct coordinate x;
  int l;
  int m;
  enum component comp;
  double (*func)(struct coordinate *x);
};

/* Compute the integrand in the spherical harmonic decomposition */
double lm_integrand (double phi, void * par)
{
  struct intpars *ip = (struct intpars *) par;
  double theta  = ip->x.theta;
  int l = ip->l;
  int m = ip->m;
  
  /* Compute the value of the function to be decomposed */
  double val;
  ip->x.phi = phi;
  val = ip->func(&(ip->x));

  counter_phi++;

  /* Return either the real or imaginary component of the integrand */
  if(ip->comp == Re)
  {
    return val * gsl_sf_legendre_sphPlm(l, m, cos(theta)) * sin(theta) * cos(m*phi);
  } else {
    return -val * gsl_sf_legendre_sphPlm(l, m, cos(theta)) * sin(theta) * sin(m*phi);
  }
}

double lm_decompose_phi (double theta, void * pars)
{
  double res, err;
  struct intpars ip;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  ip = *((struct intpars *)pars);
  ip.x.theta = theta;

  F.params   = &ip;
  F.function = &lm_integrand;

  gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-7, 1000, w, &res, &err);
  gsl_integration_workspace_free (w);

  counter_theta++;

  return res;
}

void lm_decompose(int l, int m, const double r,
  double (*func)(struct coordinate * x), double * res_re, double * res_im)
{
  double err_re, err_im;
  struct intpars ip;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  ip.x.r  = r;
  ip.l    = l;
  ip.m    = m;
  ip.func = func;
  
  F.params   = &ip;
  F.function = &lm_decompose_phi;

  ip.comp = Re;
  gsl_integration_qags (&F, 0., M_PI, 0, 1e-7, 1000, w, res_re, &err_re);

  ip.comp = Im;
  gsl_integration_qags (&F, 0., M_PI, 0, 1e-7, 1000, w, res_im, &err_im);

  gsl_integration_workspace_free (w);

  return;
}

/* Compute the integrand in the spherical harmonic decomposition */
double m_integrand (double phi, void * par)
{
  struct intpars *ip = (struct intpars *) par;
  int m = ip->m;
  
  /* Compute the value of the function to be decomposed */
  double val;
  ip->x.phi = phi;
  val = ip->func(&(ip->x));

  counter_phi++;

  /* Return either the real or imaginary component of the integrand */
  if(ip->comp == Re)
  {
    return val * cos(m*phi);
  } else {
    return -val * sin(m*phi);
  }
}

void m_decompose(int m, struct coordinate x,
  double (*func)(struct coordinate * x), double * res)
{
  double err_re, err_im;
  struct intpars ip;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  ip.x    = x;
  ip.l    = 0;
  ip.m    = m;
  ip.func = func;
  
  F.params   = &ip;
  F.function = &m_integrand;

  ip.comp = Re;
  gsl_integration_qags (&F, -M_PI, M_PI, 0, 1e-10, 1000, w, &res[0], &err_re);

  ip.comp = Im;
  gsl_integration_qags (&F, -M_PI, M_PI, 0, 1e-10, 1000, w, &res[1], &err_im);

  gsl_integration_workspace_free (w);

  return;
}

