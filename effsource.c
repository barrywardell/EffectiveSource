/*******************************************************************************
 * Copyright (C) 2011 Barry Wardell
 *
 * Authors: Barry Wardell <effectivesource@barrywardell.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA 02110-1301, USA.
 *
 * Version 2.2 - 24 February 2011
 *
 * This code compute the singular field and effective source for a point scalar
 * particle following a circular geodesic orbit in Kerr.
 *
 * Definitions: x is the 4-dimensional location of the field point, xp is the
 *        location of the particle. M is the mass of the background Kerr black
 *        hole. a = J/M is the spin of the black hole in units of M.
 *
 * Usage: Call effsource_init(M, a, x_p) every time M, a or x_p changes.
 *        In this circular orbit case, we only use x_p.r and this doesn't
 *        change so these will be a constant for the whole simulation and
 *        effsouce_init only needs to be called once when the program starts.
 *
 *        effsource_phis(x, x_p) returns the value of the singular field.
 *
 *        effsource_calc(x, phis, dphis_dr, dphis_dth, dphis_dph, box_phis)
 *        computes the singular field, its spatial derivatives and its
 *        d'Alembertian and stores them in the variables phis, dphis_dr,
 *        dphis_dth, dphis_dph and box_phis.
 *
 * Options: There is currently only one option, set by the static variable
 *        'periodic'. Setting this to one rewrites powers of dphi in terms
 *        of trigonemetric functions. This enforces periodicity in phi
 *        and avoids problems which could crop up when evaluating the effective
 *        source at the opposite side of the black hole, for example.
 */

#include <math.h>
#include "effsource.h"

/* Whether to enforce periodicity in phi. 1=yes, 0=no */
static int periodic = 1;

/* The particle's coordinate location and 4-velocity */
static struct coordinate xp, up;

/* Mass and spin of the Kerr black hole */
static double M, a;

/* Static variables used to store the coefficients of the series expansions */
static double A_0_0_2, A_0_0_4, A_0_2_0, A_0_2_2, A_0_4_0, A_1_0_2, A_1_0_4, A_1_2_0, A_1_2_2, A_1_4_0, A_2_0_0, A_2_0_2, A_2_2_0, A_3_0_0, A_3_0_2, A_3_2_0, A_4_0_0, A_5_0_0;
static double s2_0_0_2, s2_0_0_4, s2_0_2_0, s2_0_2_2, s2_0_4_0, s2_1_0_2, s2_1_0_4, s2_1_2_0, s2_1_2_2, s2_1_4_0, s2_2_0_0, s2_2_0_2, s2_2_2_0, s2_3_0_0, s2_3_0_2, s2_3_2_0, s2_4_0_0, s2_5_0_0;

/* Compute the singular field at the point x for the particle at xp */
void effsource_phis(struct coordinate * x, double * phis)
{
  double A, s2;

  double dr       = x->r - xp.r;
  double dtheta   = x->theta - xp.theta;
  double dphi     = x->phi - xp.phi;
  
  double dr2      = dr*dr;
  double dr3      = dr2*dr;
  double dr4      = dr2*dr2;
  double dr5      = dr3*dr2;
  double dtheta2  = dtheta*dtheta;
  double dtheta4  = dtheta2*dtheta2;
  double dphi2, dphi4;
  
  if(periodic)
  {
    double cosdphi  = cos(dphi);
    double cos2dphi = cos(2*dphi);
    
    dphi2 = 2.5 - (8*cosdphi)/3. + cos2dphi/6.;
    dphi4 = 6 - 8*cosdphi + 2*cos2dphi;
  } else {
    dphi2 = dphi*dphi;
    dphi4 = dphi2*dphi2;
  }
  
  A = A_0_0_2*dphi2 + A_0_0_4*dphi4 + A_0_2_0*dtheta2
    + A_0_2_2*dphi2*dtheta2 + A_0_4_0*dtheta4 + A_1_0_2*dr*dphi2
    + A_1_0_4*dr*dphi4 + A_1_2_0*dr*dtheta2 + A_1_2_2*dr*dtheta2*dphi2 
    + A_1_4_0*dr*dtheta4 + A_2_0_0*dr2 + A_2_0_2*dr2*dphi2
    + A_2_2_0*dr2*dtheta2 + A_3_0_0*dr3 + A_3_0_2*dr3*dphi2
    + A_3_2_0*dr3*dtheta2 + A_4_0_0*dr4 + A_5_0_0*dr5;

  s2 = s2_0_0_2*dphi2 + s2_0_0_4*dphi4 + s2_0_2_0*dtheta2 
    + s2_0_2_2*dtheta2*dphi2 + s2_0_4_0*dtheta4 + s2_1_0_2*dr*dphi2
    + s2_1_0_4*dr*dphi4 + s2_1_2_0*dr*dtheta2 + s2_1_2_2*dr*dtheta2*dphi2
    + s2_1_4_0*dr*dtheta4 + s2_2_0_0*dr2 + s2_2_0_2*dr2*dphi2
    + s2_2_2_0*dr2*dtheta2 + s2_3_0_0*dr3 + s2_3_0_2*dr3*dphi2
    + s2_3_2_0*dr3*dtheta2 + s2_4_0_0*dr4 + s2_5_0_0*dr5;

  *phis = A/(24.*pow(s2, 1.5));
}

/* Compute the singular field, its derivatives and its d'Alembertian */
void effsource_calc(struct coordinate * x, double *phis, double *dphis_dr,
      double *dphis_dth, double *dphis_dph, double *dphis_dt, double *box_phis)
{
  double A, dA_dr, d2A_dr2, dA_dth, d2A_dth2, dA_dph, d2A_dph2;
  double s2, sqrts2, s2_15, s2_25, s2_35, ds2_dr, d2s2_dr2, ds2_dth, d2s2_dth2, ds2_dph, d2s2_dph2;
  double d2phis_dr2, d2phis_dth2, d2phis_dph2, d2phis_dt2, d2phis_dphdt;

  double r      = x->r;
  double theta  = x->theta;
  double phi    = x->phi;
  double rp     = xp.r;
  double thetap = xp.theta;
  double phip   = xp.phi;
  
  double dr     = r-rp;
  double dtheta = theta - thetap;
  double dphi   = phi - phip;
  double om     = M / (a*M + sqrt(M*pow(rp,3)));
  double om2    = om*om;
  double a2     = a*a;

  double dr2      = dr*dr;
  double dr3      = dr2*dr;
  double dr4      = dr2*dr2;
  double dr5      = dr3*dr2;
  double dtheta2  = dtheta*dtheta;
  double dtheta3  = dtheta2*dtheta;
  double dtheta4  = dtheta2*dtheta2;
  double dphi2, dphi4;
  double dphi2_d, dphi2_d_d, dphi4_d, dphi4_d_d;
  
  if(periodic)
  {
    double cosdphi  = cos(dphi);
    double cos2dphi = cos(2*dphi);
    double sindphi  = sin(dphi);
    double sin2dphi = sin(2*dphi);
    
    dphi2 = 2.5 - (8*cosdphi)/3. + cos2dphi/6.;
    dphi4 = 6 - 8*cosdphi + 2*cos2dphi;

    dphi2_d   = (8*sindphi)/3. - sin2dphi/3.;
    dphi4_d   = 8*sindphi - 4*sin2dphi;

    dphi2_d_d = (8*cosdphi)/3. - (2*cos2dphi)/3.;
    dphi4_d_d = 8*cosdphi - 8*cos2dphi;
  } else {
    dphi2     = dphi*dphi;
    dphi2_d   = 2.0*dphi;
    dphi2_d_d = 2.0;
    
    dphi4     = dphi2*dphi2;
    dphi4_d   = 4.0*dphi2*dphi;
    dphi4_d_d = 12.0*dphi2;
  }

  /* A, dA/dr, d^2A/dr^2 */
  A = A_0_0_2*dphi2 + A_0_0_4*dphi4 + A_0_2_0*dtheta2
    + A_0_2_2*dphi2*dtheta2 + A_0_4_0*dtheta4 + A_1_0_2*dr*dphi2
    + A_1_0_4*dr*dphi4 + A_1_2_0*dr*dtheta2 + A_1_2_2*dr*dtheta2*dphi2 
    + A_1_4_0*dr*dtheta4 + A_2_0_0*dr2 + A_2_0_2*dr2*dphi2
    + A_2_2_0*dr2*dtheta2 + A_3_0_0*dr3 + A_3_0_2*dr3*dphi2
    + A_3_2_0*dr3*dtheta2 + A_4_0_0*dr4 + A_5_0_0*dr5;
  dA_dr = A_1_0_2*dphi2 + A_1_0_4*dphi4 + A_1_2_0*dtheta2
    + A_1_2_2*dtheta2*dphi2 + A_1_4_0*dtheta4 + 2.*A_2_0_0*dr
    + 2.*A_2_0_2*dr*dphi2 + 2.*A_2_2_0*dr*dtheta2 + 3.*A_3_0_0*dr2
    + 3.*A_3_0_2*dr2*dphi2 + 3.*A_3_2_0*dr2*dtheta2 + 4.*A_4_0_0*dr3
    + 5.*A_5_0_0*dr4;
  d2A_dr2 = 2.*A_2_0_0 + 2.*A_2_0_2*dphi2 + 2.*A_2_2_0*dtheta2
    + 6.*A_3_0_0*dr + 6.*A_3_0_2*dr*dphi2 + 6.*A_3_2_0*dr*dtheta2 
    + 12.*A_4_0_0*dr2 + 20.*A_5_0_0*dr3;
  dA_dth = 2.*A_0_2_0*dtheta + 2.*A_0_2_2*dphi2*dtheta + 4.*A_0_4_0*dtheta3
    + 2.*A_1_2_0*dr*dtheta + 2.*A_1_2_2*dr*dtheta*dphi2 
    + 4.*A_1_4_0*dr*dtheta3 + 2.*A_2_2_0*dr2*dtheta + 2.*A_3_2_0*dr3*dtheta;
  d2A_dth2 = 2.*A_0_2_0 + 2.*A_0_2_2*dphi2 + 12.*A_0_4_0*dtheta2
    + 2.*A_1_2_0*dr + 2.*A_1_2_2*dr*dphi2 
    + 12.*A_1_4_0*dr*dtheta2 + 2.*A_2_2_0*dr2 + 2.*A_3_2_0*dr3;
  dA_dph = A_0_0_2*dphi2_d + A_0_0_4*dphi4_d + A_0_2_2*dphi2_d*dtheta2
    + A_1_0_2*dr*dphi2_d + A_1_0_4*dr*dphi4_d + A_1_2_2*dr*dtheta2*dphi2_d
    + A_2_0_2*dr2*dphi2_d + A_3_0_2*dr3*dphi2_d;
  d2A_dph2 = A_0_0_2*dphi2_d_d + A_0_0_4*dphi4_d_d + A_0_2_2*dphi2_d_d*dtheta2
    + A_1_0_2*dr*dphi2_d_d + A_1_0_4*dr*dphi4_d_d + A_1_2_2*dr*dtheta2*dphi2_d_d
    + A_2_0_2*dr2*dphi2_d_d + A_3_0_2*dr3*dphi2_d_d;

  /* s, ds/dr, d^2s/dr^2 */
  s2 = s2_0_0_2*dphi2 + s2_0_0_4*dphi4 + s2_0_2_0*dtheta2 
    + s2_0_2_2*dtheta2*dphi2 + s2_0_4_0*dtheta4 + s2_1_0_2*dr*dphi2
    + s2_1_0_4*dr*dphi4 + s2_1_2_0*dr*dtheta2 + s2_1_2_2*dr*dtheta2*dphi2
    + s2_1_4_0*dr*dtheta4 + s2_2_0_0*dr2 + s2_2_0_2*dr2*dphi2
    + s2_2_2_0*dr2*dtheta2 + s2_3_0_0*dr3 + s2_3_0_2*dr3*dphi2
    + s2_3_2_0*dr3*dtheta2 + s2_4_0_0*dr4 + s2_5_0_0*dr5;
  sqrts2 = sqrt(s2);
  s2_15  = s2*sqrts2;
  s2_25  = s2*s2_15;
  s2_35  = s2*s2_25;
  ds2_dr = s2_1_0_2*dphi2 + s2_1_0_4*dphi4 + s2_1_2_0*dtheta2
    + s2_1_2_2*dtheta2*dphi2 + s2_1_4_0*dtheta4 + 2.*s2_2_0_0*dr
    + 2.*s2_2_0_2*dr*dphi2 + 2.*s2_2_2_0*dr*dtheta2 + 3.*s2_3_0_0*dr2
    + 3.*s2_3_0_2*dr2*dphi2 + 3.*s2_3_2_0*dr2*dtheta2 + 4.*s2_4_0_0*dr3
    + 5.*s2_5_0_0*dr4;
  d2s2_dr2 = 2.*s2_2_0_0 + 2.*s2_2_0_2*dphi2 + 2.*s2_2_2_0*dtheta2
    + 6.*s2_3_0_0*dr + 6.*s2_3_0_2*dr*dphi2 + 6.*s2_3_2_0*dr*dtheta2 
    + 12.*s2_4_0_0*dr2 + 20.*s2_5_0_0*dr3;
  ds2_dth = 2.*s2_0_2_0*dtheta + 2.*s2_0_2_2*dphi2*dtheta + 4.*s2_0_4_0*dtheta3
    + 2.*s2_1_2_0*dr*dtheta + 2.*s2_1_2_2*dr*dtheta*dphi2 
    + 4.*s2_1_4_0*dr*dtheta3 + 2.*s2_2_2_0*dr2*dtheta + 2.*s2_3_2_0*dr3*dtheta;
  d2s2_dth2 = 2.*s2_0_2_0 + 2.*s2_0_2_2*dphi2 + 12.*s2_0_4_0*dtheta2
    + 2.*s2_1_2_0*dr + 2.*s2_1_2_2*dr*dphi2 
    + 12.*s2_1_4_0*dr*dtheta2 + 2.*s2_2_2_0*dr2 + 2.*s2_3_2_0*dr3;
  ds2_dph = s2_0_0_2*dphi2_d + s2_0_0_4*dphi4_d + s2_0_2_2*dphi2_d*dtheta2
    + s2_1_0_2*dr*dphi2_d + s2_1_0_4*dr*dphi4_d + s2_1_2_2*dr*dtheta2*dphi2_d
    + s2_2_0_2*dr2*dphi2_d + s2_3_0_2*dr3*dphi2_d;
  d2s2_dph2 = s2_0_0_2*dphi2_d_d + s2_0_0_4*dphi4_d_d + s2_0_2_2*dphi2_d_d*dtheta2
    + s2_1_0_2*dr*dphi2_d_d + s2_1_0_4*dr*dphi4_d_d + s2_1_2_2*dr*dtheta2*dphi2_d_d
    + s2_2_0_2*dr2*dphi2_d_d + s2_3_0_2*dr3*dphi2_d_d;

  /* phis */
  *phis = A/(24.*s2_15);

  /* First derivatives of phis */
  *dphis_dr  = (-3*ds2_dr*A + 2*dA_dr*s2) /(48.*s2_25);
  *dphis_dth = (-3*ds2_dth*A + 2*dA_dth*s2)/(48.*s2_25);
  *dphis_dph = (-3*ds2_dph*A + 2*dA_dph*s2)/(48.*s2_25);
  *dphis_dt  = -om * (*dphis_dph);
  
  /* Second derivatives of phis */
  d2phis_dr2 =
    (15*ds2_dr*ds2_dr*A - 6*s2*(2*dA_dr*ds2_dr + d2s2_dr2*A) + 4*d2A_dr2*s2*s2)/(96.*s2_35);
  d2phis_dth2 = 
    (15*ds2_dth*ds2_dth*A - 6*s2*(2*dA_dth*ds2_dth + d2s2_dth2*A) + 4*d2A_dth2*s2*s2)/(96.*s2_35);
  d2phis_dph2 = 
    (15*ds2_dph*ds2_dph*A - 6*s2*(2*dA_dph*ds2_dph + d2s2_dph2*A) + 4*d2A_dph2*s2*s2)/(96.*s2_35);
  d2phis_dt2  = om*om*d2phis_dph2;
  d2phis_dphdt = -om*d2phis_dph2;
  
  double sinth = sin(theta);
  double costh = cos(theta);
  double cotanth = costh/sinth;
  double r2 = r*r;

  /* Box[phis]. Note that because we have a circular geodeisc, d/dth = - om d/dph. */
  *box_phis = -(d2phis_dph2*om2) 
    + (d2phis_dth2 + (*dphis_dr)*(-2*M + 2*r) + d2phis_dr2*(a2 + r*(-2*M + r)) + (*dphis_dth)*cotanth + 
        (d2phis_dph2*(a2*cotanth*cotanth + r*(-2*M*om*(-2*a + a2*om + om*r2)
        + (-2*M + r)/(sinth*sinth))))/(a2 + r*(-2*M + r))
      ) / (r2 + a2*costh*costh);

  /* We can't currently get an accurate result for dx < 0.1.
   * Just fudge it for now and set the source to 0 in this case.
   */
  if (sqrt(dr*dr + dtheta*dtheta + dphi*dphi) < 0.1)
  {
    *box_phis = 0.;
  }
}

/* Initialize array of coefficients of powers of dr, dtheta and dphi. */
void effsource_init(double mass, double spin)
{
  M = mass;
  a = spin;
}

void effsource_set_particle(struct coordinate * x_p, struct coordinate * u_p)
{
  xp = *x_p;
  up = *u_p;
  double r1 = xp.r;

  /* Compute A ccoefficients */
  {
    double v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31, v32, v33, v34, v35, v36, v37, v38, v39, v40, v41, v42, v43, v44, v45, v46, v47, v48, v49, v50, v51, v52, v53, v54, v55, v56;
    v0 = pow(M,5);
    v1 = pow(a,5);
    v2 = pow(r1,4);
    v3 = pow(M,4);
    v4 = pow(a,6);
    v5 = pow(a,7);
    v6 = pow(M,3);
    v7 = pow(M,2);
    v8 = 2*r1;
    v9 = pow(a,4);
    v10 = 5*M;
    v11 = pow(r1,5);
    v12 = pow(r1,-2);
    v13 = 3*M;
    v14 = pow(a,3);
    v15 = pow(r1,3);
    v16 = pow(a,2);
    v17 = M*r1;
    v18 = pow(r1,2);
    v19 = 1/r1;
    v20 = 2*v6;
    v21 = 2*v7;
    v22 = -17*v7;
    v23 = 8*v15;
    v24 = pow(r1,4.5);
    v25 = pow(r1,-2.5);
    v26 = pow(r1,5.5);
    v27 = 9*v18;
    v28 = M - r1;
    v29 = -v18;
    v30 = pow(r1,2.5);
    v31 = pow(M,2.5);
    v32 = -M + r1;
    v33 = pow(r1,3.5);
    v34 = r1 + v13;
    v35 = sqrt(r1);
    v36 = pow(M,1.5);
    v37 = pow(r1,1.5);
    v38 = sqrt(M);
    v39 = pow(r1,-5.5);
    v40 = -2*M + r1;
    v41 = M*v16;
    v42 = 2*M*v18;
    v43 = pow(v17,2.5);
    v44 = pow(v17,1.5);
    v45 = v13 + v8;
    v46 = sqrt(v17);
    v47 = v16 + r1*v40;
    v48 = 2*v30*v38;
    v49 = pow(v47,-2);
    v50 = 1/v47;
    v51 = 2*M*v45*v9;
    v52 = -3*v17 + v18 + 2*a*v46;
    v53 = pow(v52,-2);
    v54 = v15 + v41 + a*v46*v8;
    v55 = 1/v52;
    v56 = 1/(r1*v13 + v29 - 2*a*v46);
    A_0_0_2 = 24*v19*v47*v54*v55;
    A_0_0_4 = 2*pow(v37 + a*v38,2)*v39*(-(v33*v34) + 2*v14*v36 + 6*a*v18*v36 - 3*v34*v35*v41)*v47*v55;
    A_0_2_0 = 24*v18;
    A_0_2_2 = -12*v12*((r1 + v10)*v15*v16 + v11*v32 + 2*a*(-4*M + r1)*v33*v38 + v14*(-8*v31*v35 - 4*v44 + v48) + v51)*v55;
    A_0_4_0 = 2*v56*(3*r1*(-r1 + v10)*v16 - 6*v14*v46 + a*(-12*v44 + v48) + v18*(-3*v17 + v18 + 6*v7));
    A_1_0_2 = -24*v19*v28*v54*v55;
    A_1_0_4 = (v39*(-60*pow(a,9)*pow(M,3.5) + 12*pow(a,8)*(36*M - 5*r1)*v35*v6 - 3*v37*v4*v7*(-v15 - 4*M*v18 + 72*v6 + 125*r1*v7) + v26*v41*(143*M*v15 - 56*v2 + 306*v3 - 243*r1*v6 - 54*v18*v7) - 6*a*pow(r1,7)*v38*(4*M*v15 + 2*v2 + 2*v3 + 59*r1*v6 - 49*v18*v7) - 2*v14*v2*v36*(-181*M*v15 + 70*v2 + 168*v3 + 78*r1*v6 - 39*v18*v7) + pow(r1,8.5)*(-7*M*v15 - 2*v2 + 240*v3 - 309*r1*v6 + 114*v18*v7) + 2*v1*v18*v36*(-13*M*v15 - 12*v2 + 738*v3 - 558*r1*v6 + 169*v18*v7) + v31*v5*(227*v17 + v27 - 380*v7)*v8 - M*v30*(894*v0 + 9*v11 + 110*M*v2 - 1518*r1*v3 + 588*v18*v6 - 155*v15*v7)*v9))/pow(v52,3);
    A_1_2_0 = 24*r1;
    A_1_2_2 = -3*v25*v53*(-6*(10*M + 3*r1)*v1*v36 + v37*v41*(33*v15 - 74*M*v18 + 42*v6 - 21*r1*v7) + 2*v24*(2*v15 + M*v29 + 12*v6 - 18*r1*v7) + 2*a*v15*v38*(-10*M*v18 + v23 + 6*v6 + 9*r1*v7) + v14*v36*(v17 + 21*v18 - 60*v7)*v8 + M*v35*(16*v17 - 9*v18 + 151*v7)*v9);
    A_1_4_0 = (v13*v16 + v18*v45 + 4*a*(-3*M + r1)*v46)*v56;
    A_2_0_0 = 24*v18*v50;
    A_2_0_2 = 12*v12*v50*v55*(2*M*v2*v40 - 2*a*(v33*v36 + v43) - 2*v14*(2*v31*v35 - v30*v38 + v44) + v51 + r1*v16*(v15 + v10*v18 - 7*v6 + r1*v7));
    A_2_2_0 = 12*v16*v50;
    A_3_0_0 = 24*r1*(v16 - v17)*v49;
    A_3_0_2 = -3*v25*v49*v53*(6*(6*M + r1)*v36*v5 + M*v35*v4*(4*v17 + 3*v18 - 89*v7) - 2*a*v2*v36*(v15 + v42 + 24*v6 - 26*r1*v7) + M*v26*(v15 + v42 + 30*v6 - 23*r1*v7) + v16*v30*(-132*v0 + 4*v11 - 11*M*v2 + 236*r1*v3 - 155*v18*v6 + 20*v15*v7) + 2*v14*v18*v38*(-19*M*v15 + 8*v2 - 32*v3 + 8*r1*v6 + 26*v18*v7) + v1*v36*(22*v17 + v27 - 42*v7)*v8 + M*v37*(27*v15 - 26*M*v18 + 306*v6 - 197*r1*v7)*v9);
    A_3_2_0 = 3*v49*v55*(r1*v16*(17*v17 - 4*v18 + v22) + M*v18*(v17 + v21 + v29) + 8*v14*v28*v46 + 2*M*v9);
    A_4_0_0 = (6*v55*(-4*v16*v18*(-4*v17 + v18 + v21) - 8*v14*v30*v38 + 2*a*(2*v33*v36 + v43) + 2*v1*v46 + M*v15*(-13*v17 + 4*v18 + 5*v7) + r1*v32*v9))/pow(v47,3);
    A_5_0_0 = (3*v55*(v13*v4 + 4*a*(-(v24*v36) + M*v43 - 2*r1*v43) + 4*v14*(4*v33*v38 + M*v44 - r1*v44) + 4*(M - 3*r1)*v1*v46 - M*v15*(6*v15 - 19*M*v18 + v20 + 6*r1*v7) + v16*v18*(-31*M*v18 + v20 + v23 + 14*r1*v7) + r1*(22*v17 - 6*v18 + v22)*v9))/pow(v47,4);
  }

  /* Compute s2 ccoefficients */
  {
    double v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31, v32, v33, v34, v35, v36, v37, v38, v39, v40, v41, v42, v43, v44, v45, v46;
    v0 = -6*r1;
    v1 = pow(r1,-3);
    v2 = 5*M;
    v3 = pow(a,5);
    v4 = 19*M;
    v5 = pow(a,6);
    v6 = pow(M,3);
    v7 = 9*r1;
    v8 = 11*M;
    v9 = pow(r1,-2);
    v10 = 3*r1;
    v11 = pow(a,3);
    v12 = pow(M,2);
    v13 = pow(r1,4);
    v14 = 2*M;
    v15 = pow(a,4);
    v16 = -2*M;
    v17 = pow(r1,3);
    v18 = pow(a,2);
    v19 = M*r1;
    v20 = pow(r1,2);
    v21 = 1/r1;
    v22 = pow(r1,1.5);
    v23 = pow(M,1.5);
    v24 = 12*v12;
    v25 = pow(r1,4.5);
    v26 = 3*v18;
    v27 = sqrt(r1);
    v28 = pow(M,2.5);
    v29 = 3*v17;
    v30 = 3*v20;
    v31 = pow(r1,2.5);
    v32 = pow(r1,3.5);
    v33 = sqrt(M);
    v34 = 2*v12;
    v35 = M*v18;
    v36 = pow(v19,1.5);
    v37 = -r1 + v14;
    v38 = sqrt(v19);
    v39 = 5*v36;
    v40 = r1*(r1 + v16) + v18;
    v41 = -3*v31*v33;
    v42 = 8*v27*v28;
    v43 = pow(v40,-2);
    v44 = 1/v40;
    v45 = v17 + v35 + 2*a*r1*v38;
    v46 = 1/(-3*v19 + v20 + 2*a*v38);
    s2_0_0_2 = v21*v40*v45*v46;
    s2_0_0_4 = ((M*v15 - r1*v18*(-v19 + v20 + v34) + v13*v37)*((M + r1)*v13 + (-5*M + v10)*v18*v19 + a*(2*v32*v33 - 4*r1*v36) + v11*v14*v38)*v46)/(12.*pow(r1,6));
    s2_0_2_0 = v20;
    s2_0_2_2 = -(v46*(2*a*v31*v33*(-11*v19 + v30 + v34) + v13*(-5*v19 + v30 + v34) - 2*v11*(v39 + v41 + v42) + r1*v18*(-3*r1*v12 + 10*M*v20 + v29 + 2*v6) + M*v15*(v7 + v8))*v9)/6.;
    s2_0_4_0 = (v26 + r1*v37)/12.;
    s2_1_0_2 = (-M + r1)*v21*v45*v46;
    s2_1_0_4 = (v46*(-(pow(r1,7.5)*(-9*v12 + 4*v19 + v20)) - 2*v11*v17*(6*v12 + r1*v16 + v20)*v23 - 6*pow(a,7)*v28 + 2*r1*(2*r1 + v2)*v28*v3 + (r1*v2 - 9*v20 + v24)*v25*v35 + 2*a*pow(r1,7)*v33*v37 + v12*v27*(-3*r1 + v4)*v5 - M*v15*v22*(-17*r1*v12 + 9*M*v20 + v29 + 33*v6)))/(12.*pow(r1,6.5));
    s2_1_2_0 = r1;
    s2_1_2_2 = (v1*v46*(v18*v19*(M*v10 - 17*v20 + v24) - 2*v13*(-3*v12 + r1*v14 + v30) - 2*v11*(16*v27*v28 + v39) + M*v15*(29*M + v7) + 2*a*v32*v33*(v0 + v8)))/12.;
    s2_1_4_0 = (M - r1)/12.;
    s2_2_0_0 = v20*v44;
    s2_2_0_2 = (v44*v46*(M*(-11*M + 5*r1)*v13 - 2*a*(5*pow(v19,2.5) + v23*v32) - 2*v11*(4*v27*v28 + v36 + v41) + M*v15*(13*M + v7) + r1*v18*(v10*v12 + v29 - 17*v6 + v20*v8))*v9)/6.;
    s2_2_2_0 = ((-v19 + v26)*v44)/6.;
    s2_3_0_0 = r1*(v18 - v19)*v43;
    s2_3_0_2 = (v1*v43*v46*(v15*v19*(70*v12 - 27*v19 - 11*v20) + M*pow(r1,5)*(8*v12 - 5*v19 + v20) + a*(-20*pow(v19,3.5) + 20*v25*v28 + 2*v13*v36) + 2*v3*(v36 + v42) - M*(v10 + v4)*v5 - 4*v11*v22*v33*(v0*v12 - 4*M*v20 + v29 + 12*v6) + v18*v20*(-46*pow(M,4) - 6*v13 + v17*v2 - 25*v12*v20 + 58*r1*v6)))/12.;
    s2_3_2_0 = ((v18*(v0 + v2) + M*v20)*v43)/12.;
    s2_4_0_0 = (3*v15 + 2*r1*(M + v0)*v18 - M*(M - 8*r1)*v20)/(12.*pow(v40,3));
    s2_5_0_0 = ((4*M - 9*r1)*v15 - M*v20*(v12 + r1*v16 + 6*v20) + r1*v18*(3*v12 - 5*v19 + 12*v20))/(12.*pow(v40,4));
  }
}