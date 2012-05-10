/*******************************************************************************
 * Copyright (C) 2012 Barry Wardell
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
 * Version 3.0 alpha 1 - 10 May 2012
 *
 * This code compute the singular field and effective source for a point scalar
 * particle following an equatorial geodesic orbit in Kerr.
 *
 * Definitions: x is the 4-dimensional location of the field point, xp is the
 *        location of the particle. M is the mass of the background Kerr black
 *        hole. a = J/M is the spin of the black hole in units of M.
 *
 * Usage: Call effsource_init(M, a) at startup to set the mass and spin of the
 *        central black hole.
 *
 *        Call either effsource_set_particle(x_p, u_p) or
 *        effsource_set_particle_el(x_p, e, l, ur_p) every time the particle's
 *        position (x_p), velocity (u_p, ur_p) or constants of motion (e, l) change.
 *
 *        effsource_phis(x, *phis) stores the value of the singular field in phis
 *
 *        effsource_calc(x, *phis, *dphis_dr, *dphis_dth, *dphis_dph, *box_phis)
 *        computes the singular field, its spatial derivatives and its
 *        d'Alembertian and stores them in the variables phis, dphis_dr,
 *        dphis_dth, dphis_dph and box_phis.
 *
 */

#include <math.h>
#include <assert.h>
#include "effsource.h"
#include <stdio.h>

/* The particle's coordinate location and 4-velocity */
static struct coordinate xp;

/* Mass and spin of the Kerr black hole */
static double M, a;

/* Static variables used to store the coefficients of the series expansions */
static double A006, A008, A024, A026, A042, A044, A060, A062, A080, A106, A108, A124, A126, A142, A144, A160, A162, A180, A204, A206, A222, A224, A240, A242, A260, A304, A306, A322, A324, A340, A342, A360, A402, A404, A420, A422, A440, A502, A504, A520, A522, A540, A600, A602, A620, A700, A702, A720, A800, A900;
static double alpha20, alpha02, beta;

/* Compute the singular field at the point x for the particle at xp */
void effsource_phis(struct coordinate * x, double * phis)
{

  double A, alpha, rho2;

  double r      = x->r;
  double theta  = x->theta;
  double phi    = x->phi;
  double rp     = xp.r;
  double thetap = xp.theta;
  double phip   = xp.phi;

  double dr     = r - rp;
  double dtheta = theta - thetap;
  double dphi   = phi - phip;

  double dr2 = dr*dr;
  double dr3 = dr2*dr;
  double dr4 = dr2*dr2;
  double dr6 = dr3*dr3;
  double dr8 = dr4*dr4;

  double dtheta2  = dtheta*dtheta;
  double dtheta4  = dtheta2*dtheta2;
  double dtheta8  = dtheta4*dtheta4;
  
  double sindphi  = sin(dphi);
  double sindphi2 = sindphi*sindphi;
  double sindphi4 = sindphi2*sindphi2;
  double sindphi8 = sindphi4*sindphi4;

  A = dr6*(A600 + A700*dr) + dr8*(A800 + A900*dr) + dr4*(A420 + A520*dr + dr2*(A620 + A720*dr))*dtheta2 + 
   (A080 + A180*dr)*dtheta8 + dtheta4*(dr2*(A240 + A340*dr) + dr4*(A440 + A540*dr) + 
      (A060 + A160*dr + dr2*(A260 + A360*dr))*dtheta2) + 
   (dr4*(A402 + A502*dr + dr2*(A602 + A702*dr)) + (dr2*(A222 + A322*dr) + dr4*(A422 + A522*dr))*dtheta2 + 
      dtheta4*(A042 + A142*dr + dr2*(A242 + A342*dr) + (A062 + A162*dr)*dtheta2))*sindphi2 + 
   (A008 + A108*dr)*sindphi8 + sindphi4*(dr2*(A204 + A304*dr) + dr4*(A404 + A504*dr) + 
      (A024 + A124*dr + dr2*(A224 + A324*dr))*dtheta2 + (A044 + A144*dr)*dtheta4 + 
      (A006 + A106*dr + dr2*(A206 + A306*dr) + (A026 + A126*dr)*dtheta2)*sindphi2);

  alpha = alpha20*dr2 + alpha02*dtheta2;
  rho2 = alpha + beta*sindphi2;

  *phis = A/pow(rho2, 3.5);
}

/* TODO: This function is just copied from the old version and needs to be updated.
   Compute the singular field, its derivatives and its d'Alembertian */
void effsource_calc(struct coordinate * x, double *phis, double *dphis_dr,
      double *dphis_dth, double *dphis_dph, double *dphis_dt, double *box_phis)
{
  double A, dA_dr, d2A_dr2, dA_dth, d2A_dth2, dA_dQ, dA_dR, dA_dph, d2A_dQ2, d2A_dR2, d2A_dQdR, d2A_dph2, dA_dt, d2A_dt2, d2A_dQdt, d2A_dRdt, d2A_dphdt;
  double s2, sqrts2, s2_15, s2_25, s2_35, ds2_dr, d2s2_dr2, ds2_dth, d2s2_dth2, ds2_dQ, ds2_dR, ds2_dph, d2s2_dQ2, d2s2_dR2, d2s2_dQdR, d2s2_dph2, ds2_dt, d2s2_dt2, d2s2_dQdt, d2s2_dRdt, d2s2_dphdt;
  double d2phis_dr2, d2phis_dth2, d2phis_dph2, d2phis_dt2, d2phis_dphdt;

  double r      = x->r;
  double theta  = x->theta;
  double phi    = x->phi;
  double rp     = xp.r;
  double thetap = xp.theta;
  double phip   = xp.phi;

  double dr     = r - rp;
  double dtheta = theta - thetap;
  double dphi   = phi - phip;

  /* A, dA/dx, d^2A/dx^2 */
  A         = 0*dr + 0*dtheta + 0*dphi;
  dA_dr     = 0;
  dA_dth    = 0;
  dA_dQ     = 0;
  dA_dR     = 0;
  dA_dph    = 0;
  dA_dt     = 0;
  d2A_dr2   = 0;
  d2A_dth2  = 0;
  d2A_dQ2   = 0;
  d2A_dR2   = 0;
  d2A_dQdR  = 0;
  d2A_dph2  = 0;
  d2A_dt2   = 0;
  d2A_dQdt  = 0;
  d2A_dRdt  = 0;
  d2A_dphdt = 0;

  /* s, ds/dr, d^2s/dr^2 */
  s2         = 0;
  ds2_dr     = 0;
  ds2_dth    = 0;
  ds2_dQ     = 0;
  ds2_dR     = 0;
  ds2_dph    = 0;
  ds2_dt     = 0;
  d2s2_dr2   = 0;
  d2s2_dth2  = 0;
  d2s2_dQ2   = 0;
  d2s2_dR2   = 0;
  d2s2_dQdR  = 0;
  d2s2_dph2  = 0;
  d2s2_dt2   = 0;
  d2s2_dQdt  = 0;
  d2s2_dRdt  = 0;
  d2s2_dphdt = 0;
  sqrts2     = sqrt(s2);
  s2_15      = s2*sqrts2;
  s2_25      = s2*s2_15;
  s2_35      = s2*s2_25;

  /* phis */
  *phis = A/(24.*s2_15);

  /* First derivatives of phis */
  *dphis_dr  = (-3*ds2_dr*A + 2*dA_dr*s2) /(48.*s2_25);
  *dphis_dth = (-3*ds2_dth*A + 2*dA_dth*s2)/(48.*s2_25);
  *dphis_dph = (-3*ds2_dph*A + 2*dA_dph*s2)/(48.*s2_25);
  *dphis_dt  = (-3*ds2_dt*A + 2*dA_dt*s2)/(48.*s2_25);
  
  /* Second derivatives of phis */
  d2phis_dr2 =
    (15*ds2_dr*ds2_dr*A - 6*s2*(2*dA_dr*ds2_dr + d2s2_dr2*A) + 4*d2A_dr2*s2*s2)/(96.*s2_35);
  d2phis_dth2 =
    (15*ds2_dth*ds2_dth*A - 6*s2*(2*dA_dth*ds2_dth + d2s2_dth2*A) + 4*d2A_dth2*s2*s2)/(96.*s2_35);
  d2phis_dph2 =
    (15*ds2_dph*ds2_dph*A - 6*s2*(2*dA_dph*ds2_dph + d2s2_dph2*A) + 4*d2A_dph2*s2*s2)/(96.*s2_35);
  d2phis_dt2 =
    (15*ds2_dt*ds2_dt*A - 6*s2*(2*dA_dt*ds2_dt + d2s2_dt2*A) + 4*d2A_dt2*s2*s2)/(96.*s2_35);
  d2phis_dphdt =
    (15*ds2_dph*ds2_dt*A - 6*s2*(dA_dt*ds2_dph + dA_dph*ds2_dt + d2s2_dphdt*A) + 4*d2A_dphdt*s2*s2)/(96.*s2_35);
  
  
  /* Box[phis] */
  double sinth  = sin(theta);
  double sinth2 = sinth*sinth;
  double sin2th = sin(2*theta);
  double cos2th = cos(2*theta);
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double a2 = a*a;
  double a4 = a2*a2;

  *box_phis = -((2*a2*(*dphis_dr) - a2*d2phis_dph2 - a4*d2phis_dr2 - a2*d2phis_dth2 - 4*(*dphis_dr)*r - 2*a2*(*dphis_dr)*r + 
         4*d2phis_dph2*r + 2*a*d2phis_dphdt*r + 4*a2*d2phis_dr2*r + 2*d2phis_dth2*r + 6*(*dphis_dr)*r2 - 2*d2phis_dph2*r2 - 4*d2phis_dr2*r2 - 
         2*a2*d2phis_dr2*r2 - d2phis_dth2*r2 - 2*(*dphis_dr)*r3 + 4*d2phis_dr2*r3 - d2phis_dr2*r4 + 
         (a4*d2phis_dt2 + 4*a*d2phis_dphdt*r + 2*d2phis_dt2*r4 + a2*d2phis_dt2*r*(2 + 3*r))*sinth2 + 
         cos2th*(a4*d2phis_dr2 - 2*a*d2phis_dphdt*r + (-2 + r)*r*(d2phis_dth2 + 2*(*dphis_dr)*(-1 + r) + d2phis_dr2*(-2 + r)*r) + 
            a2*(-d2phis_dph2 + d2phis_dth2 + 2*(*dphis_dr)*(-1 + r) - 4*d2phis_dr2*r + 2*d2phis_dr2*r2) + 
            a2*d2phis_dt2*(a2 + (-2 + r)*r)*sinth2) - a2*(*dphis_dth)*sin2th + 2*(*dphis_dth)*r*sin2th - 
         (*dphis_dth)*r2*sin2th))/((sinth2*(a2 + (-2 + r)*r)*(a2 + 2*r2 + a2*cos2th)));
}

/* Initialize array of coefficients of pows of dr, dtheta and dphi. */
void effsource_init(double mass, double spin)
{
  M = mass;
  a = spin;
}

/* This is not implemented or needed */
void effsource_set_particle(struct coordinate * x_p, struct coordinate * u_p)
{
  assert(0);
}

/* Initialize array of coefficients of pows of dr, dtheta and dphi. */
void effsource_set_particle_el(struct coordinate * x_p, double e, double l, double ur_p)
{
  xp = *x_p;
  double r1 = xp.r;

  /* Compute A coefficients */
  {
    A006 = pow(pow(l,2) + pow(r1,2) + (pow(a,2)*(2*M + r1))/r1,3);
    A008 = -(pow(pow(a,2)*(2*M + r1)
    	+ r1*(pow(l,2) + pow(r1,2)),2)*(-(pow(a,8)*pow(M,2)) + 4*pow(a,7)*e*l*pow(M,2) +
    	4*a*e*l*M*(2*M - 3*r1)*pow(r1,6) - 4*pow(a,3)*e*l*M*pow(r1,3)*(22*pow(M,2) + 3*M*r1 +
    	4*pow(r1,2)) - 4*pow(a,5)*e*l*M*(16*pow(M,3) + 18*pow(M,2)*r1 + pow(r1,3)) +
    	2*pow(a,6)*M*(-2*pow(l,2)*M + 2*pow(e,2)*pow(2*M + r1,3) + r1*(2*pow(M,2) - M*r1 +
    	pow(r1,2))) + pow(r1,6)*(8*pow(l,2)*M*(-2*M + r1) + pow(r1,2)*(-4*pow(M,2) + 4*(-1 +
    	pow(e,2))*M*r1 + 3*pow(r1,2))) + 2*pow(a,2)*pow(r1,3)*(4*pow(l,2)*M*(4*pow(M,2) - 3*M*r1 +
    	2*pow(r1,2)) + pow(r1,2)*(4*pow(M,3) + 12*(-1 + pow(e,2))*pow(M,2)*r1 + 3*(1 +
    	2*pow(e,2))*M*pow(r1,2) + 3*pow(r1,3))) + pow(a,4)*(4*pow(l,2)*M*(8*pow(M,3) +
    	6*pow(M,2)*r1 - 3*M*pow(r1,2) + 2*pow(r1,3)) + pow(r1,2)*(-4*pow(M,4) + 3*(-3 +
    	16*pow(e,2))*pow(M,2)*pow(r1,2) + 12*(1 + pow(e,2))*M*pow(r1,3) + 3*pow(r1,4) +
    	4*pow(M,3)*(r1 + 12*pow(e,2)*r1)))))/(24.*pow(r1,8)*(pow(a,2) + r1*(-2*M + r1)));
    A024 =
    	3*pow(pow(a,2)*(2*M + r1) + r1*(pow(l,2) + pow(r1,2)),2);
    A026 = -((pow(a,2)*(2*M + r1) +
    	r1*(pow(l,2) + pow(r1,2)))*(-2*pow(a,7)*e*l*M*(24*pow(M,2) + 16*M*r1 + 3*pow(r1,2)) +
    	2*a*e*l*M*pow(r1,5)*(pow(l,2)*(10*M - 9*r1) + (14*M - 15*r1)*pow(r1,2)) -
    	2*pow(a,5)*e*l*M*r1*(16*pow(M,3) + 56*pow(M,2)*r1 + 50*M*pow(r1,2) + 17*pow(r1,3) +
    	3*pow(l,2)*(4*M + r1)) + pow(a,8)*M*(6*pow(e,2)*pow(2*M + r1,2) - r1*(7*M + 3*r1)) +
    	pow(r1,5)*(4*M*pow(r1,4)*(-2*M + r1 + 2*pow(e,2)*r1) + pow(l,2)*pow(r1,2)*(-36*pow(M,2) +
    	4*(9 + pow(e,2))*M*r1 - 9*pow(r1,2)) - 2*pow(l,4)*(8*pow(M,2) - 10*M*r1 + 3*pow(r1,2))) +
    	2*pow(a,3)*e*l*M*pow(r1,2)*(2*pow(l,2)*(4*pow(M,2) - 7*M*r1 - 6*pow(r1,2)) -
    	pow(r1,2)*(16*pow(M,2) + 28*M*r1 + 29*pow(r1,2))) +
    	pow(a,2)*pow(r1,2)*(pow(l,2)*pow(r1,2)*(-4*pow(M,3) + 4*(6 + pow(e,2))*pow(M,2)*r1 + (33 +
    	14*pow(e,2))*M*pow(r1,2) - 18*pow(r1,3)) + pow(l,4)*(-8*pow(M,3) + 12*pow(M,2)*r1 +
    	8*M*pow(r1,2) - 6*pow(r1,3)) + pow(r1,4)*(4*pow(M,3) + 36*pow(e,2)*pow(M,2)*r1 + 3*(3 +
    	10*pow(e,2))*M*pow(r1,2) - 3*pow(r1,3))) + pow(a,4)*r1*(12*pow(l,4)*pow(M,2) +
    	pow(l,2)*(16*pow(M,4) - 8*(-5 + pow(e,2))*pow(M,3)*r1 + 16*(2 +
    	pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-3 + 8*pow(e,2))*M*pow(r1,3) - 9*pow(r1,4)) +
    	pow(r1,2)*(4*pow(M,4) + 16*(2 + 3*pow(e,2))*pow(M,3)*r1 + 3*(3 +
    	32*pow(e,2))*pow(M,2)*pow(r1,2) + 3*(-3 + 14*pow(e,2))*M*pow(r1,3) - 6*pow(r1,4))) +
    	pow(a,6)*(pow(l,2)*M*(24*pow(M,2) + 4*(2 + 3*pow(e,2))*M*r1 + 3*(-1 +
    	2*pow(e,2))*pow(r1,2)) + r1*(-(r1*(-12*pow(M,3) + 18*pow(M,2)*r1 + 17*M*pow(r1,2) +
    	3*pow(r1,3))) + 2*pow(e,2)*M*(8*pow(M,3) + 36*pow(M,2)*r1 + 42*M*pow(r1,2) +
    	13*pow(r1,3))))))/(12.*pow(r1,6)*(pow(a,2) + r1*(-2*M + r1)));
    A042 =
    	3*pow(r1,3)*(pow(a,2)*(2*M + r1) + r1*(pow(l,2) + pow(r1,2)));
    A044 =
    	-(-4*pow(a,7)*e*l*M*(72*pow(M,2) + 59*M*r1 + 12*pow(r1,2)) +
    	4*a*e*l*M*pow(r1,3)*(pow(l,4)*(8*M - 6*r1) + 6*pow(l,2)*(6*M - 5*r1)*pow(r1,2) + 3*(10*M
    	- 9*r1)*pow(r1,4)) + pow(a,8)*(-13*pow(M,2)*r1 + 3*pow(r1,3) + 36*pow(e,2)*M*pow(2*M +
    	r1,2)) - 4*pow(a,5)*e*l*M*r1*(-48*pow(M,3) + 50*pow(M,2)*r1 + 140*M*pow(r1,2) +
    	49*pow(r1,3) + 6*pow(l,2)*(8*M + 3*r1)) + pow(r1,3)*(8*pow(l,6)*M*(-2*M + r1) +
    	pow(l,4)*pow(r1,2)*(-100*pow(M,2) + 4*(25 + pow(e,2))*M*r1 - 25*pow(r1,2)) +
    	2*pow(r1,6)*(-12*pow(M,2) + 4*(4 + 3*pow(e,2))*M*r1 - 5*pow(r1,2)) -
    	2*pow(l,2)*pow(r1,4)*(60*pow(M,2) - 4*(17 + 3*pow(e,2))*M*r1 + 19*pow(r1,2))) -
    	4*pow(a,3)*e*l*M*pow(r1,2)*(6*pow(l,4) + pow(l,2)*(-48*pow(M,2) + 36*M*r1 +
    	48*pow(r1,2)) + pow(r1,2)*(-66*pow(M,2) + 55*M*r1 + 64*pow(r1,2))) +
    	2*pow(a,2)*pow(r1,2)*(6*pow(l,6)*M + pow(l,2)*pow(r1,2)*(-96*pow(M,3) + 88*pow(M,2)*r1 +
    	6*(9 + 8*pow(e,2))*M*pow(r1,2) - 35*pow(r1,3)) + pow(r1,4)*(-12*pow(M,3) + 4*(10 +
    	9*pow(e,2))*pow(M,2)*r1 + 9*(1 + 6*pow(e,2))*M*pow(r1,2) - 15*pow(r1,3)) +
    	pow(l,4)*(-48*pow(M,3) - 8*(-3 + pow(e,2))*pow(M,2)*r1 + 2*(17 + 4*pow(e,2))*M*pow(r1,2)
    	- 11*pow(r1,3))) + pow(a,4)*r1*(3*pow(l,4)*(32*pow(M,2) + 4*(2 + pow(e,2))*M*r1 +
    	pow(r1,2)) + pow(r1,2)*(12*pow(M,4) + 116*pow(M,3)*r1 + (55 +
    	288*pow(e,2))*pow(M,2)*pow(r1,2) + 36*(-2 + 5*pow(e,2))*M*pow(r1,3) - 27*pow(r1,4)) -
    	2*pow(l,2)*(48*pow(M,4) + 4*(-7 + 12*pow(e,2))*pow(M,3)*r1 - 2*(53 +
    	24*pow(e,2))*pow(M,2)*pow(r1,2) + 4*(2 - 15*pow(e,2))*M*pow(r1,3) + 13*pow(r1,4))) +
    	2*pow(a,6)*(pow(l,2)*(72*pow(M,3) + 2*(23 + 24*pow(e,2))*pow(M,2)*r1 + 6*(1 +
    	4*pow(e,2))*M*pow(r1,2) + 3*pow(r1,3)) + r1*(-(r1*(-10*pow(M,3) + 37*pow(M,2)*r1 +
    	29*M*pow(r1,2) + 2*pow(r1,3))) + 6*pow(e,2)*M*(-8*pow(M,3) + 12*pow(M,2)*r1 +
    	30*M*pow(r1,2) + 11*pow(r1,3)))))/(24.*pow(r1,3)*(pow(a,2) + r1*(-2*M + r1)));
    A060 =
    	pow(r1,6);
    A062 = (6*pow(a,5)*e*l*M*(12*M + 5*r1) + 4*pow(a,3)*e*l*M*r1*(6*pow(l,2) -
    	20*pow(M,2) + 11*M*r1 + 18*pow(r1,2)) - 3*pow(a,6)*(r1*(M + r1) + 6*pow(e,2)*M*(2*M + r1))
    	+ 2*a*e*l*M*pow(r1,2)*(-4*pow(l,2)*(4*M - 3*r1) + pow(r1,2)*(-26*M + 21*r1)) +
    	pow(r1,2)*(8*pow(l,4)*M*(2*M - r1) + 4*pow(r1,4)*(2*pow(M,2) - (3 + 2*pow(e,2))*M*r1 +
    	pow(r1,2)) + pow(l,2)*pow(r1,2)*(36*pow(M,2) - 4*(8 + pow(e,2))*M*r1 + 7*pow(r1,2))) +
    	pow(a,2)*r1*(-12*pow(l,4)*M + 4*pow(l,2)*(10*pow(M,3) + (-3 + 4*pow(e,2))*pow(M,2)*r1 -
    	2*(3 + 2*pow(e,2))*M*pow(r1,2) + pow(r1,3)) + pow(r1,2)*(4*pow(M,3) + 4*(-4 +
    	pow(e,2))*pow(M,2)*r1 + (1 - 34*pow(e,2))*M*pow(r1,2) + 5*pow(r1,3))) -
    	pow(a,4)*(3*pow(l,2)*(12*pow(M,2) + 4*(1 + pow(e,2))*M*r1 + pow(r1,2)) +
    	2*r1*(r1*(-2*pow(M,2) - 5*M*r1 + pow(r1,2)) + 2*pow(e,2)*M*(-10*pow(M,2) + 8*M*r1 +
    	11*pow(r1,2)))))/(12.*(pow(a,2) + r1*(-2*M + r1)));
    A080 = (pow(r1,3)*(r1*(-3*pow(a,2) +
    	r1*(-2*M + r1)) + (4*M*(-3*pow(a,4)*pow(e,2) + 6*pow(a,3)*e*l + 2*a*e*l*r1*(-4*M +
    	3*r1) + pow(a,2)*(-3*pow(l,2) + 4*pow(e,2)*(M - r1)*r1) - r1*(pow(e,2)*pow(r1,3) +
    	pow(l,2)*(-4*M + 2*r1))))/(pow(a,2) + r1*(-2*M + r1))))/24.;
    A106 = -(pow(pow(l,2) +
    	pow(r1,2) + (pow(a,2)*(2*M + r1))/r1,2)*(-(pow(a,2)*M) + pow(r1,3) + (2*l*(pow(a,3)*e*M -
    	pow(a,2)*l*M + 3*a*e*M*pow(r1,2) + l*pow(r1,2)*(-2*M + r1)))/(pow(a,2) + r1*(-2*M +
    	r1))))/(2.*pow(r1,2));
    A108 = ((pow(a,2)*(2*M + r1) + r1*(pow(l,2) +
    	pow(r1,2)))*(-2*pow(a,11)*e*l*pow(M,2)*(80*pow(M,2) + 73*M*r1 + 12*pow(r1,2)) -
    	2*pow(a,9)*e*l*M*r1*(32*(-1 + 6*pow(e,2))*pow(M,4) + 12*(23 + 24*pow(e,2))*pow(M,3)*r1 +
    	(353 + 144*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(11 + 12*pow(e,2))*M*pow(r1,3) - 12*pow(r1,4)
    	+ 4*pow(l,2)*M*(16*M + 3*r1)) - 2*a*e*l*M*pow(r1,5)*(pow(r1,6)*(-68*pow(M,2) + 100*M*r1 -
    	33*pow(r1,2)) + 6*pow(l,4)*pow(r1,2)*(20*pow(M,2) + 4*pow(e,2)*M*r1 - 5*pow(r1,2)) +
    	2*pow(l,2)*pow(r1,4)*(68*pow(M,2) + 2*(-14 + 15*pow(e,2))*M*r1 - 3*pow(r1,2)) +
    	12*pow(l,6)*(4*pow(M,2) - pow(r1,2))) + pow(a,12)*pow(M,2)*(-3*r1*(3*M + 2*r1) +
    	8*pow(e,2)*(10*pow(M,2) + 11*M*r1 + 3*pow(r1,2))) -
    	2*pow(a,3)*e*l*M*pow(r1,3)*(12*pow(l,6)*(12*pow(M,2) - 2*M*r1 - pow(r1,2)) +
    	pow(r1,5)*(-768*pow(M,3) + 360*pow(M,2)*r1 + 6*(29 + 4*pow(e,2))*M*pow(r1,2) -
    	97*pow(r1,3)) + 12*pow(l,4)*r1*(20*pow(M,3) + 4*(10 + pow(e,2))*pow(M,2)*r1 + 2*(-5 +
    	2*pow(e,2))*M*pow(r1,2) - 5*pow(r1,3)) - 2*pow(l,2)*pow(r1,3)*(340*pow(M,3) - 6*(61 +
    	20*pow(e,2))*pow(M,2)*r1 + 10*(8 - 9*pow(e,2))*M*pow(r1,2) + 29*pow(r1,3))) -
    	2*pow(a,5)*e*l*M*pow(r1,2)*(6*pow(l,4)*(2*M + r1)*(60*pow(M,2) + 2*(-5 +
    	2*pow(e,2))*M*r1 - 5*pow(r1,2)) + pow(r1,3)*(-420*pow(M,4) - 876*pow(M,3)*r1 + (919 +
    	144*pow(e,2))*pow(M,2)*pow(r1,2) + 6*(17 + 12*pow(e,2))*M*pow(r1,3) - 107*pow(r1,4)) +
    	2*pow(l,2)*r1*(-304*pow(M,4) + 4*(13 + 30*pow(e,2))*pow(M,3)*r1 + 24*(17 +
    	10*pow(e,2))*pow(M,2)*pow(r1,2) + 18*(-7 + 5*pow(e,2))*M*pow(r1,3) - 41*pow(r1,4))) +
    	pow(a,2)*pow(r1,3)*(48*pow(l,8)*pow(M,2)*(2*M - r1) + 24*pow(l,6)*M*pow(r1,2)*(2*(5 +
    	4*pow(e,2))*pow(M,2) + (-5 + 4*pow(e,2))*M*r1 - 2*pow(e,2)*pow(r1,2)) -
    	4*pow(l,4)*M*pow(r1,3)*(116*pow(M,3) - 18*(9 + 10*pow(e,2))*pow(M,2)*r1 + 4*(17 -
    	15*pow(e,2))*M*pow(r1,2) + (-8 + 45*pow(e,2))*pow(r1,3)) +
    	8*pow(l,2)*pow(r1,5)*(-94*pow(M,4) + (103 + 33*pow(e,2))*pow(M,3)*r1 + (-28 +
    	11*pow(e,2))*pow(M,2)*pow(r1,2) - 6*(1 + 3*pow(e,2))*M*pow(r1,3) + 3*pow(r1,4)) +
    	pow(r1,7)*(-120*pow(M,4) - 4*(-47 + 50*pow(e,2))*pow(M,3)*r1 + 2*(-51 +
    	86*pow(e,2))*pow(M,2)*pow(r1,2) + (1 - 46*pow(e,2))*M*pow(r1,3) + 9*pow(r1,4))) -
    	2*pow(a,7)*e*l*M*r1*(2*pow(l,2)*(288*pow(M,4) + 8*(29 + 15*pow(e,2))*pow(M,3)*r1 + 2*(71
    	+ 60*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-26 + 15*pow(e,2))*M*pow(r1,3) - 15*pow(r1,4)) +
    	r1*(-256*pow(M,5) - 556*pow(M,4)*r1 + 8*(-7 + 36*pow(e,2))*pow(M,3)*pow(r1,2) + (827 +
    	288*pow(e,2))*pow(M,2)*pow(r1,3) + 2*(19 + 36*pow(e,2))*M*pow(r1,4) - 55*pow(r1,5))) +
    	pow(a,4)*pow(r1,2)*(24*pow(l,6)*M*(20*pow(M,3) + 12*pow(e,2)*pow(M,2)*r1 + (-5 +
    	2*pow(e,2))*M*pow(r1,2) - pow(e,2)*pow(r1,3)) - 4*pow(l,4)*M*r1*(152*pow(M,4) - 4*(11 +
    	60*pow(e,2))*pow(M,3)*r1 - 30*(5 + 14*pow(e,2))*pow(M,2)*pow(r1,2) - 30*(-3 +
    	pow(e,2))*M*pow(r1,3) + (-4 + 45*pow(e,2))*pow(r1,4)) +
    	2*pow(l,2)*pow(r1,3)*(-148*pow(M,5) - 8*(50 + 57*pow(e,2))*pow(M,4)*r1 + (639 +
    	452*pow(e,2))*pow(M,3)*pow(r1,2) + (-221 + 12*pow(e,2))*pow(M,2)*pow(r1,3) -
    	114*pow(e,2)*M*pow(r1,4) + 6*pow(r1,5)) + pow(r1,5)*(48*pow(M,5) - 4*(37 +
    	168*pow(e,2))*pow(M,4)*r1 - 4*(-71 + 42*pow(e,2))*pow(M,3)*pow(r1,2) + (-205 +
    	448*pow(e,2))*pow(M,2)*pow(r1,3) + 3*(13 - 28*pow(e,2))*M*pow(r1,4) + 9*pow(r1,5))) +
    	pow(a,6)*r1*(4*pow(l,4)*M*(96*pow(M,4) + 4*(11 + 90*pow(e,2))*pow(M,3)*r1 + 16*(2 +
    	15*pow(e,2))*pow(M,2)*pow(r1,2) - 38*M*pow(r1,3) - 15*pow(e,2)*pow(r1,4)) -
    	2*pow(l,2)*M*r1*(128*pow(M,5) + 4*(47 + 76*pow(e,2))*pow(M,4)*r1 - 4*(-3 +
    	76*pow(e,2))*pow(M,3)*pow(r1,2) - (415 + 716*pow(e,2))*pow(M,2)*pow(r1,3) + (149 -
    	8*pow(e,2))*M*pow(r1,4) + 4*(-1 + 20*pow(e,2))*pow(r1,5)) + pow(r1,3)*(56*pow(M,6) - 4*(9
    	+ 128*pow(e,2))*pow(M,5)*r1 - 2*(41 + 520*pow(e,2))*pow(M,4)*pow(r1,2) + (181 +
    	464*pow(e,2))*pow(M,3)*pow(r1,3) + (-163 + 568*pow(e,2))*pow(M,2)*pow(r1,4) + (49 -
    	76*pow(e,2))*M*pow(r1,5) + 3*pow(r1,6))) + (2*M -
    	r1)*pow(r1,8)*(2*pow(e,2)*M*(12*pow(l,6) + 30*pow(l,4)*pow(r1,2) + 17*pow(l,2)*pow(r1,4)
    	+ 5*pow(r1,6)) + (2*M - r1)*r1*(16*pow(l,4)*M + pow(r1,4)*(4*M + 3*r1) +
    	2*pow(l,2)*pow(r1,2)*(M + 6*r1))) + pow(a,10)*M*(2*pow(l,2)*M*(40*pow(M,2) + (29 +
    	32*pow(e,2))*M*r1 + 3*(-1 + 4*pow(e,2))*pow(r1,2)) + r1*(r1*(50*pow(M,3) + 11*pow(M,2)*r1
    	- 19*M*pow(r1,2) + 6*pow(r1,3)) + pow(e,2)*(-32*pow(M,4) + 336*pow(M,3)*r1 +
    	504*pow(M,2)*pow(r1,2) + 156*M*pow(r1,3) - 6*pow(r1,4)))) +
    	pow(a,8)*M*r1*(64*pow(l,4)*pow(M,2) + 2*pow(l,2)*(16*(-1 + 36*pow(e,2))*pow(M,4) + 4*(27
    	+ 166*pow(e,2))*pow(M,3)*r1 + (141 + 364*pow(e,2))*pow(M,2)*pow(r1,2) + 3*(-13 +
    	6*pow(e,2))*M*pow(r1,3) + 3*(1 - 7*pow(e,2))*pow(r1,4)) + r1*(r1*(-92*pow(M,4) +
    	20*pow(M,3)*r1 + 65*pow(M,2)*pow(r1,2) - 69*M*pow(r1,3) + 25*pow(r1,4)) -
    	2*pow(e,2)*(128*pow(M,5) + 368*pow(M,4)*r1 + 120*pow(M,3)*pow(r1,2) -
    	424*pow(M,2)*pow(r1,3) - 202*M*pow(r1,4) + 17*pow(r1,5))))))/(48.*pow(r1,10)*pow(pow(a,2) +
    	r1*(-2*M + r1),2));
    A124 = -((pow(a,2)*(2*M + r1) + r1*(pow(l,2) +
    	pow(r1,2)))*(4*pow(a,3)*e*l*M + pow(a,4)*r1 + 12*a*e*l*M*pow(r1,2) - (2*M -
    	r1)*pow(r1,2)*(5*pow(l,2) + 3*pow(r1,2)) + pow(a,2)*(-2*(M - 2*r1)*pow(r1,2) +
    	pow(l,2)*(-4*M + r1))))/(2.*r1*(pow(a,2) + r1*(-2*M + r1)));
    A126 =
    	(2*pow(a,11)*e*l*pow(M,2)*(16*(-16 + 9*pow(e,2))*pow(M,2) + (-293 + 144*pow(e,2))*M*r1 +
    	6*(-13 + 6*pow(e,2))*pow(r1,2)) - 2*a*e*l*M*pow(r1,5)*(pow(r1,6)*(-236*pow(M,2) +
    	232*M*r1 - 57*pow(r1,2)) + 2*pow(l,4)*pow(r1,2)*(328*pow(M,2) + 2*(-79 + 21*pow(e,2))*M*r1
    	- 3*pow(r1,2)) + 42*pow(l,6)*(4*pow(M,2) - pow(r1,2)) + 24*pow(l,2)*pow(r1,4)*(16*pow(M,2)
    	+ 5*(-2 + pow(e,2))*M*r1 + pow(r1,2))) + (2*M - r1)*pow(r1,7)*(8*M*pow(r1,6)*(4*M + (-2 +
    	5*pow(e,2))*r1) + pow(l,2)*pow(r1,4)*(20*pow(M,2) + 4*(-4 + 33*pow(e,2))*M*r1 +
    	3*pow(r1,2)) + 4*pow(l,6)*(16*pow(M,2) + (-20 + 21*pow(e,2))*M*r1 + 6*pow(r1,2)) +
    	4*pow(l,4)*pow(r1,2)*(22*pow(M,2) + (-29 + 47*pow(e,2))*M*r1 + 9*pow(r1,2))) -
    	2*pow(a,3)*e*l*M*pow(r1,3)*(42*pow(l,6)*(12*pow(M,2) - 2*M*r1 - pow(r1,2)) +
    	pow(r1,5)*(-2208*pow(M,3) + 4*(155 + 54*pow(e,2))*pow(M,2)*r1 + 4*(154 -
    	9*pow(e,2))*M*pow(r1,2) - 185*pow(r1,3)) + 4*pow(l,4)*r1*(256*pow(M,3) + 2*(127 +
    	21*pow(e,2))*pow(M,2)*r1 + 6*(-24 + 7*pow(e,2))*M*pow(r1,2) - 21*pow(r1,3)) +
    	2*pow(l,2)*pow(r1,3)*(-204*pow(M,3) + 4*(92 + 87*pow(e,2))*pow(M,2)*r1 + (-173 +
    	150*pow(e,2))*M*pow(r1,2) - 16*pow(r1,3))) + pow(a,12)*M*(-3*r1*(21*pow(M,2) + 23*M*r1 +
    	6*pow(r1,2)) + 2*pow(e,2)*(128*pow(M,3) + 172*pow(M,2)*r1 + 72*M*pow(r1,2) + 9*pow(r1,3)))
    	+ 2*pow(a,9)*e*l*M*(2*pow(l,2)*M*(216*pow(M,2) + 4*(2 + 9*pow(e,2))*M*r1 + 3*(-17 +
    	6*pow(e,2))*pow(r1,2)) + r1*(-288*(1 + 3*pow(e,2))*pow(M,4) - 4*(275 +
    	144*pow(e,2))*pow(M,3)*r1 + 3*(-403 + 72*pow(e,2))*pow(M,2)*pow(r1,2) + 6*(-43 +
    	24*pow(e,2))*M*pow(r1,3) + 36*pow(r1,4))) + pow(a,2)*pow(r1,3)*(168*pow(l,8)*pow(M,2)*(2*M
    	- r1) + pow(l,2)*pow(r1,5)*(-2112*pow(M,4) + 8*(223 + 162*pow(e,2))*pow(M,3)*r1 + 4*(-91 +
    	4*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(15 - 232*pow(e,2))*M*pow(r1,3) - 15*pow(r1,4)) +
    	2*pow(l,6)*r1*(56*pow(M,4) + 12*(23 + 28*pow(e,2))*pow(M,3)*r1 + 14*(-13 +
    	12*pow(e,2))*pow(M,2)*pow(r1,2) + 3*(13 - 28*pow(e,2))*M*pow(r1,3) - 12*pow(r1,4)) +
    	pow(r1,7)*(-264*pow(M,4) + (444 - 512*pow(e,2))*pow(M,3)*r1 + 2*(-149 +
    	246*pow(e,2))*pow(M,2)*pow(r1,2) + (65 - 158*pow(e,2))*M*pow(r1,3) + 3*pow(r1,4)) -
    	2*pow(l,4)*pow(r1,3)*(592*pow(M,4) - 4*(139 + 329*pow(e,2))*pow(M,3)*r1 + 4*(47 -
    	11*pow(e,2))*pow(M,2)*pow(r1,2) + (-101 + 265*pow(e,2))*M*pow(r1,3) + 36*pow(r1,4))) -
    	2*pow(a,5)*e*l*M*pow(r1,2)*(2*pow(l,4)*(1224*pow(M,3) + 28*(-1 + 3*pow(e,2))*pow(M,2)*r1
    	+ 2*(-95 + 21*pow(e,2))*M*pow(r1,2) - 39*pow(r1,3)) + pow(r1,3)*(-3300*pow(M,4) + 4*(-517 +
    	216*pow(e,2))*pow(M,3)*r1 + (2423 + 504*pow(e,2))*pow(M,2)*pow(r1,2) + 6*(121 -
    	24*pow(e,2))*M*pow(r1,3) - 235*pow(r1,4)) - 2*pow(l,2)*r1*(96*pow(M,4) - 4*(169 +
    	114*pow(e,2))*pow(M,3)*r1 - 36*(4 + 17*pow(e,2))*pow(M,2)*pow(r1,2) + (131 -
    	102*pow(e,2))*M*pow(r1,3) + 64*pow(r1,4))) - 2*pow(a,7)*e*l*M*r1*(-24*pow(l,4)*M*(9*M -
    	r1) + 2*pow(l,2)*(1296*pow(M,4) + 8*(35 + 51*pow(e,2))*pow(M,3)*r1 + 4*(-52 +
    	57*pow(e,2))*pow(M,2)*pow(r1,2) - (19 + 6*pow(e,2))*M*pow(r1,3) - 36*pow(r1,4)) +
    	r1*(-1600*pow(M,5) + 4*(-515 + 216*pow(e,2))*pow(M,4)*r1 + 72*(13 +
    	22*pow(e,2))*pow(M,3)*pow(r1,2) + 3*(857 + 72*pow(e,2))*pow(M,2)*pow(r1,3) + 18*(29 -
    	12*pow(e,2))*M*pow(r1,4) - 143*pow(r1,5))) +
    	pow(a,4)*pow(r1,2)*(4*pow(l,6)*M*(408*pow(M,3) + 4*(-32 + 63*pow(e,2))*pow(M,2)*r1 +
    	2*(-26 + 21*pow(e,2))*M*pow(r1,2) - (5 + 21*pow(e,2))*pow(r1,3)) +
    	pow(r1,5)*(-384*pow(M,5) - 24*(-9 + 106*pow(e,2))*pow(M,4)*r1 + 12*(35 -
    	22*pow(e,2))*pow(M,3)*pow(r1,2) + 2*(-213 + 608*pow(e,2))*pow(M,2)*pow(r1,3) + (103 -
    	214*pow(e,2))*M*pow(r1,4) + 9*pow(r1,5)) - 2*pow(l,4)*r1*(528*pow(M,5) - 40*(12 +
    	47*pow(e,2))*pow(M,4)*r1 - 4*(60 + 491*pow(e,2))*pow(M,3)*pow(r1,2) + 6*(12 +
    	65*pow(e,2))*pow(M,2)*pow(r1,3) + (35 + 239*pow(e,2))*M*pow(r1,4) + 18*pow(r1,5)) -
    	pow(l,2)*pow(r1,3)*(3096*pow(M,5) + 12*(55 - 192*pow(e,2))*pow(M,4)*r1 - 2*(1355 +
    	788*pow(e,2))*pow(M,3)*pow(r1,2) + (407 + 1116*pow(e,2))*pow(M,2)*pow(r1,3) + 12*(9 +
    	47*pow(e,2))*M*pow(r1,4) + 21*pow(r1,5))) - pow(a,6)*r1*(6*pow(l,6)*M*(24*pow(M,2) -
    	10*M*r1 + 3*pow(r1,2)) - 2*pow(l,4)*M*(864*pow(M,4) + 8*(-1 + 306*pow(e,2))*pow(M,3)*r1 +
    	4*(-29 + 164*pow(e,2))*pow(M,2)*pow(r1,2) - 32*(-2 + 9*pow(e,2))*M*pow(r1,3) - (69 +
    	59*pow(e,2))*pow(r1,4)) - pow(r1,3)*(152*pow(M,6) - 4*(23 + 800*pow(e,2))*pow(M,5)*r1 +
    	2*(31 - 1672*pow(e,2))*pow(M,4)*pow(r1,2) + (163 + 1792*pow(e,2))*pow(M,3)*pow(r1,3) +
    	(-181 + 1704*pow(e,2))*pow(M,2)*pow(r1,4) + (43 - 76*pow(e,2))*M*pow(r1,5) + 9*pow(r1,6))
    	+ pow(l,2)*r1*(1600*pow(M,6) - 8*(-167 + 300*pow(e,2))*pow(M,5)*r1 - 4*(243 +
    	1228*pow(e,2))*pow(M,4)*pow(r1,2) + 6*(-333 + 40*pow(e,2))*pow(M,3)*pow(r1,3) + (47 +
    	1196*pow(e,2))*pow(M,2)*pow(r1,4) + 4*(62 + 57*pow(e,2))*M*pow(r1,5) + 9*pow(r1,6))) +
    	pow(a,10)*M*(pow(l,2)*(-32*(-8 + 27*pow(e,2))*pow(M,3) + (242 - 448*pow(e,2))*pow(M,2)*r1
    	+ 3*(-19 + 20*pow(e,2))*M*pow(r1,2) + 18*(-3 + 2*pow(e,2))*pow(r1,3)) +
    	r1*(r1*(290*pow(M,3) + 165*pow(M,2)*r1 - 141*M*pow(r1,2) - 68*pow(r1,3)) +
    	2*pow(e,2)*(144*pow(M,4) + 680*pow(M,3)*r1 + 864*pow(M,2)*pow(r1,2) + 358*M*pow(r1,3) +
    	37*pow(r1,4)))) + pow(a,8)*(-2*pow(l,4)*M*(144*pow(M,3) + 8*(-8 + 27*pow(e,2))*pow(M,2)*r1
    	+ 6*(-6 + 7*pow(e,2))*M*pow(r1,2) - 9*(-3 + pow(e,2))*pow(r1,3)) + pow(l,2)*M*r1*(288*(1
    	+ 18*pow(e,2))*pow(M,4) + 8*(105 + 286*pow(e,2))*pow(M,3)*r1 - 2*(-589 +
    	516*pow(e,2))*pow(M,2)*pow(r1,2) + (15 - 268*pow(e,2))*M*pow(r1,3) + 2*(-93 +
    	20*pow(e,2))*pow(r1,4)) + pow(r1,2)*(-4*pow(e,2)*M*(400*pow(M,5) + 696*pow(M,4)*r1 -
    	12*pow(M,3)*pow(r1,2) - 732*pow(M,2)*pow(r1,3) - 368*M*pow(r1,4) - 19*pow(r1,5)) +
    	r1*(-404*pow(M,5) - 8*pow(M,4)*r1 + 375*pow(M,3)*pow(r1,2) - 61*pow(M,2)*pow(r1,3) -
    	61*M*pow(r1,4) + 3*pow(r1,5)))))/(48.*pow(r1,7)*pow(pow(a,2) + r1*(-2*M + r1),2));
    A142 =
    	(pow(r1,2)*(pow(a,2)*M - pow(r1,3) - (2*l*(pow(a,3)*e*M - pow(a,2)*l*M +
    	3*a*e*M*pow(r1,2) + l*pow(r1,2)*(-2*M + r1)))/(pow(a,2) + r1*(-2*M + r1)) -
    	2*(pow(a,2)*(2*M + r1) + r1*(pow(l,2) + pow(r1,2)))))/2.;
    A144 =
    	(2*pow(a,9)*e*l*M*(16*(-16 + 9*pow(e,2))*pow(M,2) + 12*(-13 + 6*pow(e,2))*M*r1 -
    	3*pow(r1,2)) + (2*M - r1)*pow(r1,4)*(16*pow(l,6)*M*(2*M - r1) + 2*pow(r1,6)*(24*pow(M,2) +
    	2*(-11 + 15*pow(e,2))*M*r1 + 5*pow(r1,2)) + 2*pow(l,4)*pow(r1,2)*(34*pow(M,2) + (-33 +
    	47*pow(e,2))*M*r1 + 8*pow(r1,2)) + pow(l,2)*pow(r1,4)*(-12*pow(M,2) + 4*(-4 +
    	33*pow(e,2))*M*r1 + 11*pow(r1,2))) - 2*a*e*l*M*pow(r1,4)*(2*pow(l,4)*(200*pow(M,2) -
    	94*M*r1 - 3*pow(r1,2)) + 3*pow(l,2)*pow(r1,2)*(188*pow(M,2) + 4*(-27 + 5*pow(e,2))*M*r1 +
    	7*pow(r1,2)) - 3*pow(r1,4)*(100*pow(M,2) - 104*M*r1 + 27*pow(r1,2))) +
    	pow(a,10)*M*(-9*r1*(13*M + 7*r1) + pow(e,2)*(256*pow(M,2) + 264*M*r1 + 66*pow(r1,2))) -
    	2*pow(a,3)*e*l*M*pow(r1,2)*(4*pow(l,4)*(180*pow(M,2) - 92*M*r1 - 3*pow(r1,2)) +
    	pow(l,2)*r1*(768*pow(M,3) + 4*(67 + 84*pow(e,2))*pow(M,2)*r1 + 4*(-149 +
    	9*pow(e,2))*M*pow(r1,2) + 15*pow(r1,3)) - 4*pow(r1,3)*(426*pow(M,3) - (85 +
    	108*pow(e,2))*pow(M,2)*r1 + 9*(-17 + 4*pow(e,2))*M*pow(r1,2) + 46*pow(r1,3))) -
    	2*pow(a,5)*e*l*M*r1*(-6*pow(l,4)*(18*M + r1) + pow(l,2)*(1728*pow(M,3) + 24*(-33 +
    	10*pow(e,2))*pow(M,2)*r1 - 4*(83 + 15*pow(e,2))*M*pow(r1,2) - 9*pow(r1,3)) +
    	2*r1*(-672*pow(M,4) + 4*(-107 + 108*pow(e,2))*pow(M,3)*r1 + 6*(79 +
    	36*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(161 - 90*pow(e,2))*M*pow(r1,3) - 65*pow(r1,4))) -
    	pow(a,2)*pow(r1,2)*(-2*pow(l,6)*M*(240*pow(M,2) - 154*M*r1 + 17*pow(r1,2)) +
    	2*pow(l,4)*r1*(48*pow(M,4) - 4*(35 + 176*pow(e,2))*pow(M,3)*r1 + 8*(8 +
    	15*pow(e,2))*pow(M,2)*pow(r1,2) + (-35 + 73*pow(e,2))*M*pow(r1,3) + 16*pow(r1,4)) +
    	pow(r1,5)*(264*pow(M,4) + 4*(-35 + 144*pow(e,2))*pow(M,3)*r1 + 2*(41 -
    	222*pow(e,2))*pow(M,2)*pow(r1,2) + (-101 + 138*pow(e,2))*M*pow(r1,3) + 31*pow(r1,4)) +
    	pow(l,2)*pow(r1,3)*(1656*pow(M,4) - 12*(65 + 166*pow(e,2))*pow(M,3)*r1 + 2*(-149 +
    	374*pow(e,2))*pow(M,2)*pow(r1,2) + (47 + 256*pow(e,2))*M*pow(r1,3) + 45*pow(r1,4))) -
    	pow(a,4)*r1*(6*pow(l,6)*M*(12*M + r1) + 2*pow(l,4)*(-576*pow(M,4) - 36*(-9 +
    	20*pow(e,2))*pow(M,3)*r1 + 2*(-19 + 137*pow(e,2))*pow(M,2)*pow(r1,2) + (47 +
    	29*pow(e,2))*M*pow(r1,3) + 8*pow(r1,4)) + pow(r1,3)*(192*pow(M,5) + 4*(85 +
    	468*pow(e,2))*pow(M,4)*r1 + 228*(-1 + 2*pow(e,2))*pow(M,3)*pow(r1,2) - 3*(117 +
    	268*pow(e,2))*pow(M,2)*pow(r1,3) - 12*(-8 + pow(e,2))*M*pow(r1,4) + 39*pow(r1,5)) +
    	pow(l,2)*r1*(1344*pow(M,5) + 32*(13 - 105*pow(e,2))*pow(M,4)*r1 - 4*(243 +
    	280*pow(e,2))*pow(M,3)*pow(r1,2) + 28*(-29 + 75*pow(e,2))*pow(M,2)*pow(r1,3) + (379 +
    	56*pow(e,2))*M*pow(r1,4) + 57*pow(r1,5))) + 2*pow(a,7)*e*l*M*(4*r1*(-8*(5 +
    	18*pow(e,2))*pow(M,3) + 3*(-59 + 12*pow(e,2))*pow(M,2)*r1 + (-121 +
    	72*pow(e,2))*M*pow(r1,2) + 6*pow(r1,3)) + 3*pow(l,2)*(144*pow(M,2) + pow(r1,2) + 4*M*(r1 +
    	3*pow(e,2)*r1))) + pow(a,8)*(pow(l,2)*M*(-32*(-8 + 27*pow(e,2))*pow(M,2) + 12*(4 -
    	19*pow(e,2))*M*r1 + 3*(-41 + 20*pow(e,2))*pow(r1,2)) + r1*(-(r1*(-420*pow(M,3) +
    	115*pow(M,2)*r1 + 256*M*pow(r1,2) + 7*pow(r1,3))) + 4*pow(e,2)*M*(40*pow(M,3) +
    	178*pow(M,2)*r1 + 209*M*pow(r1,2) + 60*pow(r1,3)))) - pow(a,6)*(6*pow(l,4)*M*(48*pow(M,2) +
    	2*(-5 + 18*pow(e,2))*M*r1 + (11 + pow(e,2))*pow(r1,2)) + pow(l,2)*r1*(-32*(5 +
    	108*pow(e,2))*pow(M,4) + 8*(-88 + 153*pow(e,2))*pow(M,3)*r1 + 2*(-247 +
    	658*pow(e,2))*pow(M,2)*pow(r1,2) + (417 - 128*pow(e,2))*M*pow(r1,3) + 23*pow(r1,4)) +
    	pow(r1,2)*(4*pow(e,2)*M*(336*pow(M,4) + 324*pow(M,3)*r1 - 128*pow(M,2)*pow(r1,2) -
    	263*M*pow(r1,3) - 66*pow(r1,4)) + r1*(276*pow(M,4) - 652*pow(M,3)*r1 -
    	323*pow(M,2)*pow(r1,2) + 326*M*pow(r1,3) + 25*pow(r1,4)))))/(48.*pow(r1,4)*pow(pow(a,2) +
    	r1*(-2*M + r1),2));
    A160 = -pow(r1,5)/2.;
    A162 = (6*pow(a,7)*e*l*M*(12*(-4 + pow(e,2))*M
    	- 11*r1) + 3*pow(a,8)*(3*r1*(-3*M + r1) + 2*pow(e,2)*M*(24*M + 13*r1)) + 2*a*e*l*M*(2*M -
    	r1)*pow(r1,3)*((82*M - 51*r1)*pow(r1,2) + pow(l,2)*(-124*M + 30*r1)) + (2*M -
    	r1)*pow(r1,3)*(8*pow(l,4)*M*(2*M - r1) + pow(l,2)*pow(r1,2)*(-36*pow(M,2) + 4*(8 +
    	11*pow(e,2))*M*r1 - 7*pow(r1,2)) + 8*pow(r1,4)*(4*pow(M,2) + (-4 + 5*pow(e,2))*M*r1 +
    	pow(r1,2))) - 2*pow(a,5)*e*l*M*(-12*pow(l,2)*(9*M - 2*r1) + r1*(8*(-11 +
    	15*pow(e,2))*pow(M,2) - 24*(-9 + 5*pow(e,2))*M*r1 + 63*pow(r1,2))) -
    	2*pow(a,3)*e*l*M*r1*(pow(l,2)*(360*pow(M,2) - 352*M*r1 + 54*pow(r1,2)) - r1*(400*pow(M,3)
    	+ 4*(19 - 54*pow(e,2))*pow(M,2)*r1 + 4*(-40 + 21*pow(e,2))*M*pow(r1,2) + 21*pow(r1,3))) -
    	pow(a,2)*r1*(-8*pow(l,4)*M*(30*pow(M,2) - 29*M*r1 + 7*pow(r1,2)) + pow(r1,3)*(88*pow(M,4) +
    	4*(27 + 64*pow(e,2))*pow(M,3)*r1 - 2*(57 + 50*pow(e,2))*pow(M,2)*pow(r1,2) + (-11 +
    	26*pow(e,2))*M*pow(r1,3) + 15*pow(r1,4)) + pow(l,2)*r1*(400*pow(M,4) + 8*(-26 +
    	63*pow(e,2))*pow(M,2)*pow(r1,2) + 4*(17 + 4*pow(e,2))*M*pow(r1,3) + pow(r1,4) -
    	64*pow(M,3)*(r1 + 14*pow(e,2)*r1))) + pow(a,6)*(-3*pow(l,2)*(24*(-2 +
    	3*pow(e,2))*pow(M,2) + (4 - 8*pow(e,2))*M*r1 - 3*pow(r1,2)) + r1*(r1*(86*pow(M,2) -
    	139*M*r1 + 19*pow(r1,2)) + 2*pow(e,2)*M*(-44*pow(M,2) + 98*M*r1 + 105*pow(r1,2)))) +
    	pow(a,4)*(24*pow(l,4)*M*(-3*M + r1) + pow(l,2)*r1*(8*(-11 + 90*pow(e,2))*pow(M,3) + 4*(59
    	- 178*pow(e,2))*pow(M,2)*r1 + 2*(-41 + 26*pow(e,2))*M*pow(r1,2) + pow(r1,3)) +
    	pow(r1,2)*(r1*(-20*pow(M,3) + 296*pow(M,2)*r1 - 149*M*pow(r1,2) + 3*pow(r1,3)) +
    	2*pow(e,2)*M*(-200*pow(M,3) - 108*pow(M,2)*r1 + 36*M*pow(r1,2) +
    	73*pow(r1,3)))))/(48.*r1*pow(pow(a,2) + r1*(-2*M + r1),2));
    A180 =
    	-(pow(r1,2)*(60*pow(a,3)*e*l*M + 4*a*e*l*M*r1*(8*M + 3*r1) - 3*pow(a,4)*(10*pow(e,2)*M
    	+ 3*r1) + r1*(8*pow(l,2)*M*(-2*M + r1) + pow(r1,2)*(8*pow(M,2) + 2*(-3 + 5*pow(e,2))*M*r1
    	+ pow(r1,2))) - 2*pow(a,2)*(15*pow(l,2)*M + r1*(r1*(-7*M + 4*r1) + 2*pow(e,2)*M*(4*M +
    	5*r1)))))/(48.*(pow(a,2) + r1*(-2*M + r1)));
    A204 = (3*pow(pow(a,2)*(2*M + r1) +
    	r1*(pow(l,2) + pow(r1,2)),2))/(pow(a,2) + r1*(-2*M + r1));
    A206 = ((pow(a,2)*(2*M + r1) +
    	r1*(pow(l,2) + pow(r1,2)))*(-4*pow(a,7)*e*l*M*(12*pow(M,2) + 23*M*r1 + 3*pow(r1,2)) -
    	4*a*e*l*M*pow(r1,5)*(pow(l,2)*(67*M - 36*r1) + (17*M - 15*r1)*pow(r1,2)) -
    	4*pow(a,5)*e*l*M*r1*(-64*pow(M,3) - 104*pow(M,2)*r1 + 55*M*pow(r1,2) + 7*pow(r1,3) +
    	3*pow(l,2)*(8*M + r1)) + pow(a,8)*M*(-(r1*(M + 6*r1)) + 12*pow(e,2)*(2*pow(M,2) + 3*M*r1 +
    	pow(r1,2))) + 4*pow(a,3)*e*l*M*pow(r1,2)*(pow(l,2)*(16*pow(M,2) - 73*M*r1 + 3*pow(r1,2))
    	+ pow(r1,2)*(196*pow(M,2) - 17*M*r1 + 11*pow(r1,2))) + pow(r1,5)*(8*pow(l,4)*(10*pow(M,2) -
    	11*M*r1 + 3*pow(r1,2)) + pow(r1,4)*(40*pow(M,2) - 2*(13 + 2*pow(e,2))*M*r1 + 3*pow(r1,2)) +
    	2*pow(l,2)*pow(r1,2)*(54*pow(M,2) + (-51 + 2*pow(e,2))*M*r1 + 12*pow(r1,2))) -
    	2*pow(a,2)*pow(r1,2)*(M*pow(r1,4)*(58*pow(M,2) + 3*(-27 + 26*pow(e,2))*M*r1 + 33*pow(r1,2))
    	+ 2*pow(l,4)*(8*pow(M,3) - 33*pow(M,2)*r1 + 10*M*pow(r1,2) + 3*pow(r1,3)) -
    	pow(l,2)*pow(r1,2)*(-194*pow(M,3) + (171 + 104*pow(e,2))*pow(M,2)*r1 + (-57 +
    	10*pow(e,2))*M*pow(r1,2) + 3*pow(r1,3))) + pow(a,4)*r1*(48*pow(l,4)*pow(M,2) +
    	pow(r1,2)*(4*pow(M,4) + 20*(1 - 18*pow(e,2))*pow(M,3)*r1 - 3*(-53 +
    	76*pow(e,2))*pow(M,2)*pow(r1,2) + 24*(-4 + pow(e,2))*M*pow(r1,3) - 9*pow(r1,4)) -
    	2*pow(l,2)*(64*pow(M,4) + 8*(11 + 2*pow(e,2))*pow(M,3)*r1 - (109 +
    	80*pow(e,2))*pow(M,2)*pow(r1,2) + (51 - 14*pow(e,2))*M*pow(r1,3) + 9*pow(r1,4))) +
    	2*pow(a,6)*(pow(l,2)*M*(12*pow(M,2) + 4*(7 + 6*pow(e,2))*M*r1 + 3*(-1 +
    	2*pow(e,2))*pow(r1,2)) - r1*(pow(r1,3)*(31*M + 3*r1) + 2*pow(e,2)*M*(32*pow(M,3) +
    	60*pow(M,2)*r1 + 9*M*pow(r1,2) - 8*pow(r1,3))))))/(24.*pow(r1,6)*pow(pow(a,2) + r1*(-2*M +
    	r1),2));
    A222 = (6*pow(r1,3)*(pow(a,2)*(2*M + r1) + r1*(pow(l,2) + pow(r1,2))))/(pow(a,2) +
    	r1*(-2*M + r1));
    A224 = (4*pow(a,7)*e*l*M*(24*pow(M,2) + 19*M*r1 + 9*pow(r1,2)) -
    	4*a*e*l*M*pow(r1,3)*(4*pow(l,4)*M + 9*pow(l,2)*(12*M - 7*r1)*pow(r1,2) + 6*(9*M -
    	7*r1)*pow(r1,4)) + 4*pow(a,5)*e*l*M*r1*(48*pow(M,3) + 142*pow(M,2)*r1 + 76*M*pow(r1,2) +
    	50*pow(r1,3) + pow(l,2)*(-6*M + 9*r1)) - pow(a,8)*(24*pow(e,2)*pow(M,2)*(2*M + r1) +
    	r1*(49*pow(M,2) + 42*M*r1 + 6*pow(r1,2))) + pow(r1,3)*(4*pow(l,6)*M*(2*M - r1) +
    	3*pow(r1,6)*(40*pow(M,2) - 2*(21 + 2*pow(e,2))*M*r1 + 11*pow(r1,2)) +
    	6*pow(l,2)*pow(r1,4)*(56*pow(M,2) - 62*M*r1 + 17*pow(r1,2)) +
    	2*pow(l,4)*pow(r1,2)*(118*pow(M,2) + (-131 + 2*pow(e,2))*M*r1 + 36*pow(r1,2))) +
    	4*pow(a,3)*e*l*M*pow(r1,3)*(-6*pow(l,2)*(9*M - 8*r1) + r1*(114*pow(M,2) + 23*M*r1 +
    	83*pow(r1,2))) + 2*pow(a,2)*pow(r1,3)*(pow(l,4)*(4*(18 + pow(e,2))*pow(M,2) + (-71 +
    	2*pow(e,2))*M*r1 + 15*pow(r1,2)) - 3*pow(r1,3)*(12*pow(M,3) + 6*(-3 +
    	8*pow(e,2))*pow(M,2)*r1 + (31 + 6*pow(e,2))*M*pow(r1,2) - 13*pow(r1,3)) +
    	3*pow(l,2)*r1*(-32*pow(M,3) + 2*(19 + 9*pow(e,2))*pow(M,2)*r1 - 61*M*pow(r1,2) +
    	24*pow(r1,3))) + pow(a,4)*r1*(6*pow(l,4)*(2*pow(M,2) - 6*M*r1 - pow(r1,2)) -
    	3*pow(r1,2)*(76*pow(M,4) + 12*(-1 + 16*pow(e,2))*pow(M,3)*r1 + 3*(-19 +
    	56*pow(e,2))*pow(M,2)*pow(r1,2) + 4*(8 + 3*pow(e,2))*M*pow(r1,3) - 17*pow(r1,4)) -
    	2*pow(l,2)*(48*pow(M,4) + 68*pow(M,3)*r1 - 4*(13 + 9*pow(e,2))*pow(M,2)*pow(r1,2) +
    	84*M*pow(r1,3) - 15*pow(r1,4))) - 2*pow(a,6)*(pow(l,2)*(24*pow(M,3) + (26 -
    	6*pow(e,2))*pow(M,2)*r1 + 39*M*pow(r1,2) + 6*pow(r1,3)) + M*r1*(r1*(-106*pow(M,2) - 37*M*r1
    	+ 39*pow(r1,2)) + 6*pow(e,2)*(8*pow(M,3) + 36*pow(M,2)*r1 + 20*M*pow(r1,2) +
    	pow(r1,3)))))/(24.*pow(r1,3)*pow(pow(a,2) + r1*(-2*M + r1),2));
    A240 =
    	(3*pow(r1,6))/(pow(a,2) + r1*(-2*M + r1));
    A242 = (24*pow(a,5)*e*l*M*(7*M + 4*r1) +
    	12*pow(a,3)*e*l*M*r1*(4*pow(l,2) - 16*pow(M,2) + 9*M*r1 + 23*pow(r1,2)) -
    	12*a*e*l*M*pow(r1,2)*(pow(l,2)*(8*M - 4*r1) + (19*M - 13*r1)*pow(r1,2)) -
    	3*pow(a,6)*(4*pow(e,2)*M*(7*M + 3*r1) + r1*(10*M + 3*r1)) + pow(r1,2)*(24*pow(l,4)*M*(2*M
    	- r1) + pow(r1,4)*(120*pow(M,2) - 2*(71 + 6*pow(e,2))*M*r1 + 41*pow(r1,2)) +
    	pow(l,2)*pow(r1,2)*(228*pow(M,2) - 232*M*r1 + 59*pow(r1,2))) + pow(a,2)*r1*(-24*pow(l,4)*M
    	+ pow(r1,2)*(-36*pow(M,3) + (22 - 84*pow(e,2))*pow(M,2)*r1 - 6*(19 +
    	10*pow(e,2))*M*pow(r1,2) + 67*pow(r1,3)) + 2*pow(l,2)*(48*pow(M,3) - 6*(13 +
    	2*pow(e,2))*M*pow(r1,2) + 25*pow(r1,3) + 6*pow(M,2)*(r1 + 4*pow(e,2)*r1))) -
    	pow(a,4)*(3*pow(l,2)*(28*pow(M,2) + 4*(5 + 2*pow(e,2))*M*r1 + 3*pow(r1,2)) +
    	r1*(r1*(-78*pow(M,2) + 14*M*r1 - 17*pow(r1,2)) + 12*pow(e,2)*M*(-8*pow(M,2) + 10*M*r1 +
    	7*pow(r1,2)))))/(24.*pow(pow(a,2) + r1*(-2*M + r1),2));
    A260 =
    	(pow(r1,3)*(48*pow(a,3)*e*l*M - 3*pow(a,4)*(8*pow(e,2)*M + r1) + 16*a*e*l*M*r1*(-5*M +
    	3*r1) + r1*(20*pow(l,2)*M*(2*M - r1) + pow(r1,2)*(40*pow(M,2) - 2*(21 + 2*pow(e,2))*M*r1 +
    	11*pow(r1,2))) - 2*pow(a,2)*(12*pow(l,2)*M + r1*((7*M - 4*r1)*r1 + 2*pow(e,2)*M*(-10*M +
    	7*r1)))))/(24.*pow(pow(a,2) + r1*(-2*M + r1),2));
    A304 = -((pow(a,2)*(2*M + r1) +
    	r1*(pow(l,2) + pow(r1,2)))*(4*pow(a,3)*e*l*M + pow(a,4)*r1 + 12*a*e*l*M*pow(r1,2) +
    	pow(r1,2)*(pow(r1,2)*(-5*M + 2*r1) + pow(l,2)*(-9*M + 4*r1)) + pow(a,2)*(pow(l,2)*(-4*M +
    	r1) + r1*(2*pow(M,2) - 3*M*r1 + 3*pow(r1,2)))))/(2.*r1*pow(pow(a,2) + r1*(-2*M + r1),2));
    	A306 = (-4*pow(a,11)*e*l*pow(M,2)*(4*(35 + 9*pow(e,2))*pow(M,2) + (127 +
    	54*pow(e,2))*M*r1 + 3*(5 + 6*pow(e,2))*pow(r1,2)) -
    	4*a*e*l*M*pow(r1,5)*(pow(r1,6)*(-92*pow(M,2) + 3*(31 + 18*pow(e,2))*M*r1 - 24*pow(r1,2))
    	+ 21*pow(l,6)*(4*pow(M,2) - pow(r1,2)) + pow(l,2)*pow(r1,4)*(106*pow(M,2) + (-55 +
    	114*pow(e,2))*M*r1 + 3*pow(r1,2)) + pow(l,4)*pow(r1,2)*(486*pow(M,2) + 2*(-170 +
    	21*pow(e,2))*M*r1 + 51*pow(r1,2))) -
    	4*pow(a,3)*e*l*M*pow(r1,3)*(2*pow(l,4)*pow(r1,2)*((829 + 42*pow(e,2))*pow(M,2) + 2*(-107
    	+ 21*pow(e,2))*M*r1 - 84*pow(r1,2)) + 21*pow(l,6)*(12*pow(M,2) - 2*M*r1 - pow(r1,2)) +
    	pow(l,2)*pow(r1,3)*(-3230*pow(M,3) + 3*(895 + 224*pow(e,2))*pow(M,2)*r1 + 16*(-7 +
    	24*pow(e,2))*M*pow(r1,2) - 208*pow(r1,3)) + pow(r1,5)*(-1722*pow(M,3) + (943 -
    	270*pow(e,2))*pow(M,2)*r1 + (209 + 270*pow(e,2))*M*pow(r1,2) - 118*pow(r1,3))) +
    	pow(r1,7)*(2*pow(l,6)*(2*M - r1)*(54*pow(M,2) + (-71 + 42*pow(e,2))*M*r1 + 24*pow(r1,2)) +
    	pow(l,2)*pow(r1,4)*(220*pow(M,3) + 340*(-1 + pow(e,2))*pow(M,2)*r1 + (211 -
    	184*pow(e,2))*M*pow(r1,2) - 48*pow(r1,3)) + 2*pow(l,4)*pow(r1,2)*(136*pow(M,3) + 2*(-125 +
    	111*pow(e,2))*pow(M,2)*r1 + (163 - 113*pow(e,2))*M*pow(r1,2) - 36*pow(r1,3)) +
    	pow(r1,6)*(116*pow(M,3) + 8*(-20 + 11*pow(e,2))*pow(M,2)*r1 + 3*(29 -
    	18*pow(e,2))*M*pow(r1,2) - 18*pow(r1,3))) + pow(a,12)*M*(3*r1*(-8*pow(M,2) - 3*M*r1 +
    	2*pow(r1,2)) + pow(e,2)*(280*pow(M,3) + 308*pow(M,2)*r1 + 72*M*pow(r1,2) - 6*pow(r1,3))) -
    	4*pow(a,9)*e*l*M*(pow(l,2)*M*(108*pow(M,2) + 2*(137 + 24*pow(e,2))*M*r1 + 3*(13 +
    	6*pow(e,2))*pow(r1,2)) + r1*(24*(-1 + 12*pow(e,2))*pow(M,4) + 2*(215 +
    	234*pow(e,2))*pow(M,3)*r1 + 12*(49 + 39*pow(e,2))*pow(M,2)*pow(r1,2) + 9*(-5 +
    	18*pow(e,2))*M*pow(r1,3) - 30*pow(r1,4))) -
    	pow(a,2)*pow(r1,3)*(-168*pow(l,8)*pow(M,2)*(2*M - r1) + 2*pow(l,6)*r1*(240*pow(M,4) -
    	4*(233 + 84*pow(e,2))*pow(M,3)*r1 - 6*(-83 + 28*pow(e,2))*pow(M,2)*pow(r1,2) + (23 +
    	84*pow(e,2))*M*pow(r1,3) - 36*pow(r1,4)) + pow(l,2)*pow(r1,5)*(3640*pow(M,4) + 24*(-208 +
    	29*pow(e,2))*pow(M,3)*r1 + 2*(1091 - 1004*pow(e,2))*pow(M,2)*pow(r1,2) + 8*(-27 +
    	100*pow(e,2))*M*pow(r1,3) - 27*pow(r1,4)) + 4*pow(l,4)*pow(r1,3)*(1042*pow(M,4) - (1441 +
    	1086*pow(e,2))*pow(M,3)*r1 + 2*(293 + 3*pow(e,2))*pow(M,2)*pow(r1,2) + (-19 +
    	183*pow(e,2))*M*pow(r1,3) - 21*pow(r1,4)) + pow(r1,7)*(696*pow(M,4) + 8*(-165 +
    	137*pow(e,2))*pow(M,3)*r1 + 4*(204 - 241*pow(e,2))*pow(M,2)*pow(r1,2) + (-181 +
    	268*pow(e,2))*M*pow(r1,3) + 15*pow(r1,4))) -
    	4*pow(a,5)*e*l*M*pow(r1,2)*(pow(l,4)*(1008*pow(M,3) + 4*(257 + 21*pow(e,2))*pow(M,2)*r1
    	+ 2*(-104 + 21*pow(e,2))*M*pow(r1,2) - 87*pow(r1,3)) + pow(l,2)*r1*(-944*pow(M,4) +
    	10*(-161 + 24*pow(e,2))*pow(M,3)*r1 + 144*(21 + 8*pow(e,2))*pow(M,2)*pow(r1,2) + 6*(-89 +
    	74*pow(e,2))*M*pow(r1,3) - 307*pow(r1,4)) - 2*pow(r1,3)*(192*pow(M,4) + (823 +
    	378*pow(e,2))*pow(M,3)*r1 + 2*(-544 + 9*pow(e,2))*pow(M,2)*pow(r1,2) -
    	252*pow(e,2)*M*pow(r1,3) + 97*pow(r1,4))) - 4*pow(a,7)*e*l*M*r1*(24*pow(l,4)*M*(6*M +
    	r1) + pow(l,2)*(864*pow(M,4) + 4*(163 + 84*pow(e,2))*pow(M,3)*r1 + (1271 +
    	624*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-101 + 96*pow(e,2))*M*pow(r1,3) - 96*pow(r1,4)) -
    	2*r1*(224*pow(M,5) + 460*pow(M,4)*r1 + 9*(17 + 2*pow(e,2))*pow(M,3)*pow(r1,2) - (769 +
    	324*pow(e,2))*pow(M,2)*pow(r1,3) + 8*(11 - 27*pow(e,2))*M*pow(r1,4) + 65*pow(r1,5))) +
    	pow(a,4)*pow(r1,2)*(2*pow(l,6)*M*(672*pow(M,3) + 4*(101 + 126*pow(e,2))*pow(M,2)*r1 +
    	28*(-11 + 3*pow(e,2))*M*pow(r1,2) - (31 + 42*pow(e,2))*pow(r1,3)) +
    	4*pow(l,4)*r1*(-472*pow(M,5) + 2*(-239 + 180*pow(e,2))*pow(M,4)*r1 + 8*(163 +
    	293*pow(e,2))*pow(M,3)*pow(r1,2) + 3*(-223 + 58*pow(e,2))*pow(M,2)*pow(r1,3) + (41 -
    	198*pow(e,2))*M*pow(r1,4) + 39*pow(r1,5)) + pow(r1,5)*(180*pow(M,5) - 12*(-9 +
    	238*pow(e,2))*pow(M,4)*r1 + 3*(491 - 180*pow(e,2))*pow(M,3)*pow(r1,2) + 8*(-234 +
    	287*pow(e,2))*pow(M,2)*pow(r1,3) + (403 - 538*pow(e,2))*M*pow(r1,4) + 63*pow(r1,5)) +
    	pow(l,2)*pow(r1,3)*(-364*pow(M,5) - 4*(763 + 2976*pow(e,2))*pow(M,4)*r1 + (6105 +
    	4072*pow(e,2))*pow(M,3)*pow(r1,2) + 4*(-989 + 831*pow(e,2))*pow(M,2)*pow(r1,3) + 3*(181 -
    	436*pow(e,2))*M*pow(r1,4) + 198*pow(r1,5))) + pow(a,6)*r1*(6*pow(l,6)*M*(32*pow(M,2) +
    	2*M*r1 + pow(r1,2)) + 4*pow(l,4)*M*(288*pow(M,4) + 4*(23 + 252*pow(e,2))*pow(M,3)*r1 +
    	(389 + 1450*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-125 + 69*pow(e,2))*M*pow(r1,3) - (6 +
    	73*pow(e,2))*pow(r1,4)) + pow(r1,3)*(184*pow(M,6) - 4*(11 + 268*pow(e,2))*pow(M,5)*r1 -
    	6*(37 + 632*pow(e,2))*pow(M,4)*pow(r1,2) + 5*(-11 + 416*pow(e,2))*pow(M,3)*pow(r1,3) +
    	(-1303 + 2360*pow(e,2))*pow(M,2)*pow(r1,4) + 3*(219 - 184*pow(e,2))*M*pow(r1,5) +
    	99*pow(r1,6)) + pow(l,2)*r1*(-896*pow(M,6) - 32*(38 + 59*pow(e,2))*pow(M,5)*r1 - 16*(21 +
    	292*pow(e,2))*pow(M,4)*pow(r1,2) + 2*(1643 + 4488*pow(e,2))*pow(M,3)*pow(r1,3) + (-1979 +
    	2372*pow(e,2))*pow(M,2)*pow(r1,4) + 4*(154 - 241*pow(e,2))*M*pow(r1,5) + 123*pow(r1,6))) +
    	pow(a,10)*M*(pow(l,2)*(8*(35 + 54*pow(e,2))*pow(M,3) + 8*(25 + 109*pow(e,2))*pow(M,2)*r1
    	+ 3*(-7 + 76*pow(e,2))*M*pow(r1,2) + 6*(3 - 2*pow(e,2))*pow(r1,3)) + r1*(r1*(134*pow(M,3)
    	- 21*pow(M,2)*r1 - 69*M*pow(r1,2) + 58*pow(r1,3)) + pow(e,2)*(-48*pow(M,4) +
    	1072*pow(M,3)*r1 + 1704*pow(M,2)*pow(r1,2) + 452*M*pow(r1,3) - 76*pow(r1,4)))) +
    	pow(a,8)*(2*pow(l,4)*M*(72*pow(M,3) + 4*(55 + 72*pow(e,2))*pow(M,2)*r1 +
    	78*pow(e,2)*M*pow(r1,2) - 3*(-3 + pow(e,2))*pow(r1,3)) + pow(l,2)*M*r1*(48*(-1 +
    	72*pow(e,2))*pow(M,4) + 8*(81 + 514*pow(e,2))*pow(M,3)*r1 + (869 +
    	5400*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-285 + 472*pow(e,2))*M*pow(r1,3) + 4*(24 -
    	71*pow(e,2))*pow(r1,4)) + pow(r1,2)*(r1*(-264*pow(M,5) + 92*pow(M,4)*r1 +
    	135*pow(M,3)*pow(r1,2) - 147*pow(M,2)*pow(r1,3) + 400*M*pow(r1,4) + 39*pow(r1,5)) -
    	2*pow(e,2)*M*(448*pow(M,5) + 1232*pow(M,4)*r1 + 680*pow(M,3)*pow(r1,2) -
    	1460*pow(M,2)*pow(r1,3) - 660*M*pow(r1,4) + 149*pow(r1,5)))))/(48.*pow(r1,7)*pow(pow(a,2) +
    	r1*(-2*M + r1),3));
    A322 = -((pow(r1,2)*(2*pow(a,3)*e*l*M + 6*a*e*l*M*pow(r1,2) +
    	pow(a,4)*(3*M + 2*r1) + pow(r1,2)*(pow(r1,2)*(-5*M + 2*r1) + pow(l,2)*(-7*M + 3*r1)) -
    	2*pow(a,2)*(pow(l,2)*(M - r1) + r1*(2*pow(M,2) + M*r1 - 2*pow(r1,2)))))/pow(pow(a,2) +
    	r1*(-2*M + r1),2));
    A324 = (4*pow(a,9)*e*l*M*(2*(-113 + 9*pow(e,2))*pow(M,2) - 57*M*r1 +
    	3*pow(r1,2)) - 4*a*e*l*M*pow(r1,4)*(pow(l,4)*(240*pow(M,2) - 22*M*r1 - 51*pow(r1,2)) +
    	pow(l,2)*pow(r1,2)*(314*pow(M,2) + (-83 + 114*pow(e,2))*M*r1 - 42*pow(r1,2)) +
    	3*pow(r1,4)*(-84*pow(M,2) + 4*(23 + 9*pow(e,2))*M*r1 - 27*pow(r1,2))) +
    	pow(a,10)*M*(-21*r1*(7*M + 2*r1) + 4*pow(e,2)*(113*pow(M,2) + 81*M*r1 + 9*pow(r1,2))) -
    	4*pow(a,7)*e*l*M*(3*pow(l,2)*(-18*pow(M,2) + 8*M*r1 + pow(r1,2)) + r1*(4*(31 +
    	108*pow(e,2))*pow(M,3) + 9*(109 + 30*pow(e,2))*pow(M,2)*r1 + (37 +
    	108*pow(e,2))*M*pow(r1,2) - 69*pow(r1,3))) -
    	4*pow(a,3)*e*l*M*pow(r1,2)*(pow(l,4)*(612*pow(M,2) - 92*M*r1 - 45*pow(r1,2)) +
    	pow(r1,3)*(-2250*pow(M,3) + (955 - 54*pow(e,2))*pow(M,2)*r1 + 3*(143 +
    	108*pow(e,2))*M*pow(r1,2) - 187*pow(r1,3)) + pow(l,2)*r1*(16*pow(M,3) + 2*(409 +
    	114*pow(e,2))*pow(M,2)*r1 + 2*(-185 + 108*pow(e,2))*M*pow(r1,2) - 45*pow(r1,3))) +
    	pow(r1,4)*(2*pow(l,6)*M*(12*pow(M,2) - 16*M*r1 + 5*pow(r1,2)) +
    	2*pow(l,4)*pow(r1,2)*(128*pow(M,3) + 2*(-113 + 111*pow(e,2))*pow(M,2)*r1 + (129 -
    	113*pow(e,2))*M*pow(r1,2) - 24*pow(r1,3)) + 4*pow(l,2)*pow(r1,4)*(76*pow(M,3) + 2*(-66 +
    	85*pow(e,2))*pow(M,2)*r1 + (77 - 92*pow(e,2))*M*pow(r1,2) - 15*pow(r1,3)) +
    	3*pow(r1,6)*(116*pow(M,3) + 8*(-21 + 11*pow(e,2))*pow(M,2)*r1 + (83 -
    	54*pow(e,2))*M*pow(r1,2) - 14*pow(r1,3))) - 4*pow(a,5)*e*l*M*r1*(6*pow(l,4)*r1 +
    	pow(l,2)*(1296*pow(M,3) + 6*(47 + 34*pow(e,2))*pow(M,2)*r1 + (-371 +
    	102*pow(e,2))*M*pow(r1,2) - 60*pow(r1,3)) + r1*(-1008*pow(M,4) + 2*(-641 +
    	216*pow(e,2))*pow(M,3)*r1 + 6*(335 + 63*pow(e,2))*pow(M,2)*pow(r1,2) + (209 +
    	324*pow(e,2))*M*pow(r1,3) - 172*pow(r1,4))) +
    	pow(a,2)*pow(r1,2)*(2*pow(l,6)*M*(408*pow(M,2) - 194*M*r1 - 5*pow(r1,2)) -
    	8*pow(l,4)*M*r1*(104*pow(M,3) - 3*(73 + 77*pow(e,2))*pow(M,2)*r1 + (115 -
    	101*pow(e,2))*M*pow(r1,2) + 17*(-1 + 3*pow(e,2))*pow(r1,3)) +
    	2*pow(l,2)*pow(r1,3)*(-2348*pow(M,4) + 4*(672 + 127*pow(e,2))*pow(M,3)*r1 + 3*(-355 +
    	324*pow(e,2))*pow(M,2)*pow(r1,2) + 3*(45 - 176*pow(e,2))*M*pow(r1,3) + 9*pow(r1,4)) -
    	3*pow(r1,5)*(464*pow(M,4) + 8*(-100 + 97*pow(e,2))*pow(M,3)*r1 + 4*(139 -
    	166*pow(e,2))*pow(M,2)*pow(r1,2) + 4*(-39 + 49*pow(e,2))*M*pow(r1,3) + 11*pow(r1,4))) +
    	pow(a,4)*r1*(12*pow(l,6)*M*r1 + 3*pow(r1,3)*(60*pow(M,5) - 4*(23 +
    	326*pow(e,2))*pow(M,4)*r1 + (751 - 436*pow(e,2))*pow(M,3)*pow(r1,2) + 4*(-208 +
    	293*pow(e,2))*pow(M,2)*pow(r1,3) + (205 - 252*pow(e,2))*M*pow(r1,4) + 30*pow(r1,5)) +
    	2*pow(l,2)*r1*(-1008*pow(M,5) + 4*(-227 + 328*pow(e,2))*pow(M,4)*r1 + 2*(1325 +
    	758*pow(e,2))*pow(M,3)*pow(r1,2) + 14*(-95 + 6*pow(e,2))*pow(M,2)*pow(r1,3) + 3*(43 -
    	160*pow(e,2))*M*pow(r1,4) + 78*pow(r1,5)) + 2*pow(l,4)*(864*pow(M,4) + 2*(-193 +
    	107*pow(e,2))*pow(M,2)*pow(r1,2) + 17*(1 - 5*pow(e,2))*M*pow(r1,3) + 24*pow(r1,4) +
    	12*pow(M,3)*(r1 + 102*pow(e,2)*r1))) + pow(a,8)*(2*pow(l,2)*M*((226 -
    	108*pow(e,2))*pow(M,2) + 24*(-2 + pow(e,2))*M*r1 + 3*(-15 + 8*pow(e,2))*pow(r1,2)) +
    	r1*(2*pow(e,2)*M*(124*pow(M,3) + 982*pow(M,2)*r1 + 686*M*pow(r1,2) + 3*pow(r1,3)) +
    	3*r1*(179*pow(M,3) - 70*pow(M,2)*r1 + 30*M*pow(r1,2) + 10*pow(r1,3)))) +
    	pow(a,6)*(-12*pow(l,4)*M*(6*pow(M,2) - 4*M*r1 - (-3 + pow(e,2))*pow(r1,2)) +
    	2*pow(l,2)*r1*(4*(31 + 648*pow(e,2))*pow(M,4) + 28*(35 + 39*pow(e,2))*pow(M,3)*r1 - (361
    	+ 140*pow(e,2))*pow(M,2)*pow(r1,2) + (43 - 112*pow(e,2))*M*pow(r1,3) + 39*pow(r1,4)) +
    	pow(r1,2)*(3*r1*(-192*pow(M,4) + 210*pow(M,3)*r1 - 381*pow(M,2)*pow(r1,2) + 176*M*pow(r1,3)
    	+ 37*pow(r1,4)) - 4*pow(e,2)*M*(504*pow(M,4) + 828*pow(M,3)*r1 - 433*pow(M,2)*pow(r1,2) -
    	709*M*pow(r1,3) + 90*pow(r1,4)))))/(48.*pow(r1,4)*pow(pow(a,2) + r1*(-2*M + r1),3));
    A340 =
    	-(pow(r1,5)*(3*pow(a,2) + r1*(-5*M + 2*r1)))/(2.*pow(pow(a,2) + r1*(-2*M + r1),2));
    A342 =
    	(24*pow(a,7)*e*l*M*((-25 + 3*pow(e,2))*M - 7*r1) + 3*pow(a,8)*(r1*(8*M + 9*r1) +
    	2*pow(e,2)*M*(50*M + 19*r1)) - 4*a*e*l*M*pow(r1,3)*(3*pow(r1,2)*(-76*pow(M,2) + (103 +
    	18*pow(e,2))*M*r1 - 36*pow(r1,2)) + pow(l,2)*(208*pow(M,2) - 64*M*r1 - 27*pow(r1,2))) -
    	12*pow(a,5)*e*l*M*(pow(l,2)*(-18*M + 8*r1) + r1*((4 + 48*pow(e,2))*pow(M,2) + (71 -
    	6*pow(e,2))*M*r1 + 20*pow(r1,2))) - 12*pow(a,3)*e*l*M*r1*(pow(l,2)*(144*pow(M,2) -
    	64*M*r1 - pow(r1,2)) + 2*r1*(-88*pow(M,3) + 2*(-5 + 18*pow(e,2))*pow(M,2)*r1 + 9*(4 +
    	pow(e,2))*M*pow(r1,2) - 13*pow(r1,3))) + pow(r1,3)*(-8*pow(l,4)*M*(2*pow(M,2) - 3*M*r1 +
    	pow(r1,2)) + pow(r1,4)*(348*pow(M,3) + 8*(-64 + 33*pow(e,2))*pow(M,2)*r1 + (245 -
    	162*pow(e,2))*M*pow(r1,2) - 38*pow(r1,3)) + pow(l,2)*pow(r1,2)*(-52*pow(M,3) + 4*(19 +
    	85*pow(e,2))*pow(M,2)*r1 - (33 + 184*pow(e,2))*M*pow(r1,2) + 4*pow(r1,3))) -
    	pow(a,2)*r1*(-12*pow(l,4)*M*(48*pow(M,2) - 29*M*r1 + pow(r1,2)) +
    	pow(l,2)*r1*(1056*pow(M,4) - 16*(33 + 107*pow(e,2))*pow(M,3)*r1 + 4*(-35 +
    	16*pow(e,2))*pow(M,2)*pow(r1,2) + 4*(9 + 64*pow(e,2))*M*pow(r1,3) + 25*pow(r1,4)) +
    	pow(r1,3)*(696*pow(M,4) + 8*(-103 + 171*pow(e,2))*pow(M,3)*r1 + 4*(166 -
    	273*pow(e,2))*pow(M,2)*pow(r1,2) + 3*(-111 + 124*pow(e,2))*M*pow(r1,3) + 57*pow(r1,4))) +
    	pow(a,6)*(pow(l,2)*((300 - 216*pow(e,2))*pow(M,2) + 6*(9 + 8*pow(e,2))*M*r1 +
    	27*pow(r1,2)) + r1*(12*pow(e,2)*M*(2*pow(M,2) + 61*M*r1 + 15*pow(r1,2)) + r1*(-312*pow(M,2)
    	- 25*M*r1 + 61*pow(r1,2)))) + pow(a,4)*(24*pow(l,4)*M*(-3*M + 2*r1) + pow(l,2)*r1*(24*(1 +
    	72*pow(e,2))*pow(M,3) + 12*(10 - 41*pow(e,2))*pow(M,2)*r1 - 3*(23 +
    	8*pow(e,2))*M*pow(r1,2) - 2*pow(r1,3)) + pow(r1,2)*(-12*pow(e,2)*M*(88*pow(M,3) +
    	64*pow(M,2)*r1 - 97*M*pow(r1,2) + 12*pow(r1,3)) + r1*(876*pow(M,3) - 512*pow(M,2)*r1 +
    	63*M*pow(r1,2) + 15*pow(r1,3)))))/(48.*r1*pow(pow(a,2) + r1*(-2*M + r1),3));
    A360 =
    	(pow(r1,2)*(-144*pow(a,5)*e*l*M + 28*pow(a,3)*e*l*M*(4*M - 3*r1)*r1 +
    	3*pow(a,6)*(24*pow(e,2)*M + 7*r1) + 4*a*e*l*M*pow(r1,2)*(68*pow(M,2) - 54*M*r1 +
    	15*pow(r1,2)) + pow(r1,2)*(-2*pow(l,2)*M*(68*pow(M,2) - 76*M*r1 + 21*pow(r1,2)) +
    	pow(r1,2)*(116*pow(M,3) + 8*(-21 + 11*pow(e,2))*pow(M,2)*r1 + (83 -
    	54*pow(e,2))*M*pow(r1,2) - 14*pow(r1,3))) + pow(a,4)*(72*pow(l,2)*M + r1*(r1*(-53*M +
    	28*r1) + pow(e,2)*(-56*pow(M,2) + 74*M*r1))) - pow(a,2)*r1*(2*pow(l,2)*M*(28*M - 5*r1) +
    	r1*(r1*(36*pow(M,2) - 30*M*r1 + 7*pow(r1,2)) + 4*pow(e,2)*M*(34*pow(M,2) - 16*M*r1 +
    	13*pow(r1,2))))))/(48.*pow(pow(a,2) + r1*(-2*M + r1),3));
    A402 = (3*pow(r1,3)*(pow(a,2)*(2*M
    	+ r1) + r1*(pow(l,2) + pow(r1,2))))/pow(pow(a,2) + r1*(-2*M + r1),2);
    A404 =
    	(-4*pow(a,7)*e*l*M*(48*pow(M,2) + 40*M*r1 + 3*pow(r1,2)) -
    	4*pow(a,5)*e*l*M*r1*(-96*pow(M,3) - 110*pow(M,2)*r1 + 73*M*pow(r1,2) + 8*pow(r1,3) +
    	9*pow(l,2)*(6*M + r1)) + pow(a,8)*(12*pow(e,2)*M*(8*pow(M,2) + 10*M*r1 + 3*pow(r1,2)) -
    	r1*(62*pow(M,2) + 42*M*r1 + 3*pow(r1,2))) + pow(r1,3)*(4*pow(l,6)*M*(-2*M + r1) +
    	3*pow(r1,6)*(23*pow(M,2) + 2*(-9 + 2*pow(e,2))*M*r1 + 2*pow(r1,2)) +
    	6*pow(l,2)*pow(r1,4)*(21*pow(M,2) + (-17 + 4*pow(e,2))*M*r1 + 2*pow(r1,2)) +
    	pow(l,4)*pow(r1,2)*(73*pow(M,2) + 4*(-17 + 2*pow(e,2))*M*r1 + 12*pow(r1,2))) -
    	4*pow(a,3)*e*l*M*pow(r1,2)*(6*pow(l,4) + pow(l,2)*(-48*pow(M,2) + 81*M*r1 + 9*pow(r1,2))
    	+ pow(r1,2)*(-234*pow(M,2) + 50*M*r1 + 17*pow(r1,2))) +
    	4*a*e*l*M*pow(r1,3)*(pow(l,4)*(4*M - 6*r1) + 3*(M - 4*r1)*pow(r1,4) +
    	pow(l,2)*(-45*M*pow(r1,2) + 6*pow(r1,3))) + 2*pow(a,2)*pow(r1,2)*(6*pow(l,6)*M +
    	3*pow(r1,4)*(-16*pow(M,3) - 4*(-7 + 9*pow(e,2))*pow(M,2)*r1 + (-25 +
    	12*pow(e,2))*M*pow(r1,2) + 5*pow(r1,3)) + 3*pow(l,2)*pow(r1,2)*(-76*pow(M,3) + 2*(35 +
    	9*pow(e,2))*pow(M,2)*r1 + (-33 + 16*pow(e,2))*M*pow(r1,2) + 9*pow(r1,3)) +
    	pow(l,4)*(-48*pow(M,3) + (78 - 4*pow(e,2))*pow(M,2)*r1 + 5*(-5 + 2*pow(e,2))*M*pow(r1,2)
    	+ 12*pow(r1,3))) + pow(a,4)*r1*(3*pow(l,4)*(36*pow(M,2) + 4*(-1 + pow(e,2))*M*r1 -
    	pow(r1,2)) - 2*pow(l,2)*(96*pow(M,4) + 4*(19 + 12*pow(e,2))*pow(M,3)*r1 - (143 +
    	84*pow(e,2))*pow(M,2)*pow(r1,2) + 15*(3 - 4*pow(e,2))*M*pow(r1,3) - 18*pow(r1,4)) +
    	3*pow(r1,2)*(-36*pow(M,4) - 4*(-5 + 48*pow(e,2))*pow(M,3)*r1 - 9*(-5 +
    	8*pow(e,2))*pow(M,2)*pow(r1,2) + 8*(-5 + 6*pow(e,2))*M*pow(r1,3) + 13*pow(r1,4))) +
    	2*pow(a,6)*(pow(l,2)*(48*pow(M,3) + 2*(10 + 27*pow(e,2))*pow(M,2)*r1 + 3*(-11 +
    	8*pow(e,2))*M*pow(r1,2) - 3*pow(r1,3)) + r1*(-12*pow(e,2)*M*(8*pow(M,3) + 12*pow(M,2)*r1 -
    	5*M*pow(r1,2) - 5*pow(r1,3)) + r1*(74*pow(M,3) + 17*pow(M,2)*r1 - 33*M*pow(r1,2) +
    	6*pow(r1,3)))))/(24.*pow(r1,3)*pow(pow(a,2) + r1*(-2*M + r1),3));
    A420 =
    	(3*pow(r1,6))/pow(pow(a,2) + r1*(-2*M + r1),2);
    A422 = (2*pow(a,5)*e*l*M*(-8*M + r1) -
    	2*pow(a,3)*e*l*M*r1*(4*pow(l,2) - 8*pow(M,2) + M*r1 - 7*pow(r1,2)) -
    	2*a*e*l*M*pow(r1,3)*(4*pow(l,2) + r1*(3*M + 4*r1)) + pow(a,6)*M*(-7*r1 + pow(e,2)*(8*M +
    	6*r1)) + pow(r1,4)*(pow(r1,2)*(23*pow(M,2) + 2*(-11 + 2*pow(e,2))*M*r1 + 4*pow(r1,2)) +
    	pow(l,2)*(25*pow(M,2) + (-23 + 4*pow(e,2))*M*r1 + 4*pow(r1,2))) +
    	pow(a,2)*r1*(4*pow(l,4)*M + pow(r1,2)*(-16*pow(M,3) - 4*(-5 + 8*pow(e,2))*pow(M,2)*r1 +
    	7*(-5 + 2*pow(e,2))*M*pow(r1,2) + 14*pow(r1,3)) + pow(l,2)*(-8*pow(M,3) + 10*pow(M,2)*r1 +
    	(-23 + 8*pow(e,2))*M*pow(r1,2) + 15*pow(r1,3))) + pow(a,4)*(4*pow(l,2)*M*(2*M + (-2 +
    	pow(e,2))*r1) + r1*(-8*pow(e,2)*M*(pow(M,2) + M*r1 - 2*pow(r1,2)) + r1*(17*pow(M,2) -
    	8*M*r1 + 10*pow(r1,2)))))/(4.*pow(pow(a,2) + r1*(-2*M + r1),3));
    A440 =
    	(pow(r1,4)*(9*pow(a,4) - 48*a*e*l*pow(M,2) + 12*pow(l,2)*M*(2*M - r1) +
    	pow(r1,2)*(69*pow(M,2) + 2*(-31 + 6*pow(e,2))*M*r1 + 10*pow(r1,2)) +
    	2*pow(a,2)*(6*pow(e,2)*M*(2*M + r1) + r1*(-30*M + 17*r1))))/(24.*pow(pow(a,2) + r1*(-2*M +
    	r1),3));
    A502 = -(pow(r1,2)*(2*pow(a,3)*e*l*M + 6*a*e*l*M*pow(r1,2) + pow(a,4)*(3*M +
    	2*r1) + pow(r1,2)*(pow(r1,2)*(-4*M + r1) + pow(l,2)*(-6*M + 2*r1)) +
    	pow(a,2)*(-2*pow(l,2)*(M - r1) + r1*(-2*pow(M,2) - 3*M*r1 +
    	3*pow(r1,2)))))/(2.*pow(pow(a,2) + r1*(-2*M + r1),3));
    A504 = (-2*pow(a,9)*e*l*M*(4*(49 +
    	27*pow(e,2))*pow(M,2) + 6*(-7 + 12*pow(e,2))*M*r1 - 9*pow(r1,2)) +
    	pow(a,10)*M*(3*r1*(-10*M + 7*r1) + 2*pow(e,2)*(98*pow(M,2) + 30*M*r1 - 15*pow(r1,2))) -
    	2*a*e*l*M*pow(r1,4)*(pow(l,2)*pow(r1,2)*(-725*pow(M,2) + 8*(172 + 21*pow(e,2))*M*r1 -
    	570*pow(r1,2)) + 3*pow(r1,4)*(-135*pow(M,2) + 8*(23 + 9*pow(e,2))*M*r1 - 76*pow(r1,2)) +
    	16*pow(l,4)*(2*pow(M,2) + 15*M*r1 - 9*pow(r1,2))) -
    	2*pow(a,7)*e*l*M*(pow(l,2)*(324*pow(M,2) + 12*(5 + 3*pow(e,2))*M*r1 + 9*pow(r1,2)) +
    	2*r1*(16*(2 + 9*pow(e,2))*pow(M,3) + 6*(94 + 57*pow(e,2))*pow(M,2)*r1 + 2*(-65 +
    	126*pow(e,2))*M*pow(r1,2) - 39*pow(r1,3))) + pow(r1,4)*(-2*pow(l,6)*M*(44*pow(M,2) -
    	52*M*r1 + 15*pow(r1,2)) + 3*pow(r1,6)*(52*pow(M,3) + (-67 + 44*pow(e,2))*pow(M,2)*r1 +
    	2*(14 - 15*pow(e,2))*M*pow(r1,2) - 2*pow(r1,3)) + 2*pow(l,4)*pow(r1,2)*(-153*pow(M,3) +
    	(271 + 140*pow(e,2))*pow(M,2)*r1 - 2*(82 + 39*pow(e,2))*M*pow(r1,2) + 36*pow(r1,3)) +
    	pow(l,2)*pow(r1,4)*(34*pow(M,3) + (21 + 428*pow(e,2))*pow(M,2)*r1 - 2*(35 +
    	124*pow(e,2))*M*pow(r1,2) + 36*pow(r1,3))) -
    	2*pow(a,3)*e*l*M*pow(r1,2)*(2*pow(l,4)*(252*pow(M,2) + 116*M*r1 - 63*pow(r1,2)) +
    	pow(l,2)*r1*(-784*pow(M,3) + 3*(219 + 40*pow(e,2))*pow(M,2)*r1 + 18*(39 +
    	22*pow(e,2))*M*pow(r1,2) - 168*pow(r1,3)) - 2*pow(r1,3)*(765*pow(M,3) + 2*(8 +
    	135*pow(e,2))*pow(M,2)*r1 - 9*(29 + 44*pow(e,2))*M*pow(r1,2) + 107*pow(r1,3))) -
    	2*pow(a,5)*e*l*M*r1*(18*pow(l,4)*(6*M + r1) + pow(l,2)*(864*pow(M,3) + 84*(15 +
    	2*pow(e,2))*pow(M,2)*r1 + 2*(-145 + 132*pow(e,2))*M*pow(r1,2) - 123*pow(r1,3)) +
    	r1*(-576*pow(M,4) - 1390*pow(M,3)*r1 + 3*(667 + 108*pow(e,2))*pow(M,2)*pow(r1,2) + 4*(89 +
    	252*pow(e,2))*M*pow(r1,3) - 55*pow(r1,4))) -
    	pow(a,2)*pow(r1,2)*(4*pow(l,6)*M*(-84*pow(M,2) + 8*M*r1 + 23*pow(r1,2)) +
    	3*pow(r1,5)*(210*pow(M,4) + (-375 + 484*pow(e,2))*pow(M,3)*r1 - 4*(-71 +
    	99*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-61 + 65*pow(e,2))*M*pow(r1,3) + 22*pow(r1,4)) +
    	2*pow(l,4)*r1*(392*pow(M,4) - (409 + 196*pow(e,2))*pow(M,3)*r1 - (39 +
    	584*pow(e,2))*pow(M,2)*pow(r1,2) + (-41 + 167*pow(e,2))*M*pow(r1,3) + 84*pow(r1,4)) +
    	pow(l,2)*pow(r1,3)*(1618*pow(M,4) + (-1559 + 2272*pow(e,2))*pow(M,3)*r1 + (707 -
    	4036*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-211 + 424*pow(e,2))*M*pow(r1,3) + 168*pow(r1,4)))
    	+ pow(a,4)*r1*(18*pow(l,6)*M*(4*M + r1) + 8*pow(l,4)*(72*pow(M,4) + 18*(4 +
    	7*pow(e,2))*pow(M,3)*r1 + (-91 + 128*pow(e,2))*pow(M,2)*pow(r1,2) - 4*(-4 +
    	5*pow(e,2))*M*pow(r1,3) + 9*pow(r1,4)) + 3*pow(r1,3)*(44*pow(M,5) - 56*(-1 +
    	8*pow(e,2))*pow(M,4)*r1 + 2*(233 - 192*pow(e,2))*pow(M,3)*pow(r1,2) + (-525 +
    	736*pow(e,2))*pow(M,2)*pow(r1,3) - 4*(-38 + 55*pow(e,2))*M*pow(r1,4) - 23*pow(r1,5)) -
    	pow(l,2)*r1*(576*pow(M,5) + 4*(263 + 196*pow(e,2))*pow(M,4)*r1 - 13*(215 +
    	88*pow(e,2))*pow(M,3)*pow(r1,2) - 6*(-313 + 524*pow(e,2))*pow(M,2)*pow(r1,3) + 4*(-85 +
    	241*pow(e,2))*M*pow(r1,4) + 87*pow(r1,5))) + pow(a,8)*(pow(l,2)*M*(4*(49 +
    	162*pow(e,2))*pow(M,2) + 12*(-12 + 23*pow(e,2))*M*r1 + 3*(11 - 4*pow(e,2))*pow(r1,2)) +
    	r1*(2*pow(e,2)*M*(32*pow(M,3) + 608*pow(M,2)*r1 + 286*M*pow(r1,2) - 105*pow(r1,3)) +
    	3*r1*(24*pow(M,3) + 2*pow(M,2)*r1 + 130*M*pow(r1,2) + 15*pow(r1,3)))) +
    	pow(a,6)*(6*pow(l,4)*M*(36*pow(M,2) + 2*(-1 + 18*pow(e,2))*M*r1 + (5 +
    	3*pow(e,2))*pow(r1,2)) + pow(l,2)*r1*(64*(1 + 27*pow(e,2))*pow(M,4) + 16*(65 +
    	207*pow(e,2))*pow(M,3)*r1 + 68*(-14 + 17*pow(e,2))*pow(M,2)*pow(r1,2) + (595 -
    	376*pow(e,2))*M*pow(r1,3) + 117*pow(r1,4)) + pow(r1,2)*(3*r1*(-58*pow(M,4) - 45*pow(M,3)*r1
    	- 420*pow(M,2)*pow(r1,2) + 181*M*pow(r1,3) + 12*pow(r1,4)) - 4*pow(e,2)*M*(144*pow(M,4) +
    	432*pow(M,3)*r1 - 218*pow(M,2)*pow(r1,2) - 416*M*pow(r1,3) +
    	135*pow(r1,4)))))/(48.*pow(r1,4)*pow(pow(a,2) + r1*(-2*M + r1),4));
    A520 =
    	-(pow(r1,5)*(3*pow(a,2) + r1*(-4*M + r1)))/(2.*pow(pow(a,2) + r1*(-2*M + r1),3));
    A522 =
    	(-6*pow(a,7)*e*l*M*(4*(14 + 3*pow(e,2))*M + 23*r1) + 3*pow(a,8)*(2*pow(e,2)*M*(28*M -
    	r1) + r1*(43*M + 9*r1)) - 2*a*e*l*M*pow(r1,3)*(pow(l,2)*(64*pow(M,2) + 344*M*r1 -
    	222*pow(r1,2)) + 3*pow(r1,2)*(-199*pow(M,2) + 4*(79 + 18*pow(e,2))*M*r1 - 134*pow(r1,2))) -
    	6*pow(a,5)*e*l*M*(4*pow(l,2)*(9*M + 2*r1) + r1*(12*(7 + 6*pow(e,2))*pow(M,2) + 6*(17 +
    	16*pow(e,2))*M*r1 + 7*pow(r1,2))) - 6*pow(a,3)*e*l*M*r1*(2*pow(l,2)*(108*pow(M,2) +
    	52*M*r1 - 33*pow(r1,2)) + r1*(-272*pow(M,3) + (-73 + 72*pow(e,2))*pow(M,2)*r1 + 2*(103 +
    	78*pow(e,2))*M*pow(r1,2) - 116*pow(r1,3))) + pow(r1,3)*(-4*pow(l,4)*M*(38*pow(M,2) -
    	45*M*r1 + 13*pow(r1,2)) + 6*pow(r1,4)*(52*pow(M,3) + (-63 + 44*pow(e,2))*pow(M,2)*r1 -
    	6*(-3 + 5*pow(e,2))*M*pow(r1,2) + 2*pow(r1,3)) + pow(l,2)*pow(r1,2)*(-158*pow(M,3) + (369
    	+ 428*pow(e,2))*pow(M,2)*r1 - 4*(73 + 62*pow(e,2))*M*pow(r1,2) + 84*pow(r1,3))) +
    	pow(a,2)*r1*(24*pow(l,4)*M*(18*pow(M,2) + M*r1 - 7*pow(r1,2)) + pow(l,2)*r1*(-816*pow(M,4)
    	+ 2*(213 + 356*pow(e,2))*pow(M,3)*r1 + (103 + 1444*pow(e,2))*pow(M,2)*pow(r1,2) + 4*(32 -
    	125*pow(e,2))*M*pow(r1,3) - 156*pow(r1,4)) - 3*pow(r1,3)*(210*pow(M,4) + (-295 +
    	548*pow(e,2))*pow(M,3)*r1 + (313 - 488*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(-97 +
    	95*pow(e,2))*M*pow(r1,3) + 40*pow(r1,4))) + 3*pow(a,6)*(pow(l,2)*(8*(7 +
    	9*pow(e,2))*pow(M,2) + 8*(6 + pow(e,2))*M*r1 + 9*pow(r1,2)) +
    	r1*(2*pow(e,2)*M*(42*pow(M,2) + 148*M*r1 - 37*pow(r1,2)) + r1*(-240*pow(M,2) + 87*M*r1 +
    	11*pow(r1,2)))) + 3*pow(a,4)*(8*pow(l,4)*M*(3*M + r1) + pow(l,2)*r1*(12*(7 +
    	36*pow(e,2))*pow(M,3) + 4*(-23 + 98*pow(e,2))*pow(M,2)*r1 + 4*(3 -
    	19*pow(e,2))*M*pow(r1,2) - 13*pow(r1,3)) - pow(r1,2)*(r1*(-385*pow(M,3) + 317*pow(M,2)*r1 -
    	120*M*pow(r1,2) + 42*pow(r1,3)) + 2*pow(e,2)*M*(136*pow(M,3) + 144*pow(M,2)*r1 -
    	312*M*pow(r1,2) + 101*pow(r1,3)))))/(48.*r1*pow(pow(a,2) + r1*(-2*M + r1),4));
    A540 =
    	(pow(r1,2)*(-72*pow(a,5)*e*l*M + 9*pow(a,6)*(4*pow(e,2)*M + r1) +
    	12*pow(a,3)*e*l*M*r1*(-12*M + 11*r1) + 12*a*e*l*M*pow(r1,2)*(32*pow(M,2) - 42*M*r1 +
    	17*pow(r1,2)) + pow(a,4)*(36*pow(l,2)*M + r1*(18*pow(e,2)*M*(4*M - 3*r1) + (24*M -
    	13*r1)*r1)) + pow(r1,2)*(-6*pow(l,2)*M*(32*pow(M,2) - 38*M*r1 + 11*pow(r1,2)) +
    	pow(r1,2)*(156*pow(M,3) + (-193 + 132*pow(e,2))*pow(M,2)*r1 + 2*(32 -
    	45*pow(e,2))*M*pow(r1,2) + 2*pow(r1,3))) - pow(a,2)*r1*(pow(l,2)*(-72*pow(M,2) + 78*M*r1)
    	+ r1*(12*pow(e,2)*M*(16*pow(M,2) - 23*M*r1 + 15*pow(r1,2)) + r1*(183*pow(M,2) - 196*M*r1 +
    	62*pow(r1,2))))))/(48.*pow(pow(a,2) + r1*(-2*M + r1),4));
    A600 = pow(r1,6)/pow(pow(a,2) +
    	r1*(-2*M + r1),3);
    A602 = (-24*pow(a,5)*e*l*M*(5*M + r1) -
    	16*pow(a,3)*e*l*M*r1*(3*pow(l,2) - 8*pow(M,2) + 2*M*r1 + 3*pow(r1,2)) +
    	8*a*e*l*M*pow(r1,2)*(pow(l,2)*(4*M - 6*r1) + (11*M - 15*r1)*pow(r1,2)) +
    	3*pow(a,6)*(r1*(-6*M + r1) + 4*pow(e,2)*M*(5*M + 3*r1)) + pow(r1,2)*(8*pow(l,4)*M*(-2*M +
    	r1) + M*pow(r1,4)*(43*M + 4*(-8 + 5*pow(e,2))*r1) + pow(l,2)*pow(r1,2)*(3*pow(M,2) + 4*(3
    	+ 4*pow(e,2))*M*r1 - 12*pow(r1,2))) + pow(a,2)*r1*(24*pow(l,4)*M + pow(r1,2)*(-34*pow(M,3)
    	+ (39 - 100*pow(e,2))*pow(M,2)*r1 + 2*(-47 + 38*pow(e,2))*M*pow(r1,2) + 36*pow(r1,3)) +
    	pow(l,2)*(-64*pow(M,3) - 8*(-3 + 2*pow(e,2))*pow(M,2)*r1 + 10*(-3 +
    	4*pow(e,2))*M*pow(r1,2) + 48*pow(r1,3))) + pow(a,4)*(3*pow(l,2)*(20*pow(M,2) + 4*(-1 +
    	2*pow(e,2))*M*r1 + pow(r1,2)) + r1*(r1*(32*pow(M,2) - 14*M*r1 + 39*pow(r1,2)) +
    	pow(e,2)*(-64*pow(M,3) + 8*pow(M,2)*r1 + 92*M*pow(r1,2)))))/(24.*pow(pow(a,2) + r1*(-2*M +
    	r1),4));
    A620 = (pow(r1,3)*(-48*pow(a,3)*e*l*M + 16*a*e*l*M*(M - 3*r1)*r1 +
    	3*pow(a,4)*(8*pow(e,2)*M + 5*r1) + M*r1*(pow(l,2)*(-8*M + 4*r1) + pow(r1,2)*(43*M + 4*(-8
    	+ 5*pow(e,2))*r1)) + pow(a,2)*(24*pow(l,2)*M + 2*r1*(r1*(-31*M + 18*r1) +
    	pow(e,2)*(-4*pow(M,2) + 22*M*r1)))))/(24.*pow(pow(a,2) + r1*(-2*M + r1),4));
    A700 =
    	(pow(r1,5)*(-pow(a,2) + M*r1))/(2.*pow(pow(a,2) + r1*(-2*M + r1),4));
    A702 =
    	(-12*pow(a,7)*e*l*M*((2 + 6*pow(e,2))*M + 3*r1) + 3*pow(a,8)*(2*pow(e,2)*M*(2*M - 7*r1)
    	+ r1*(26*M + 3*r1)) - 4*a*e*l*M*pow(r1,3)*(pow(l,2)*(-52*pow(M,2) + 144*M*r1 -
    	69*pow(r1,2)) + pow(r1,2)*(-85*pow(M,2) + 2*(61 + 27*pow(e,2))*M*r1 - 51*pow(r1,2))) -
    	8*pow(a,5)*e*l*M*(27*pow(l,2)*M + r1*((35 + 12*pow(e,2))*pow(M,2) + 3*(8 +
    	17*pow(e,2))*M*r1 - 9*pow(r1,2))) + M*pow(r1,3)*(-4*pow(l,4)*(26*pow(M,2) - 31*M*r1 +
    	9*pow(r1,2)) + pow(l,2)*pow(r1,2)*(-25*pow(M,2) + 4*(3 + 44*pow(e,2))*M*r1 + 2*(5 -
    	54*pow(e,2))*pow(r1,2)) + pow(r1,4)*(91*pow(M,2) + 2*(-59 + 40*pow(e,2))*M*r1 + 2*(23 -
    	29*pow(e,2))*pow(r1,2))) - 4*pow(a,3)*e*l*M*r1*(pow(l,2)*(72*pow(M,2) + 172*M*r1 -
    	69*pow(r1,2)) + r1*(-80*pow(M,3) - 65*pow(M,2)*r1 + 2*(64 + 69*pow(e,2))*M*pow(r1,2) -
    	84*pow(r1,3))) + pow(a,2)*r1*(4*pow(l,4)*M*(24*pow(M,2) + 35*M*r1 - 31*pow(r1,2)) +
    	pow(l,2)*r1*(-160*pow(M,4) - 52*(-1 + 2*pow(e,2))*pow(M,3)*r1 + (27 +
    	1004*pow(e,2))*pow(M,2)*pow(r1,2) + 2*(51 - 130*pow(e,2))*M*pow(r1,3) - 84*pow(r1,4)) -
    	pow(r1,3)*(166*pow(M,4) + (-289 + 532*pow(e,2))*pow(M,3)*r1 + (371 -
    	472*pow(e,2))*pow(M,2)*pow(r1,2) + 16*(-17 + 14*pow(e,2))*M*pow(r1,3) + 72*pow(r1,4))) +
    	pow(a,6)*(3*pow(l,2)*((4 + 72*pow(e,2))*pow(M,2) + 26*M*r1 + 3*pow(r1,2)) +
    	r1*(4*pow(e,2)*M*(35*pow(M,2) + 88*M*r1 - 48*pow(r1,2)) + r1*(-322*pow(M,2) + 147*M*r1 -
    	9*pow(r1,2)))) + pow(a,4)*(72*pow(l,4)*pow(M,2) + pow(l,2)*r1*(4*(35 +
    	72*pow(e,2))*pow(M,3) + 4*(-40 + 239*pow(e,2))*pow(M,2)*r1 + (23 -
    	152*pow(e,2))*M*pow(r1,2) - 36*pow(r1,3)) - pow(r1,2)*(4*pow(e,2)*M*(40*pow(M,3) +
    	78*pow(M,2)*r1 - 195*M*pow(r1,2) + 79*pow(r1,3)) + r1*(-376*pow(M,3) + 365*pow(M,2)*r1 -
    	205*M*pow(r1,2) + 90*pow(r1,3)))))/(48.*r1*pow(pow(a,2) + r1*(-2*M + r1),5));
    A720 =
    	-(pow(r1,2)*(-48*pow(a,5)*e*l*M + 4*pow(a,3)*e*l*M*(68*M - 57*r1)*r1 +
    	3*pow(a,6)*(8*pow(e,2)*M + 3*r1) - 4*a*e*l*M*pow(r1,2)*(52*pow(M,2) - 98*M*r1 +
    	45*pow(r1,2)) + M*pow(r1,2)*(2*pow(l,2)*(52*pow(M,2) - 64*M*r1 + 19*pow(r1,2)) +
    	pow(r1,2)*(-91*pow(M,2) + 2*(59 - 40*pow(e,2))*M*r1 + 2*(-23 + 29*pow(e,2))*pow(r1,2))) +
    	pow(a,4)*(24*pow(l,2)*M + r1*(-2*pow(e,2)*M*(68*M - 65*r1) + r1*(-79*M + 42*r1))) +
    	pow(a,2)*r1*(2*pow(l,2)*M*(-68*M + 49*r1) + r1*(4*pow(e,2)*M*(26*pow(M,2) - 66*M*r1 +
    	41*pow(r1,2)) + r1*(187*pow(M,2) - 212*M*r1 + 72*pow(r1,2))))))/(48.*pow(pow(a,2) + r1*(-2*M
    	+ r1),5));
    A800 = (pow(r1,3)*(-12*pow(a,3)*e*l*M + 4*a*e*l*M*(2*M - 3*r1)*r1 +
    	3*pow(a,4)*(2*pow(e,2)*M + r1) + M*r1*(pow(l,2)*(-4*M + 2*r1) + pow(r1,2)*(5*M + 4*(-1 +
    	pow(e,2))*r1)) + 2*pow(a,2)*(3*pow(l,2)*M + r1*(r1*(-5*M + 3*r1) + pow(e,2)*M*(-2*M +
    	5*r1)))))/(12.*pow(pow(a,2) + r1*(-2*M + r1),5));
    A900 = -(pow(r1,2)*(-18*pow(a,5)*e*l*M +
    	2*pow(a,3)*e*l*M*(26*M - 21*r1)*r1 + 3*pow(a,6)*(3*pow(e,2)*M + r1) -
    	8*a*e*l*M*pow(r1,2)*(2*pow(M,2) - 6*M*r1 + 3*pow(r1,2)) +
    	M*pow(r1,2)*(pow(l,2)*(8*pow(M,2) - 10*M*r1 + 3*pow(r1,2)) + 2*pow(r1,2)*(-5*pow(M,2) + (7
    	- 4*pow(e,2))*M*r1 + 3*(-1 + pow(e,2))*pow(r1,2))) + pow(a,4)*(9*pow(l,2)*M +
    	r1*(-26*pow(e,2)*M*(M - r1) + r1*(-17*M + 9*r1))) + pow(a,2)*r1*(2*pow(l,2)*M*(-13*M +
    	8*r1) + r1*(r1*(30*pow(M,2) - 35*M*r1 + 12*pow(r1,2)) + pow(e,2)*M*(8*pow(M,2) - 38*M*r1 +
    	23*pow(r1,2))))))/(24.*pow(pow(a,2) + r1*(-2*M + r1),6));
  }

  /* Compute alpha, beta coefficients */
  alpha20 = r1*r1/(a*a+r1*(r1-2*M));
  alpha02 = r1*r1;
  beta = l*l + r1*r1 + a*a*(r1+2*M)/r1;
}
