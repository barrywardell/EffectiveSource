/*******************************************************************************
 * Copyright (C) 2012 Barry Wardell
 ******************************************************************************/

#include <math.h>
#include <assert.h>
#include "effsource.h"
#include <stdio.h>
#include <gsl/gsl_sf.h>

/* The particle's coordinate location and 4-velocity */
static struct coordinate xp;

/* Mass and spin of the Kerr black hole */
static double M, a;

/* Static variables used to store the coefficients of the series expansions */
static double A006, A008, A024, A026, A042, A044, A060, A062, A080, A106, A108, A124, A126, A142, A144, A160, A162, A180, A204, A206, A222, A224, A240, A242, A260, A304, A306, A322, A324, A340, A342, A360, A402, A404, A420, A422, A440, A502, A504, A520, A522, A540, A600, A602, A620, A700, A702, A720, A800, A900;
static double alpha20, alpha02, beta;

/* Legendre function for z>1, n half-integer, m integer */
double LegendreP(double l, int m, double z)
{
  double n = l < 0 ? -(l+1) : l;
  double a = -n;
  double b = 1+n;
  double h;

  if(m>=0) {
    double c = 1.0+m;
    h = gsl_sf_gamma(-a + m + 1)/gsl_sf_gamma(-a - m + 1)*pow(z-1.0,m/2.0)*pow(z+1.0,-m/2.0)*pow((1.0+z)/2.0,-b)*gsl_sf_hyperg_2F1_renorm(c-a, b, c, (z-1.0)/(z+1.0));
  } else {
    double c = 1.0-m;
    h = pow(z-1.0,-m/2.0)*pow(z+1.0,m/2.0)*pow((1.0+z)/2.0,-b)*gsl_sf_hyperg_2F1_renorm(c-a, b, c, (z-1.0)/(z+1.0));
  }
  return h;
}

/* Compute the singular field at the point x for the particle at xp */
void effsource_PhiS(struct coordinate * x, double * PhiS)
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
  
  double sindphi  = sin(0.5*dphi);
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

  *PhiS = A/pow(rho2, 3.5);
}

/* Compute the singular field at the point x for the particle at xp */
void effsource_PhiS_m(int m, struct coordinate * x, double * PhiS)
{
  double A[5], alpha;

  double r      = x->r;
  const double theta  = x->theta;
  const double rp     = xp.r;
  const double thetap = xp.theta;

  const double dr     = r - rp;
  const double dtheta = theta - thetap;

  const double dr2 = dr*dr;
  const double dr3 = dr2*dr;
  const double dr4 = dr2*dr2;
  const double dr6 = dr3*dr3;
  const double dr8 = dr4*dr4;

  const double dtheta2  = dtheta*dtheta;
  const double dtheta4  = dtheta2*dtheta2;
  const double dtheta8  = dtheta4*dtheta4;

  A[0] = dr6*(A600 + A700*dr) + dr8*(A800 + A900*dr) + dr4*(A420 + A520*dr + dr2*(A620 + A720*dr))*dtheta2 + (A080 + A180*dr)*dtheta8 + dtheta4*(dr2*(A240 + A340*dr) + dr4*(A440 + A540*dr) + (A060 + A160*dr + dr2*(A260 + A360*dr))*dtheta2);
  A[1] = dr4*(A402 + A502*dr + dr2*(A602 + A702*dr)) + (dr2*(A222 + A322*dr) + dr4*(A422 + A522*dr))*dtheta2 + dtheta4*(A042 + A142*dr + dr2*(A242 + A342*dr) + (A062 + A162*dr)*dtheta2);
  A[2] = dr2*(A204 + A304*dr) + dr4*(A404 + A504*dr) + (A024 + A124*dr + dr2*(A224 + A324*dr))*dtheta2 + (A044 + A144*dr)*dtheta4;
  A[3] = A006 + A106*dr + dr2*(A206 + A306*dr) + (A026 + A126*dr)*dtheta2;
  A[4] = A008 + A108*dr;
  
  alpha = alpha20*dr2 + alpha02*dtheta2;

  const double zc = sqrt(beta);
  const double rho = sqrt(alpha);
  const double R0 = sqrt(pow(rho,2) + pow(zc,2));
  const double z = (pow(rho,2) + pow(zc,2)/2.0)/(rho*R0);

  double R[10], ZpR[10];
  int mBlock, mp1, mm1;
  for(int i=0; i<=8; i+=2)
  {
    int n = i-7;
    if (1+n/2.0 <= 0) {
      mBlock = -abs(m);
      mm1 = -abs(m-1);
      mp1 = -abs(m+1);
    } else {
      mBlock = abs(m);
      mm1 = abs(m-1);
      mp1 = abs(m+1);
    }
    const double pch = gsl_sf_poch(1+n/2.0, mBlock);
    const double pchmp1 = gsl_sf_poch(1+n/2.0, mp1);
    const double pchmm1 = gsl_sf_poch(1+n/2.0, mm1);
    R[i]   = 2*M_PI*pow(-1,mBlock)*pow(rho,n/2.0)*pow(R0,n/2.0)*LegendreP(n/2.0, mBlock, z)/pch;
    ZpR[i] = pow(-1,m)*M_PI*zc*pow(rho,n/2.0)*pow(R0,n/2.0)*(LegendreP(n/2.0, mm1, z)/pchmm1-LegendreP(n/2.0, mp1, z)/pchmp1);
  }

  double RePhiS = 1.0/pow(beta,4)*(
	A[4]*R[8]
	+ (beta*A[3]-4*alpha*A[4])*R[6]
	+ (pow(beta,2)*A[2]-3*alpha*beta*A[3]+6*pow(alpha,2)*A[4])*R[4]
	+ (pow(beta,3)*A[1]-2*alpha*pow(beta,2)*A[2]+3*pow(alpha,2)*beta*A[3]-4*pow(alpha,3)*A[4])*R[2]
	+ (pow(beta,4)*A[0]-alpha*pow(beta,3)*A[1]+pow(alpha,2)*pow(beta,2)*A[2]-pow(alpha,3)*beta*A[3]+pow(alpha,4)*A[4])*R[0]);

  /* Store calculated quantities into the arrays provided by the caller,
     including the phase factor exp(-i*m*phi_p) */
  double cosmph = cos(m*xp.phi);
  double sinmph = sin(m*xp.phi);

  PhiS[0] = RePhiS*cosmph;
  PhiS[1] = - RePhiS*sinmph;
}

/* Compute the singular field, its derivatives and its d'Alembertian */
void effsource_calc(struct coordinate * x,
  double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src)
{
  double A, dA_dr, d2A_dr2, dA_dth, d2A_dth2, dA_dR, dA_dph,  d2A_dR2,  d2A_dph2, dA_dt, d2A_dt2, d2A_dphdt;
  double s2, sqrts2, s2_15, s2_25, s2_35, s2_45, s2_55, ds2_dr, d2s2_dr2, ds2_dth, d2s2_dth2, ds2_dR, ds2_dph, d2s2_dR2, d2s2_dph2, ds2_dt, d2s2_dt2, d2s2_dphdt;

  double dPhiS_dt, dPhiS_dr, dPhiS_dth, dPhiS_dph, d2PhiS_dt2, d2PhiS_dtr, d2PhiS_dtth;
  double d2PhiS_dtph, d2PhiS_dr2, d2PhiS_drth, d2PhiS_drph, d2PhiS_dth2, d2PhiS_dthph, d2PhiS_dph2;

  double r      = x->r;
  double theta  = x->theta;
  double phi    = x->phi;
  double rp     = xp.r;
  double thetap = xp.theta;
  double phip   = xp.phi;

  double dr     = r - rp;
  double dtheta = theta - thetap;
  double dphi   = phi - phip;

  double dr2      = dr*dr;
  double dr3      = dr2*dr;
  double dr4      = dr2*dr2;
  double dr5      = dr3*dr2;
  double dr6      = dr3*dr3;
  double dr7      = dr4*dr3;
  double dr8      = dr4*dr4;

  double dtheta2  = dtheta*dtheta;
  double dtheta3  = dtheta2*dtheta;
  double dtheta4  = dtheta2*dtheta2;
  double dtheta5  = dtheta3*dtheta2;
  double dtheta6  = dtheta3*dtheta3;
  double dtheta7  = dtheta4*dtheta3;
  double dtheta8  = dtheta4*dtheta4;

  double R        = sin(0.5*dphi);
  double R2       = R*R;
  double R3       = R2*R;
  double R4       = R2*R2;
  double R5       = R3*R2;
  double R6       = R3*R3;
  double R7       = R4*R3;
  double R8       = R4*R4;
  double dR       = 0.5*cos(0.5*dphi);

  double om       = M / (a*M + sqrt(M*pow(rp,3)));

  /* A, dA/dx, d^2A/dx^2 */
  A         = dr6*(A600 + A700*dr) + dr8*(A800 + A900*dr) + dr4*(A420 + A520*dr + dr2*(A620 + A720*dr))*dtheta2 + (A080 + A180*dr)*dtheta8 + dtheta4*(dr2*(A240 + A340*dr) + dr4*(A440 + A540*dr) + (A060 + A160*dr + dr2*(A260 + A360*dr))*dtheta2) + (dr4*(A402 + A502*dr + dr2*(A602 + A702*dr)) + (dr2*(A222 + A322*dr) + dr4*(A422 + A522*dr))*dtheta2 + dtheta4*(A042 + A142*dr + dr2*(A242 + A342*dr) + (A062 + A162*dr)*dtheta2))*R2 + (A008 + A108*dr)*R8 + R4*(dr2*(A204 + A304*dr) + dr4*(A404 + A504*dr) + (A024 + A124*dr + dr2*(A224 + A324*dr))*dtheta2 + (A044 + A144*dr)*dtheta4 + (A006 + A106*dr + dr2*(A206 + A306*dr) + (A026 + A126*dr)*dtheta2)*R2);;
  dA_dr     = 6*(A600 + A700*dr)*dr5 + A700*dr6 + 8*(A800 + A900*dr)*dr7 + A900*dr8 + 4*(A420 + A520*dr + (A620 + A720*dr)*dr2)*dr3*dtheta2 + (A520 + 2*dr*(A620 + A720*dr) + A720*dr2)*dr4*dtheta2 + (2*dr*(A240 + A340*dr) + A340*dr2 + 4*(A440 + A540*dr)*dr3 + A540*dr4 + (A160 + 2*dr*(A260 + A360*dr) + A360*dr2)*dtheta2)*dtheta4 + A180*dtheta8 + (4*(A402 + A502*dr + (A602 + A702*dr)*dr2)*dr3 + (A502 + 2*dr*(A602 + A702*dr) + A702*dr2)*dr4 + (2*dr*(A222 + A322*dr) + A322*dr2 + 4*(A422 + A522*dr)*dr3 + A522*dr4)* dtheta2 + (A142 + 2*dr*(A242 + A342*dr) + A342*dr2 + A162*dtheta2)* dtheta4)*R2 + (2*dr*(A204 + A304*dr) + A304*dr2 + 4*(A404 + A504*dr)*dr3 + A504*dr4 + (A124 + 2*dr*(A224 + A324*dr) + A324*dr2)*dtheta2 + A144*dtheta4 + (A106 + 2*dr*(A206 + A306*dr) + A306*dr2 + A126*dtheta2)*R2)* R4 + A108*R8;
  dA_dth    = 2*(A420 + A520*dr + (A620 + A720*dr)*dr2)*dr4*dtheta + 4*((A240 + A340*dr)*dr2 + (A440 + A540*dr)*dr4 + (A060 + A160*dr + (A260 + A360*dr)*dr2)*dtheta2)*dtheta3 + 2*(A060 + A160*dr + (A260 + A360*dr)*dr2)*dtheta5 + 8*(A080 + A180*dr)*dtheta7 + (2*((A222 + A322*dr)*dr2 + (A422 + A522*dr)*dr4)*dtheta + 4*(A042 + A142*dr + (A242 + A342*dr)*dr2 + (A062 + A162*dr)*dtheta2)* dtheta3 + 2*(A062 + A162*dr)*dtheta5)*R2 + (2*(A024 + A124*dr + (A224 + A324*dr)*dr2)*dtheta + 4*(A044 + A144*dr)*dtheta3 + 2*(A026 + A126*dr)*dtheta*R2)* R4;
  dA_dR     = 2*((A402 + A502*dr + (A602 + A702*dr)*dr2)*dr4 + ((A222 + A322*dr)*dr2 + (A422 + A522*dr)*dr4)*dtheta2 + (A042 + A142*dr + (A242 + A342*dr)*dr2 + (A062 + A162*dr)*dtheta2)* dtheta4)*R + 4*((A204 + A304*dr)*dr2 + (A404 + A504*dr)*dr4 + (A024 + A124*dr + (A224 + A324*dr)*dr2)*dtheta2 + (A044 + A144*dr)*dtheta4 + (A006 + A106*dr + (A206 + A306*dr)*dr2 + (A026 + A126*dr)*dtheta2)* R2)*R3 + 2*(A006 + A106*dr + (A206 + A306*dr)*dr2 + (A026 + A126*dr)*dtheta2)*R5 + 8*(A008 + A108*dr)*R7;
  dA_dph    = dA_dR*dR;
  dA_dt     = -om*dA_dph;
  d2A_dr2   = 30*(A600 + A700*dr)*dr4 + 12*A700*dr5 + 56*(A800 + A900*dr)*dr6 + 16*A900*dr7 + 12*dr2*(A420 + A520*dr + (A620 + A720*dr)*dr2)*dtheta2 + 8*(A520 + 2*dr*(A620 + A720*dr) + A720*dr2)*dr3*dtheta2 + (4*A720*dr + 2*(A620 + A720*dr))*dr4*dtheta2 + (4*A340*dr + 2*(A240 + A340*dr) + 12*(A440 + A540*dr)*dr2 + 8*A540*dr3 + (4*A360*dr + 2*(A260 + A360*dr))*dtheta2)*dtheta4 + (12*dr2*(A402 + A502*dr + (A602 + A702*dr)*dr2) + 8*(A502 + 2*dr*(A602 + A702*dr) + A702*dr2)*dr3 + (4*A702*dr + 2*(A602 + A702*dr))*dr4 + (4*A322*dr + 2*(A222 + A322*dr) + 12*(A422 + A522*dr)*dr2 + 8*A522*dr3)*dtheta2 + (4*A342*dr + 2*(A242 + A342*dr))*dtheta4)* R2 + (4*A304*dr + 2*(A204 + A304*dr) + 12*(A404 + A504*dr)*dr2 + 8*A504*dr3 + (4*A324*dr + 2*(A224 + A324*dr))*dtheta2 + (4*A306*dr + 2*(A206 + A306*dr))*R2)*R4;
  d2A_dth2  = 2*(A420 + A520*dr + (A620 + A720*dr)*dr2)*dr4 + 12*dtheta2*((A240 + A340*dr)*dr2 + (A440 + A540*dr)*dr4 + (A060 + A160*dr + (A260 + A360*dr)*dr2)*dtheta2) + 18*(A060 + A160*dr + (A260 + A360*dr)*dr2)*dtheta4 + 56*(A080 + A180*dr)*dtheta6 + (2*((A222 + A322*dr)*dr2 + (A422 + A522*dr)*dr4) + 12*dtheta2*(A042 + A142*dr + (A242 + A342*dr)*dr2 +  (A062 + A162*dr)*dtheta2) + 18*(A062 + A162*dr)*dtheta4)*R2 + (2*(A024 + A124*dr + (A224 + A324*dr)*dr2) + 12*(A044 + A144*dr)*dtheta2 + 2*(A026 + A126*dr)*R2)*R4; d2A_dR2 = 2*((A402 + A502*dr + (A602 + A702*dr)*dr2)*dr4 + ((A222 + A322*dr)*dr2 + (A422 + A522*dr)*dr4)*dtheta2 + (A042 + A142*dr + (A242 + A342*dr)*dr2 + (A062 + A162*dr)*dtheta2)* dtheta4) + 12*R2*((A204 + A304*dr)*dr2 + (A404 + A504*dr)*dr4 + (A024 + A124*dr + (A224 + A324*dr)*dr2)*dtheta2 + (A044 + A144*dr)*dtheta4 + (A006 + A106*dr + (A206 + A306*dr)*dr2 + (A026 + A126*dr)*dtheta2)* R2) + 18*(A006 + A106*dr + (A206 + A306*dr)*dr2 + (A026 + A126*dr)*dtheta2)*R4 + 56*(A008 + A108*dr)*R6;
  d2A_dph2  = - 0.25*R*dA_dR + dR*dR*d2A_dR2;
  d2A_dt2   = om*om*d2A_dph2;
  d2A_dphdt = -om*d2A_dph2;

  /* s, ds/dr, d^2s/dr^2 */
  s2         = alpha20*dr2 + alpha02*dtheta2 + beta*R2;
  ds2_dr     = 2*alpha20*dr;
  ds2_dth    = 2*alpha02*dtheta;
  ds2_dR     = 2*beta*R;
  ds2_dph    = ds2_dR*dR;
  ds2_dt     = -om*ds2_dph;
  d2s2_dr2   = 2*alpha20;
  d2s2_dth2  = 2*alpha02;
  d2s2_dR2   = 2*beta;
  d2s2_dph2  = - 0.25*R*ds2_dR + dR*dR*d2s2_dR2;
  d2s2_dt2   = om*om*d2s2_dph2;
  d2s2_dphdt = -om*d2s2_dph2;
  sqrts2     = sqrt(s2);
  s2_15      = s2*sqrts2;
  s2_25      = s2*s2_15;
  s2_35      = s2*s2_25;
  s2_45      = s2*s2_35;
  s2_55      = s2*s2_45;

  /* PhiS */
  *PhiS = A/s2_35;

  /* First derivatives of PhiS */
  dPhiS_dt  = (-7*ds2_dt*A + 2*dA_dt*s2)/(2.*s2_45);
  dPhiS_dr  = (-7*ds2_dr*A + 2*dA_dr*s2) /(2.*s2_45);
  dPhiS_dth = (-7*ds2_dth*A + 2*dA_dth*s2)/(2.*s2_45);
  dPhiS_dph = (-7*ds2_dph*A + 2*dA_dph*s2)/(2.*s2_45);
  
  /* Second derivatives of PhiS */
  d2PhiS_dr2 =
    (63*ds2_dr*ds2_dr*A - 14*s2*(2*dA_dr*ds2_dr + d2s2_dr2*A) + 4*d2A_dr2*s2*s2)/(4.*s2_55);
  d2PhiS_dth2 =
    (63*ds2_dth*ds2_dth*A - 14*s2*(2*dA_dth*ds2_dth + d2s2_dth2*A) + 4*d2A_dth2*s2*s2)/(4.*s2_55);
  d2PhiS_dph2 =
    (63*ds2_dph*ds2_dph*A - 14*s2*(2*dA_dph*ds2_dph + d2s2_dph2*A) + 4*d2A_dph2*s2*s2)/(4.*s2_55);
  d2PhiS_dt2 =
    (63*ds2_dt*ds2_dt*A - 14*s2*(2*dA_dt*ds2_dt + d2s2_dt2*A) + 4*d2A_dt2*s2*s2)/(4.*s2_55);
  d2PhiS_dtph =
    (63*ds2_dph*ds2_dt*A - 14*s2*(dA_dt*ds2_dph + dA_dph*ds2_dt + d2s2_dphdt*A) + 4*d2A_dphdt*s2*s2)/(4.*s2_55);
  d2PhiS_dtr   = NAN;
  d2PhiS_dtth  = NAN;
  d2PhiS_drth  = NAN;
  d2PhiS_drph  = NAN;
  d2PhiS_dthph = NAN;
  
  
  /* Box[PhiS] */
  double sinth  = sin(theta);
  double sinth2 = sinth*sinth;
  double sin2th = sin(2*theta);
  double cos2th = cos(2*theta);
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double a2 = a*a;
  double a4 = a2*a2;

  *src = -((2*a2*dPhiS_dr - a2*d2PhiS_dph2 - a4*d2PhiS_dr2 - a2*d2PhiS_dth2 - 4*dPhiS_dr*r - 2*a2*dPhiS_dr*r +
         4*d2PhiS_dph2*r + 2*a*d2PhiS_dtph*r + 4*a2*d2PhiS_dr2*r + 2*d2PhiS_dth2*r + 6*dPhiS_dr*r2 - 2*d2PhiS_dph2*r2 - 4*d2PhiS_dr2*r2 -
         2*a2*d2PhiS_dr2*r2 - d2PhiS_dth2*r2 - 2*dPhiS_dr*r3 + 4*d2PhiS_dr2*r3 - d2PhiS_dr2*r4 +
         (a4*d2PhiS_dt2 + 4*a*d2PhiS_dtph*r + 2*d2PhiS_dt2*r4 + a2*d2PhiS_dt2*r*(2 + 3*r))*sinth2 +
         cos2th*(a4*d2PhiS_dr2 - 2*a*d2PhiS_dtph*r + (-2 + r)*r*(d2PhiS_dth2 + 2*dPhiS_dr*(-1 + r) + d2PhiS_dr2*(-2 + r)*r) +
            a2*(-d2PhiS_dph2 + d2PhiS_dth2 + 2*dPhiS_dr*(-1 + r) - 4*d2PhiS_dr2*r + 2*d2PhiS_dr2*r2) +
            a2*d2PhiS_dt2*(a2 + (-2 + r)*r)*sinth2) - a2*dPhiS_dth*sin2th + 2*dPhiS_dth*r*sin2th -
         dPhiS_dth*r2*sin2th))/((sinth2*(a2 + (-2 + r)*r)*(a2 + 2*r2 + a2*cos2th)));

  dPhiS_dx[0] = dPhiS_dt;
  dPhiS_dx[1] = dPhiS_dr;
  dPhiS_dx[2] = dPhiS_dth;
  dPhiS_dx[3] = dPhiS_dph;

  d2PhiS_dx2[0] = d2PhiS_dt2;
  d2PhiS_dx2[1] = d2PhiS_dtr;
  d2PhiS_dx2[2] = d2PhiS_dtth;
  d2PhiS_dx2[3] = d2PhiS_dtph;
  d2PhiS_dx2[4] = d2PhiS_dr2;
  d2PhiS_dx2[5] = d2PhiS_drth;
  d2PhiS_dx2[6] = d2PhiS_drph;
  d2PhiS_dx2[7] = d2PhiS_dth2;
  d2PhiS_dx2[8] = d2PhiS_dthph;
  d2PhiS_dx2[9] = d2PhiS_dph2;
}

/* Compute the 2D singular field, its derivatives and its d'Alembertian */
void effsource_calc_m(int m, struct coordinate * x,
  double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src)
{
  double A[5], alpha;
  double dA_dr[5], dalpha_dr;
  double d2A_dr2[5], d2alpha_dr2;
  double dA_dtheta[5], dalpha_dtheta;
  double d2A_dtheta2[5], d2alpha_dtheta2;

  double s, ds_dr, d2s_dr2, ds_dtheta, d2s_dtheta2;

  double dPhiS_dt, dPhiS_dr, dPhiS_dth, dPhiS_dph, d2PhiS_dt2, d2PhiS_dtr, d2PhiS_dtth;
  double d2PhiS_dtph, d2PhiS_dr2, d2PhiS_drth, d2PhiS_drph, d2PhiS_dth2, d2PhiS_dthph, d2PhiS_dph2;

  const double r      = x->r;
  const double theta  = x->theta;
  const double rp     = xp.r;
  const double thetap = xp.theta;
  const double om       = M / (a*M + sqrt(M*pow(rp,3)));

  const double dr     = r - rp;
  const double dtheta = theta - thetap;

  const double dr2 = dr*dr;
  const double dr3 = dr2*dr;
  const double dr4 = dr2*dr2;
  const double dr5 = dr3*dr2;
  const double dr6 = dr3*dr3;
  const double dr7 = dr4*dr3;
  const double dr8 = dr4*dr4;

  const double dtheta2  = dtheta*dtheta;
  const double dtheta3  = dtheta2*dtheta;
  const double dtheta4  = dtheta2*dtheta2;
  const double dtheta5  = dtheta3*dtheta2;
  const double dtheta6  = dtheta3*dtheta3;
  const double dtheta7  = dtheta4*dtheta3;
  const double dtheta8  = dtheta4*dtheta4;

  /* Coefficients of sin(dphi) in the numerator */
  A[0] = dr6*(A600 + A700*dr) + dr8*(A800 + A900*dr) + dr4*(A420 + A520*dr + dr2*(A620 + A720*dr))*dtheta2 + (A080 + A180*dr)*dtheta8 + dtheta4*(dr2*(A240 + A340*dr) + dr4*(A440 + A540*dr) + (A060 + A160*dr + dr2*(A260 + A360*dr))*dtheta2);
  A[1] = dr4*(A402 + A502*dr + dr2*(A602 + A702*dr)) + (dr2*(A222 + A322*dr) + dr4*(A422 + A522*dr))*dtheta2 + dtheta4*(A042 + A142*dr + dr2*(A242 + A342*dr) + (A062 + A162*dr)*dtheta2);
  A[2] = dr2*(A204 + A304*dr) + dr4*(A404 + A504*dr) + (A024 + A124*dr + dr2*(A224 + A324*dr))*dtheta2 + (A044 + A144*dr)*dtheta4;
  A[3] = A006 + A106*dr + dr2*(A206 + A306*dr) + (A026 + A126*dr)*dtheta2;
  A[4] = A008 + A108*dr;

  /* r derivatives of coefficients */
  dA_dr[0] = 6*A600*dr5 + 7*A700*dr6 + 8*A800*dr7 + 9*A900*dr8 + 4*A420*dr3*dtheta2 + 5*A520*dr4*dtheta2 + 6*A620*dr5*dtheta2 + 7*A720*dr6*dtheta2 + 2*A240*dr*dtheta4 + 3*A340*dr2*dtheta4 + 4*A440*dr3*dtheta4 + 5*A540*dr4*dtheta4 + A160*dtheta6 + 2*A260*dr*dtheta6 + 3*A360*dr2*dtheta6 + A180*dtheta8;
  dA_dr[1] = 4*A402*dr3 + 5*A502*dr4 + 6*A602*dr5 + 7*A702*dr6 + 2*A222*dr*dtheta2 + 3*A322*dr2*dtheta2 + 4*A422*dr3*dtheta2 + 5*A522*dr4*dtheta2 + A142*dtheta4 + 2*A242*dr*dtheta4 + 3*A342*dr2*dtheta4 + A162*dtheta6;
  dA_dr[2] = 2*A204*dr + 3*A304*dr2 + 4*A404*dr3 + 5*A504*dr4 + A124*dtheta2 + 2*A224*dr*dtheta2 + 3*A324*dr2*dtheta2 + A144*dtheta4;
  dA_dr[3] = A106 + 2*A206*dr + 3*A306*dr2 + A126*dtheta2;
  dA_dr[4] = A108;

  /* r,r derivatives of coefficients */
  d2A_dr2[0] = 2*(15*A600*dr4 + 21*A700*dr5 + 28*A800*dr6 + 36*A900*dr7 + 6*A420*dr2*dtheta2 + 10*A520*dr3*dtheta2 + 15*A620*dr4*dtheta2 + 21*A720*dr5*dtheta2 + A240*dtheta4 + 3*A340*dr*dtheta4 + 6*A440*dr2*dtheta4 + 10*A540*dr3*dtheta4 + A260*dtheta6 + 3*A360*dr*dtheta6);
  d2A_dr2[1] = 2*(6*A402*dr2 + 10*A502*dr3 + 15*A602*dr4 + 21*A702*dr5 + A222*dtheta2 + 3*A322*dr*dtheta2 + 6*A422*dr2*dtheta2 + 10*A522*dr3*dtheta2 + A242*dtheta4 + 3*A342*dr*dtheta4);
  d2A_dr2[2] = 2*(A204 + 3*A304*dr + 6*A404*dr2 + 10*A504*dr3 + A224*dtheta2 + 3*A324*dr*dtheta2);
  d2A_dr2[3] = 2*(A206 + 3*A306*dr);
  d2A_dr2[4] = 0;

  /* theta derivatives of coefficients */
  dA_dtheta[0] = 2*(A420 + dr*(A520 + dr*(A620 + A720*dr)))*dr4*dtheta + 4*((A240 + A340*dr)*dr2 + (A440 + A540*dr)*dr4 + (A060 + dr*(A160 + dr*(A260 + A360*dr)))*dtheta2)*dtheta3 + 2*(A060 + dr*(A160 + dr*(A260 + A360*dr)))*dtheta5 + 8*(A080 + A180*dr)*dtheta7;
  dA_dtheta[1] = 2*dtheta*((A222 + A322*dr)*dr2 + (A422 + A522*dr)*dr4 + 2*dtheta2*(A042 + A142*dr + (A242 + A342*dr)*dr2 + (A062 + A162*dr)*dtheta2) + (A062 + A162*dr)*dtheta4);
  dA_dtheta[2] = 2*dtheta*(A024 + A124*dr + (A224 + A324*dr)*dr2 + 2*(A044 + A144*dr)*dtheta2);
  dA_dtheta[3] = 2*(A026 + A126*dr)*dtheta;
  dA_dtheta[4] = 0;

  /* theta,theta derivatives of coefficients */
  d2A_dtheta2[0] = 2*((A420 + dr*(A520 + dr*(A620 + A720*dr)))*dr4 + 6*dtheta2*((A240 + A340*dr)*dr2 + (A440 + A540*dr)*dr4 + (A060 + dr*(A160 + dr*(A260 + A360*dr)))*dtheta2) + 9*(A060 + dr*(A160 + dr*(A260 + A360*dr)))*dtheta4 + 28*(A080 + A180*dr)*dtheta6);
  d2A_dtheta2[1] = 2*((A222 + A322*dr)*dr2 + (A422 + A522*dr)*dr4 + 6*dtheta2*(A042 + A142*dr + (A242 + A342*dr)*dr2 + (A062 + A162*dr)*dtheta2) + 9*(A062 + A162*dr)*dtheta4);
  d2A_dtheta2[2] = 2*(A024 + A124*dr + (A224 + A324*dr)*dr2 + 6*(A044 + A144*dr)*dtheta2);
  d2A_dtheta2[3] = 2*(A026 + A126*dr);
  d2A_dtheta2[4] = 0;

  /* alpha term appearing in rho */
  alpha = alpha20*dr2 + alpha02*dtheta2;

  /* Derivatives of alpha */
  dalpha_dr       = 2*alpha20*dr;
  d2alpha_dr2     = 2*alpha20;
  dalpha_dtheta   = 2*alpha02*dtheta;
  d2alpha_dtheta2 = 2*alpha02;

  const double zc = sqrt(beta);
  const double rho = sqrt(alpha);
  const double R0 = sqrt(pow(rho,2) + pow(zc,2));
  const double z = (pow(rho,2) + pow(zc,2)/2.0)/(rho*R0);

  const double drho_dr = 0.5*dalpha_dr/rho;
  const double d2rho_dr2 = 0.5*d2alpha_dr2/rho-0.25*pow(dalpha_dr,2)/pow(rho,3);
  const double drho_dtheta = 0.5*dalpha_dtheta/rho;
  const double d2rho_dtheta2 = 0.5*d2alpha_dtheta2/rho-0.25*pow(dalpha_dtheta,2)/pow(rho,3);

  /* m-mode integrals of R */
  double R[9], DrhoR[9], D2rhoR[9];
  double ZpR[9], DrhoZpR[9], D2rhoZpR[9];
  int mBlock, mp1, mm1;
  for(int i=0; i<=8; i+=2)
  {
    int n = i-7;
    if (1+n/2.0 <= 0) {
      mBlock = -abs(m);
      mm1 = -abs(m-1);
      mp1 = -abs(m+1);
    } else {
      mBlock = abs(m);
      mm1 = abs(m-1);
      mp1 = abs(m+1);
    }
    const double pch = gsl_sf_poch(1+n/2.0, mBlock);
    const double pchmp1 = gsl_sf_poch(1+n/2.0, mp1);
    const double pchmm1 = gsl_sf_poch(1+n/2.0, mm1);
    R[i]   = 2*M_PI*pow(-1,mBlock)*pow(rho,n/2.0)*pow(R0,n/2.0)*LegendreP(n/2.0, mBlock, z)/pch;
    ZpR[i] = pow(-1,m)*M_PI*zc*pow(rho,n/2.0)*pow(R0,n/2.0)*(LegendreP(n/2.0, mm1, z)/pchmm1-LegendreP(n/2.0, mp1, z)/pchmp1);

    DrhoR[i]   = 2*pow(-1,mBlock)*M_PI*((-2 + 2*mBlock - n)*pow(R0,(-2 + n)/2.)*pow(rho,n/2.)*LegendreP(1 + n/2.,mBlock,z) + (1 + n)*pow(R0,2*(-1 + n/4.))*pow(rho,-1 + n/2.)*(2*pow(rho,2) + pow(zc,2))*LegendreP(n/2.,mBlock,z))/pch;
    D2rhoR[i]  = (2*pow(-1,mBlock)*M_PI*((-2 + 2*mBlock - n)*pow(R0,(-6 + n)/2.)*pow(rho,-1 + n/2.)*((5 + 4*n)*pow(rho,2) + (3 + 2*n)*pow(zc,2))*LegendreP(1 + n/2.,mBlock,z) + 
       (-4 + 2*mBlock - n)*(-2 + 2*mBlock - n)*pow(R0,2*(-1 + n/4.))*pow(rho,n/2.)*LegendreP(2 + n/2.,mBlock,z) + 
       (1 + n)*pow(R0,2*(-2 + n/4.))*pow(rho,-2 + n/2.)*(2*pow(rho,4) + 3*pow(rho,2)*pow(zc,2) + n*pow(2*pow(rho,2) + pow(zc,2),2))*LegendreP(n/2.,mBlock,z)))/pch;

    DrhoZpR[i] = pow(-1,m)*M_PI*zc*(((-2 + 2*mm1 - n)*pow(R0,(-2 + n)/2.)*pow(rho,n/2.)*LegendreP(1 + n/2.,mm1,z))/pchmm1
      + ((2 - 2*mp1 + n)*pow(R0,(-2 + n)/2.)*pow(rho,n/2.)*LegendreP(1 + n/2.,mp1,z))/pchmp1
      + ((1 + n)*pow(R0,2*(-1 + n/4.))*pow(rho,-1 + n/2.)*(2*pow(rho,2) + pow(zc,2))*LegendreP(n/2.,mm1,z))/pchmm1
      - ((1 + n)*pow(R0,2*(-1 + n/4.))*pow(rho,-1 + n/2.)*(2*pow(rho,2) + pow(zc,2))*LegendreP(n/2.,mp1,z))/pchmp1);
    D2rhoZpR[i] = pow(-1,m)*M_PI*zc*(((-2 + 2*mm1 - n)*pow(R0,(-6 + n)/2.)*pow(rho,-1 + n/2.)*((3 + 2*n)*pow(zc,2) + (5 + 4*n)*pow(rho,2))*LegendreP(1 + n/2.,mm1,z))/pchmm1
      - ((-2 + 2*mp1 - n)*pow(R0,(-6 + n)/2.)*pow(rho,-1 + n/2.)*((3 + 2*n)*pow(zc,2) + (5 + 4*n)*pow(rho,2))*LegendreP(1 + n/2.,mp1,z))/pchmp1
      + ((-4 + 2*mm1 - n)*(-2 + 2*mm1 - n)*pow(R0,2*(-1 + n/4.))*pow(rho,n/2.)*LegendreP(2 + n/2.,mm1,z))/pchmm1
      - ((-4 + 2*mp1 - n)*(-2 + 2*mp1 - n)*pow(R0,2*(-1 + n/4.))*pow(rho,n/2.)*LegendreP(2 + n/2.,mp1,z))/pchmp1
      + ((1 + n)*pow(R0,2*(-2 + n/4.))*pow(rho,-2 + n/2.)*(3*pow(zc,2)*pow(rho,2) + 2*pow(rho,4) + n*pow(pow(zc,2) + 2*pow(rho,2),2))*LegendreP(n/2.,mm1,z))/pchmm1
      - ((1 + n)*pow(R0,2*(-2 + n/4.))*pow(rho,-2 + n/2.)*(3*pow(zc,2)*pow(rho,2) + 2*pow(rho,4) + n*pow(pow(zc,2) + 2*pow(rho,2),2))*LegendreP(n/2.,mp1,z))/pchmp1);
  }

  /* Singular field */
  double RePhiS = 1.0/pow(beta,4)*(
	A[4]*R[8]
	+ (beta*A[3]-4*alpha*A[4])*R[6]
	+ (pow(beta,2)*A[2]-3*alpha*beta*A[3]+6*pow(alpha,2)*A[4])*R[4]
	+ (pow(beta,3)*A[1]-2*alpha*pow(beta,2)*A[2]+3*pow(alpha,2)*beta*A[3]-4*pow(alpha,3)*A[4])*R[2]
	+ (pow(beta,4)*A[0]-alpha*pow(beta,3)*A[1]+pow(alpha,2)*pow(beta,2)*A[2]-pow(alpha,3)*beta*A[3]+pow(alpha,4)*A[4])*R[0]);

  /* Derivatives of singular field - first account for derivatives of A coefficients */
  double dRePhiS_dr = 1.0/pow(beta,4)*(
	dA_dr[4]*R[8]
	+ (beta*dA_dr[3]-4*alpha*dA_dr[4])*R[6]
	+ (pow(beta,2)*dA_dr[2]-3*alpha*beta*dA_dr[3]+6*pow(alpha,2)*dA_dr[4])*R[4]
	+ (pow(beta,3)*dA_dr[1]-2*alpha*pow(beta,2)*dA_dr[2]+3*pow(alpha,2)*beta*dA_dr[3]-4*pow(alpha,3)*dA_dr[4])*R[2]
	+ (pow(beta,4)*dA_dr[0]-alpha*pow(beta,3)*dA_dr[1]+pow(alpha,2)*pow(beta,2)*dA_dr[2]-pow(alpha,3)*beta*dA_dr[3]+pow(alpha,4)*dA_dr[4])*R[0]);

  double d2RePhiS_dr2 = 1.0/pow(beta,4)*(
	d2A_dr2[4]*R[8]
	+ (beta*d2A_dr2[3]-4*alpha*d2A_dr2[4])*R[6]
	+ (pow(beta,2)*d2A_dr2[2]-3*alpha*beta*d2A_dr2[3]+6*pow(alpha,2)*d2A_dr2[4])*R[4]
	+ (pow(beta,3)*d2A_dr2[1]-2*alpha*pow(beta,2)*d2A_dr2[2]+3*pow(alpha,2)*beta*d2A_dr2[3]-4*pow(alpha,3)*d2A_dr2[4])*R[2]
	+ (pow(beta,4)*d2A_dr2[0]-alpha*pow(beta,3)*d2A_dr2[1]+pow(alpha,2)*pow(beta,2)*d2A_dr2[2]-pow(alpha,3)*beta*d2A_dr2[3]+pow(alpha,4)*d2A_dr2[4])*R[0]);

  double dRePhiS_dtheta = 1.0/pow(beta,4)*(
	dA_dtheta[4]*R[8]
	+ (beta*dA_dtheta[3]-4*alpha*dA_dtheta[4])*R[6]
	+ (pow(beta,2)*dA_dtheta[2]-3*alpha*beta*dA_dtheta[3]+6*pow(alpha,2)*dA_dtheta[4])*R[4]
	+ (pow(beta,3)*dA_dtheta[1]-2*alpha*pow(beta,2)*dA_dtheta[2]+3*pow(alpha,2)*beta*dA_dtheta[3]-4*pow(alpha,3)*dA_dtheta[4])*R[2]
	+ (pow(beta,4)*dA_dtheta[0]-alpha*pow(beta,3)*dA_dtheta[1]+pow(alpha,2)*pow(beta,2)*dA_dtheta[2]-pow(alpha,3)*beta*dA_dtheta[3]+pow(alpha,4)*dA_dtheta[4])*R[0]);

  double d2RePhiS_dtheta2 = 1.0/pow(beta,4)*(
	d2A_dtheta2[4]*R[8]
	+ (beta*d2A_dtheta2[3]-4*alpha*d2A_dtheta2[4])*R[6]
	+ (pow(beta,2)*d2A_dtheta2[2]-3*alpha*beta*d2A_dtheta2[3]+6*pow(alpha,2)*d2A_dtheta2[4])*R[4]
	+ (pow(beta,3)*d2A_dtheta2[1]-2*alpha*pow(beta,2)*d2A_dtheta2[2]+3*pow(alpha,2)*beta*d2A_dtheta2[3]-4*pow(alpha,3)*d2A_dtheta2[4])*R[2]
	+ (pow(beta,4)*d2A_dtheta2[0]-alpha*pow(beta,3)*d2A_dtheta2[1]+pow(alpha,2)*pow(beta,2)*d2A_dtheta2[2]-pow(alpha,3)*beta*d2A_dtheta2[3]+pow(alpha,4)*d2A_dtheta2[4])*R[0]);

  /* Next account for derivatives of alpha and R */
  double dRePhiS_dalpha = 1.0/pow(beta,4)*(
	+ (-4*A[4])*R[6]
	+ (-3*beta*A[3]+12*alpha*A[4])*R[4]
	+ (-2*pow(beta,2)*A[2]+6*alpha*beta*A[3]-12*pow(alpha,2)*A[4])*R[2]
	+ (-pow(beta,3)*A[1]+2*alpha*pow(beta,2)*A[2]-3*pow(alpha,2)*beta*A[3]+4*pow(alpha,3)*A[4])*R[0]);

  double dRePhiS_drho = 1.0/pow(beta,4)*(
	A[4]*DrhoR[8]
	+ (beta*A[3]-4*alpha*A[4])*DrhoR[6]
	+ (pow(beta,2)*A[2]-3*alpha*beta*A[3]+6*pow(alpha,2)*A[4])*DrhoR[4]
	+ (pow(beta,3)*A[1]-2*alpha*pow(beta,2)*A[2]+3*pow(alpha,2)*beta*A[3]-4*pow(alpha,3)*A[4])*DrhoR[2]
	+ (pow(beta,4)*A[0]-alpha*pow(beta,3)*A[1]+pow(alpha,2)*pow(beta,2)*A[2]-pow(alpha,3)*beta*A[3]+pow(alpha,4)*A[4])*DrhoR[0]);

  double d2RePhiS_dalphadr = 1.0/pow(beta,4)*(
	+ (-4*dA_dr[4])*R[6]
	+ (-3*beta*dA_dr[3]+12*alpha*dA_dr[4])*R[4]
	+ (-2*pow(beta,2)*dA_dr[2]+6*alpha*beta*dA_dr[3]-12*pow(alpha,2)*dA_dr[4])*R[2]
	+ (-pow(beta,3)*dA_dr[1]+2*alpha*pow(beta,2)*dA_dr[2]-3*pow(alpha,2)*beta*dA_dr[3]+4*pow(alpha,3)*dA_dr[4])*R[0]);

  double d2RePhiS_drhodr = 1.0/pow(beta,4)*(
	dA_dr[4]*DrhoR[8]
	+ (beta*dA_dr[3]-4*alpha*dA_dr[4])*DrhoR[6]
	+ (pow(beta,2)*dA_dr[2]-3*alpha*beta*dA_dr[3]+6*pow(alpha,2)*dA_dr[4])*DrhoR[4]
	+ (pow(beta,3)*dA_dr[1]-2*alpha*pow(beta,2)*dA_dr[2]+3*pow(alpha,2)*beta*dA_dr[3]-4*pow(alpha,3)*dA_dr[4])*DrhoR[2]
	+ (pow(beta,4)*dA_dr[0]-alpha*pow(beta,3)*dA_dr[1]+pow(alpha,2)*pow(beta,2)*dA_dr[2]-pow(alpha,3)*beta*dA_dr[3]+pow(alpha,4)*dA_dr[4])*DrhoR[0]);

  double d2RePhiS_dalphadtheta = 1.0/pow(beta,4)*(
	+ (-4*dA_dtheta[4])*R[6]
	+ (-3*beta*dA_dtheta[3]+12*alpha*dA_dtheta[4])*R[4]
	+ (-2*pow(beta,2)*dA_dtheta[2]+6*alpha*beta*dA_dtheta[3]-12*pow(alpha,2)*dA_dtheta[4])*R[2]
	+ (-pow(beta,3)*dA_dtheta[1]+2*alpha*pow(beta,2)*dA_dtheta[2]-3*pow(alpha,2)*beta*dA_dtheta[3]+4*pow(alpha,3)*dA_dtheta[4])*R[0]);

  double d2RePhiS_drhodtheta = 1.0/pow(beta,4)*(
	dA_dtheta[4]*DrhoR[8]
	+ (beta*dA_dtheta[3]-4*alpha*dA_dtheta[4])*DrhoR[6]
	+ (pow(beta,2)*dA_dtheta[2]-3*alpha*beta*dA_dtheta[3]+6*pow(alpha,2)*dA_dtheta[4])*DrhoR[4]
	+ (pow(beta,3)*dA_dtheta[1]-2*alpha*pow(beta,2)*dA_dtheta[2]+3*pow(alpha,2)*beta*dA_dtheta[3]-4*pow(alpha,3)*dA_dtheta[4])*DrhoR[2]
	+ (pow(beta,4)*dA_dtheta[0]-alpha*pow(beta,3)*dA_dtheta[1]+pow(alpha,2)*pow(beta,2)*dA_dtheta[2]-pow(alpha,3)*beta*dA_dtheta[3]+pow(alpha,4)*dA_dtheta[4])*DrhoR[0]);

  double d2RePhiS_dalpha2 = 1.0/pow(beta,4)*(
	+ (12*A[4])*R[4]
	+ (6*beta*A[3]-24*alpha*A[4])*R[2]
	+ (2*pow(beta,2)*A[2]-6*alpha*beta*A[3]+12*pow(alpha,2)*A[4])*R[0]);

  double d2RePhiS_dalphadrho = 1.0/pow(beta,4)*(
  	+ (-4*A[4])*DrhoR[6]
  	+ (-3*beta*A[3]+12*alpha*A[4])*DrhoR[4]
  	+ (-2*pow(beta,2)*A[2]+6*alpha*beta*A[3]-12*pow(alpha,2)*A[4])*DrhoR[2]
  	+ (-pow(beta,3)*A[1]+2*alpha*pow(beta,2)*A[2]-3*pow(alpha,2)*beta*A[3]+4*pow(alpha,3)*A[4])*DrhoR[0]);

  double d2RePhiS_drho2 = 1.0/pow(beta,4)*(
	A[4]*D2rhoR[8]
	+ (beta*A[3]-4*alpha*A[4])*D2rhoR[6]
	+ (pow(beta,2)*A[2]-3*alpha*beta*A[3]+6*pow(alpha,2)*A[4])*D2rhoR[4]
	+ (pow(beta,3)*A[1]-2*alpha*pow(beta,2)*A[2]+3*pow(alpha,2)*beta*A[3]-4*pow(alpha,3)*A[4])*D2rhoR[2]
	+ (pow(beta,4)*A[0]-alpha*pow(beta,3)*A[1]+pow(alpha,2)*pow(beta,2)*A[2]-pow(alpha,3)*beta*A[3]+pow(alpha,4)*A[4])*D2rhoR[0]);

  /* Total first derivatives of PhiS */
  dPhiS_dt  = - m * om * (RePhiS); // This should be interpreted as pure-imaginary
  dPhiS_dr  = dRePhiS_dr + dRePhiS_dalpha*dalpha_dr + dRePhiS_drho*drho_dr;
  dPhiS_dth = dRePhiS_dtheta + dRePhiS_dalpha*dalpha_dtheta + dRePhiS_drho*drho_dtheta;
  dPhiS_dph = m * (RePhiS); // This should be interpreted as pure-imaginary

  /* Total second derivatives of PhiS */
  d2PhiS_dr2   = d2RePhiS_dr2 + d2RePhiS_dalpha2*pow(dalpha_dr,2) + dRePhiS_dalpha*d2alpha_dr2 + d2RePhiS_drho2*pow(drho_dr,2) + dRePhiS_drho*d2rho_dr2
	  + 2*(d2RePhiS_dalphadrho*dalpha_dr*drho_dr + d2RePhiS_dalphadr*dalpha_dr + d2RePhiS_drhodr*drho_dr);
  d2PhiS_dth2  = d2RePhiS_dtheta2 + d2RePhiS_dalpha2*pow(dalpha_dtheta,2) + dRePhiS_dalpha*d2alpha_dtheta2 + d2RePhiS_drho2*pow(drho_dtheta,2) + dRePhiS_drho*d2rho_dtheta2
	  + 2*(d2RePhiS_dalphadrho*dalpha_dtheta*drho_dtheta + d2RePhiS_dalphadtheta*dalpha_dtheta + d2RePhiS_drhodtheta*drho_dtheta);
  d2PhiS_dph2  = - m*m*(RePhiS);
  d2PhiS_dt2   = - m*m*om*om*(RePhiS);
  d2PhiS_dtph = m*m*om*(RePhiS);
  d2PhiS_dtr   = NAN;
  d2PhiS_dtth  = NAN;
  d2PhiS_drth  = NAN;
  d2PhiS_drph  = NAN;
  d2PhiS_dthph = NAN;

  /* Box[PhiS] */
  double sinth  = sin(theta);
  double sinth2 = sinth*sinth;
  double sin2th = sin(2*theta);
  double cos2th = cos(2*theta);
  double r2 = r*r;
  double r3 = r2*r;
  double r4 = r2*r2;
  double a2 = a*a;
  double a4 = a2*a2;

  double effsrc = -((2*a2*dPhiS_dr - a2*d2PhiS_dph2 - a4*d2PhiS_dr2 - a2*d2PhiS_dth2 - 4*dPhiS_dr*r - 2*a2*dPhiS_dr*r +
         4*d2PhiS_dph2*r + 2*a*d2PhiS_dtph*r + 4*a2*d2PhiS_dr2*r + 2*d2PhiS_dth2*r + 6*dPhiS_dr*r2 - 2*d2PhiS_dph2*r2 - 4*d2PhiS_dr2*r2 -
         2*a2*d2PhiS_dr2*r2 - d2PhiS_dth2*r2 - 2*dPhiS_dr*r3 + 4*d2PhiS_dr2*r3 - d2PhiS_dr2*r4 +
         (a4*d2PhiS_dt2 + 4*a*d2PhiS_dtph*r + 2*d2PhiS_dt2*r4 + a2*d2PhiS_dt2*r*(2 + 3*r))*sinth2 +
         cos2th*(a4*d2PhiS_dr2 - 2*a*d2PhiS_dtph*r + (-2 + r)*r*(d2PhiS_dth2 + 2*dPhiS_dr*(-1 + r) + d2PhiS_dr2*(-2 + r)*r) +
            a2*(-d2PhiS_dph2 + d2PhiS_dth2 + 2*dPhiS_dr*(-1 + r) - 4*d2PhiS_dr2*r + 2*d2PhiS_dr2*r2) +
            a2*d2PhiS_dt2*(a2 + (-2 + r)*r)*sinth2) - a2*dPhiS_dth*sin2th + 2*dPhiS_dth*r*sin2th -
         dPhiS_dth*r2*sin2th))/((sinth2*(a2 + (-2 + r)*r)*(a2 + 2*r2 + a2*cos2th)));

  /* Store calculated quantities into the arrays provided by the caller,
     including the phase factor exp(-i*m*phi_p) */
  double cosmph = cos(m*xp.phi);
  double sinmph = sin(m*xp.phi);

  PhiS[0] = RePhiS*cosmph;
  PhiS[1] = - RePhiS*sinmph;

  dPhiS_dx[0] = dPhiS_dt*sinmph;
  dPhiS_dx[1] = dPhiS_dt*cosmph;
  dPhiS_dx[2] = dPhiS_dr*cosmph;
  dPhiS_dx[3] = - dPhiS_dr*sinmph;
  dPhiS_dx[4] = dPhiS_dth*cosmph;
  dPhiS_dx[5] = - dPhiS_dth*sinmph;
  dPhiS_dx[6] = dPhiS_dph*sinmph;
  dPhiS_dx[7] = dPhiS_dph*cosmph;

  d2PhiS_dx2[0]  = d2PhiS_dt2*cosmph;
  d2PhiS_dx2[1]  = - d2PhiS_dt2*sinmph;
  d2PhiS_dx2[2]  = d2PhiS_dtr*cosmph;
  d2PhiS_dx2[3]  = - d2PhiS_dtr*sinmph;
  d2PhiS_dx2[4]  = d2PhiS_dtth*cosmph;
  d2PhiS_dx2[5]  = - d2PhiS_dtth*sinmph;
  d2PhiS_dx2[6]  = d2PhiS_dtph*cosmph;
  d2PhiS_dx2[7]  = - d2PhiS_dtph*sinmph;
  d2PhiS_dx2[8]  = d2PhiS_dr2*cosmph;
  d2PhiS_dx2[9]  = - d2PhiS_dr2*sinmph;
  d2PhiS_dx2[10] = d2PhiS_drth*cosmph;
  d2PhiS_dx2[11] = - d2PhiS_drth*sinmph;
  d2PhiS_dx2[12] = d2PhiS_drph*cosmph;
  d2PhiS_dx2[13] = - d2PhiS_drph*sinmph;
  d2PhiS_dx2[14] = d2PhiS_dth2*cosmph;
  d2PhiS_dx2[15] = - d2PhiS_dth2*sinmph;
  d2PhiS_dx2[16] = d2PhiS_dthph*cosmph;
  d2PhiS_dx2[17] = - d2PhiS_dthph*sinmph;
  d2PhiS_dx2[18] = d2PhiS_dph2*cosmph;
  d2PhiS_dx2[19] = - d2PhiS_dph2*sinmph;

  src[0] = effsrc*cosmph;
  src[1] = - effsrc*sinmph;
}

/* Initialize array of coefficients of pows of dr, dtheta and dphi. */
void effsource_init(double mass, double spin)
{
  M = mass;
  a = spin;
}

/* Initialize array of coefficients of pows of dr, dtheta and dphi. */
void effsource_set_particle(struct coordinate * x_p, double E, double L, double ur_p)
{
  xp = *x_p;
  double r = xp.r;

  /* Compute A coefficients */
  {
	A006 = 64*pow(pow(L,2) + pow(r,2) + (pow(a,2)*(2*M + r))/r,3);
	A008 = (-32*pow(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)),2)*(-(pow(a,8)*pow(M,2)) + 4*pow(a,7)*E*L*pow(M,2) + 4*a*E*L*M*(2*M - 3*r)*pow(r,6) - 4*pow(a,3)*E*L*M*pow(r,3)*(22*pow(M,2) + 3*M*r + 4*pow(r,2)) - 4*pow(a,5)*E*L*M*(16*pow(M,3) + 18*pow(M,2)*r + pow(r,3)) + pow(r,6)*(2*M*pow(r,2)*(-2*M + r + 2*pow(E,2)*r) + pow(L,2)*(-16*pow(M,2) + 14*M*r - 3*pow(r,2))) + 2*pow(a,6)*M*(-2*pow(L,2)*M + 2*pow(E,2)*pow(2*M + r,3) + r*(2*pow(M,2) - M*r + pow(r,2))) + pow(a,2)*pow(r,3)*(2*M*pow(r,2)*(4*pow(M,2) + 6*(-1 + 2*pow(E,2))*M*r + 3*(1 + 2*pow(E,2))*pow(r,2)) + pow(L,2)*(32*pow(M,3) - 24*pow(M,2)*r + 16*M*pow(r,2) - 3*pow(r,3))) + pow(a,4)*M*(4*pow(L,2)*(8*pow(M,3) + 6*pow(M,2)*r - 3*M*pow(r,2) + 2*pow(r,3)) + pow(r,2)*(-4*pow(M,3) + 3*(-3 + 16*pow(E,2))*M*pow(r,2) + 6*(1 + 2*pow(E,2))*pow(r,3) + 4*pow(M,2)*(r + 12*pow(E,2)*r)))))/(3.*pow(r,8)*(pow(a,2) + r*(-2*M + r)));
	A024 = 48*pow(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)),2);
	A026 = (-16*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)))*(-2*pow(a,7)*E*L*M*(24*pow(M,2) + 16*M*r + 3*pow(r,2)) + 2*a*E*L*M*pow(r,5)*(pow(L,2)*(10*M - 9*r) + (14*M - 15*r)*pow(r,2)) - 2*pow(a,5)*E*L*M*r*(16*pow(M,3) + 56*pow(M,2)*r + 50*M*pow(r,2) + 17*pow(r,3) + 3*pow(L,2)*(4*M + r)) + pow(a,8)*M*(6*pow(E,2)*pow(2*M + r,2) - r*(7*M + 3*r)) + pow(r,5)*(2*pow(L,2)*pow(r,2)*(-18*pow(M,2) + (21 + 2*pow(E,2))*M*r - 6*pow(r,2)) + pow(r,4)*(-8*pow(M,2) + 2*(5 + 4*pow(E,2))*M*r - 3*pow(r,2)) - 2*pow(L,4)*(8*pow(M,2) - 10*M*r + 3*pow(r,2))) + 2*pow(a,3)*E*L*M*pow(r,2)*(2*pow(L,2)*(4*pow(M,2) - 7*M*r - 6*pow(r,2)) - pow(r,2)*(16*pow(M,2) + 28*M*r + 29*pow(r,2))) + pow(a,4)*r*(12*pow(L,4)*pow(M,2) + pow(L,2)*(16*pow(M,4) - 8*(-5 + pow(E,2))*pow(M,3)*r + 16*(2 + pow(E,2))*pow(M,2)*pow(r,2) + 2*(-3 + 8*pow(E,2))*M*pow(r,3) - 9*pow(r,4)) + pow(r,2)*(4*pow(M,4) + 16*(2 + 3*pow(E,2))*pow(M,3)*r + 3*(3 + 32*pow(E,2))*pow(M,2)*pow(r,2) + 3*(-5 + 14*pow(E,2))*M*pow(r,3) - 9*pow(r,4))) + pow(a,2)*pow(r,2)*(pow(L,2)*pow(r,2)*(-4*pow(M,3) + 4*(6 + pow(E,2))*pow(M,2)*r + (33 + 14*pow(E,2))*M*pow(r,2) - 21*pow(r,3)) + pow(L,4)*(-8*pow(M,3) + 12*pow(M,2)*r + 8*M*pow(r,2) - 6*pow(r,3)) + pow(r,4)*(4*pow(M,3) + 3*(3 + 10*pow(E,2))*M*pow(r,2) - 9*pow(r,3) + 12*pow(M,2)*(r + 3*pow(E,2)*r))) + pow(a,6)*(pow(L,2)*M*(24*pow(M,2) + 4*(2 + 3*pow(E,2))*M*r + 3*(-1 + 2*pow(E,2))*pow(r,2)) + r*(-(r*(-12*pow(M,3) + 18*pow(M,2)*r + 17*M*pow(r,2) + 3*pow(r,3))) + 2*pow(E,2)*M*(8*pow(M,3) + 36*pow(M,2)*r + 42*M*pow(r,2) + 13*pow(r,3))))))/(3.*pow(r,6)*(pow(a,2) + r*(-2*M + r)));
	A042 = 12*pow(r,3)*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)));
	A044 = (-2*(-4*pow(a,7)*E*L*M*(72*pow(M,2) + 59*M*r + 12*pow(r,2)) + 4*a*E*L*M*pow(r,3)*(pow(L,4)*(8*M - 6*r) + 6*pow(L,2)*(6*M - 5*r)*pow(r,2) + 3*(10*M - 9*r)*pow(r,4)) + pow(a,8)*(-13*pow(M,2)*r + 3*pow(r,3) + 36*pow(E,2)*M*pow(2*M + r,2)) - 4*pow(a,5)*E*L*M*r*(-48*pow(M,3) + 50*pow(M,2)*r + 140*M*pow(r,2) + 49*pow(r,3) + 6*pow(L,2)*(8*M + 3*r)) + pow(r,3)*(8*pow(L,6)*M*(-2*M + r) + pow(L,2)*pow(r,4)*(-120*pow(M,2) + 2*(71 + 12*pow(E,2))*M*r - 41*pow(r,2)) + pow(L,4)*pow(r,2)*(-100*pow(M,2) + 4*(25 + pow(E,2))*M*r - 25*pow(r,2)) + pow(r,6)*(-24*pow(M,2) + 2*(19 + 12*pow(E,2))*M*r - 13*pow(r,2))) - 4*pow(a,3)*E*L*M*pow(r,2)*(6*pow(L,4) + pow(L,2)*(-48*pow(M,2) + 36*M*r + 48*pow(r,2)) + pow(r,2)*(-66*pow(M,2) + 55*M*r + 64*pow(r,2))) + pow(a,2)*pow(r,2)*(12*pow(L,6)*M + pow(L,2)*pow(r,2)*(-192*pow(M,3) + 176*pow(M,2)*r + 12*(9 + 8*pow(E,2))*M*pow(r,2) - 73*pow(r,3)) + 2*pow(r,4)*(-12*pow(M,3) + 2*(23 + 18*pow(E,2))*pow(M,2)*r + 9*(1 + 6*pow(E,2))*M*pow(r,2) - 18*pow(r,3)) - 2*pow(L,4)*(48*pow(M,3) + 8*(-3 + pow(E,2))*pow(M,2)*r - 2*(17 + 4*pow(E,2))*M*pow(r,2) + 11*pow(r,3))) + pow(a,4)*r*(3*pow(L,4)*(32*pow(M,2) + 4*(2 + pow(E,2))*M*r + pow(r,2)) + pow(r,2)*(12*pow(M,4) + 116*pow(M,3)*r + (55 + 288*pow(E,2))*pow(M,2)*pow(r,2) + 6*(-13 + 30*pow(E,2))*M*pow(r,3) - 30*pow(r,4)) - 2*pow(L,2)*(48*pow(M,4) + 4*(-7 + 12*pow(E,2))*pow(M,3)*r - 2*(53 + 24*pow(E,2))*pow(M,2)*pow(r,2) + 4*(2 - 15*pow(E,2))*M*pow(r,3) + 13*pow(r,4))) + 2*pow(a,6)*(pow(L,2)*(72*pow(M,3) + 2*(23 + 24*pow(E,2))*pow(M,2)*r + 6*(1 + 4*pow(E,2))*M*pow(r,2) + 3*pow(r,3)) + r*(-(r*(-10*pow(M,3) + 37*pow(M,2)*r + 29*M*pow(r,2) + 2*pow(r,3))) + 6*pow(E,2)*M*(-8*pow(M,3) + 12*pow(M,2)*r + 30*M*pow(r,2) + 11*pow(r,3))))))/(3.*pow(r,3)*(pow(a,2) + r*(-2*M + r)));
	A060 = pow(r,6);
	A062 = (6*pow(a,5)*E*L*M*(12*M + 5*r) + 4*pow(a,3)*E*L*M*r*(6*pow(L,2) - 20*pow(M,2) + 11*M*r + 18*pow(r,2)) - 3*pow(a,6)*(r*(M + r) + 6*pow(E,2)*M*(2*M + r)) + 2*a*E*L*M*pow(r,2)*(-4*pow(L,2)*(4*M - 3*r) + pow(r,2)*(-26*M + 21*r)) + pow(r,2)*(8*pow(L,4)*M*(2*M - r) + 4*pow(r,4)*(2*pow(M,2) - (3 + 2*pow(E,2))*M*r + pow(r,2)) + pow(L,2)*pow(r,2)*(36*pow(M,2) - 4*(8 + pow(E,2))*M*r + 7*pow(r,2))) + pow(a,2)*r*(-12*pow(L,4)*M + 4*pow(L,2)*(10*pow(M,3) + (-3 + 4*pow(E,2))*pow(M,2)*r - 2*(3 + 2*pow(E,2))*M*pow(r,2) + pow(r,3)) + pow(r,2)*(4*pow(M,3) + 4*(-4 + pow(E,2))*pow(M,2)*r + (1 - 34*pow(E,2))*M*pow(r,2) + 5*pow(r,3))) - pow(a,4)*(3*pow(L,2)*(12*pow(M,2) + 4*(1 + pow(E,2))*M*r + pow(r,2)) + 2*r*(r*(-2*pow(M,2) - 5*M*r + pow(r,2)) + 2*pow(E,2)*M*(-10*pow(M,2) + 8*M*r + 11*pow(r,2)))))/(3.*(pow(a,2) + r*(-2*M + r)));
	A080 = (pow(r,3)*(r*(-3*pow(a,2) + r*(-2*M + r)) + (4*M*(-3*pow(a,4)*pow(E,2) + 6*pow(a,3)*E*L + 2*a*E*L*r*(-4*M + 3*r) + pow(a,2)*(-3*pow(L,2) + 4*pow(E,2)*(M - r)*r) - r*(pow(E,2)*pow(r,3) + pow(L,2)*(-4*M + 2*r))))/(pow(a,2) + r*(-2*M + r))))/24.;
	A106 = (-16*pow(pow(L,2) + pow(r,2) + (pow(a,2)*(2*M + r))/r,2)*(-2*pow(a,2)*M + 2*pow(r,3) + (4*L*(pow(a,3)*E*M - pow(a,2)*L*M + 3*a*E*M*pow(r,2) + L*pow(r,2)*(-2*M + r)))/(pow(a,2) + r*(-2*M + r))))/pow(r,2);
	A108 = (-16*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)))*(2*pow(a,11)*E*L*pow(M,2)*(80*pow(M,2) + 73*M*r + 12*pow(r,2)) + 2*pow(a,9)*E*L*M*r*(32*(-1 + 6*pow(E,2))*pow(M,4) + 12*(23 + 24*pow(E,2))*pow(M,3)*r + (353 + 144*pow(E,2))*pow(M,2)*pow(r,2) + 2*(11 + 12*pow(E,2))*M*pow(r,3) - 12*pow(r,4) + 4*pow(L,2)*M*(16*M + 3*r)) + pow(a,12)*pow(M,2)*(3*r*(3*M + 2*r) - 8*pow(E,2)*(10*pow(M,2) + 11*M*r + 3*pow(r,2))) + 2*a*E*L*M*pow(r,5)*(6*pow(L,4)*pow(r,2)*(20*pow(M,2) + 4*pow(E,2)*M*r - 5*pow(r,2)) + 12*pow(L,6)*(4*pow(M,2) - pow(r,2)) + pow(L,2)*pow(r,4)*(136*pow(M,2) + 2*(-37 + 30*pow(E,2))*M*r + 3*pow(r,2)) - 2*pow(r,6)*(34*pow(M,2) - 41*M*r + 12*pow(r,2))) + 4*pow(a,3)*E*L*M*pow(r,3)*(6*pow(L,6)*(12*pow(M,2) - 2*M*r - pow(r,2)) + pow(L,2)*pow(r,3)*(-340*pow(M,3) + 6*(61 + 20*pow(E,2))*pow(M,2)*r + (-83 + 90*pow(E,2))*M*pow(r,2) - 23*pow(r,3)) + 2*pow(r,5)*(-192*pow(M,3) + 81*pow(M,2)*r + 6*(7 + pow(E,2))*M*pow(r,2) - 19*pow(r,3)) + 6*pow(L,4)*r*(20*pow(M,3) + 4*(10 + pow(E,2))*pow(M,2)*r + 2*(-5 + 2*pow(E,2))*M*pow(r,2) - 5*pow(r,3))) + 2*pow(a,5)*E*L*M*pow(r,2)*(6*pow(L,4)*(2*M + r)*(60*pow(M,2) + 2*(-5 + 2*pow(E,2))*M*r - 5*pow(r,2)) + pow(r,3)*(-420*pow(M,4) - 876*pow(M,3)*r + (907 + 144*pow(E,2))*pow(M,2)*pow(r,2) + 24*(5 + 3*pow(E,2))*M*pow(r,3) - 92*pow(r,4)) + pow(L,2)*r*(-608*pow(M,4) + 8*(13 + 30*pow(E,2))*pow(M,3)*r + 48*(17 + 10*pow(E,2))*pow(M,2)*pow(r,2) + 36*(-7 + 5*pow(E,2))*M*pow(r,3) - 79*pow(r,4))) + pow(a,2)*pow(r,3)*(48*pow(L,8)*pow(M,2)*(-2*M + r) - 24*pow(L,6)*M*pow(r,2)*(2*(5 + 4*pow(E,2))*pow(M,2) + (-5 + 4*pow(E,2))*M*r - 2*pow(E,2)*pow(r,2)) + 2*M*pow(r,7)*(60*pow(M,3) + 4*(-22 + 25*pow(E,2))*pow(M,2)*r + (51 - 86*pow(E,2))*M*pow(r,2) + (-11 + 23*pow(E,2))*pow(r,3)) + pow(L,2)*pow(r,5)*(752*pow(M,4) - 4*(197 + 66*pow(E,2))*pow(M,3)*r + 8*(28 - 11*pow(E,2))*pow(M,2)*pow(r,2) + 3*(1 + 48*pow(E,2))*M*pow(r,3) - 6*pow(r,4)) + 2*pow(L,4)*pow(r,3)*(232*pow(M,4) - 36*(9 + 10*pow(E,2))*pow(M,3)*r - 2*(-71 + 60*pow(E,2))*pow(M,2)*pow(r,2) + 5*(-5 + 18*pow(E,2))*M*pow(r,3) + 3*pow(r,4))) + 2*pow(a,7)*E*L*M*r*(2*pow(L,2)*(288*pow(M,4) + 8*(29 + 15*pow(E,2))*pow(M,3)*r + 2*(71 + 60*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-26 + 15*pow(E,2))*M*pow(r,3) - 15*pow(r,4)) + r*(-256*pow(M,5) - 556*pow(M,4)*r + 8*(-7 + 36*pow(E,2))*pow(M,3)*pow(r,2) + (827 + 288*pow(E,2))*pow(M,2)*pow(r,3) + 4*(11 + 18*pow(E,2))*M*pow(r,4) - 52*pow(r,5))) - pow(a,4)*pow(r,2)*(24*pow(L,6)*M*(20*pow(M,3) + 12*pow(E,2)*pow(M,2)*r + (-5 + 2*pow(E,2))*M*pow(r,2) - pow(E,2)*pow(r,3)) + 2*pow(L,4)*M*r*(-304*pow(M,4) + 8*(11 + 60*pow(E,2))*pow(M,3)*r + 60*(5 + 14*pow(E,2))*pow(M,2)*pow(r,2) + 60*(-3 + pow(E,2))*M*pow(r,3) + (11 - 90*pow(E,2))*pow(r,4)) + M*pow(r,5)*(48*pow(M,4) - 4*(31 + 168*pow(E,2))*pow(M,3)*r - 8*(-34 + 21*pow(E,2))*pow(M,2)*pow(r,2) + (-199 + 448*pow(E,2))*M*pow(r,3) + 12*(4 - 7*pow(E,2))*pow(r,4)) + pow(L,2)*pow(r,3)*(-296*pow(M,5) - 16*(50 + 57*pow(E,2))*pow(M,4)*r + 2*(627 + 452*pow(E,2))*pow(M,3)*pow(r,2) + 2*(-215 + 12*pow(E,2))*pow(M,2)*pow(r,3) - 6*(-3 + 38*pow(E,2))*M*pow(r,4) + 3*pow(r,5))) + pow(a,6)*M*r*(-4*pow(L,4)*(96*pow(M,4) + 4*(11 + 90*pow(E,2))*pow(M,3)*r + 16*(2 + 15*pow(E,2))*pow(M,2)*pow(r,2) - 38*M*pow(r,3) - 15*pow(E,2)*pow(r,4)) + pow(r,3)*(-56*pow(M,5) + 4*(9 + 128*pow(E,2))*pow(M,4)*r + 2*(41 + 520*pow(E,2))*pow(M,3)*pow(r,2) - (157 + 464*pow(E,2))*pow(M,2)*pow(r,3) + (163 - 568*pow(E,2))*M*pow(r,4) + 4*(-13 + 19*pow(E,2))*pow(r,5)) + pow(L,2)*r*(256*pow(M,5) + 8*(47 + 76*pow(E,2))*pow(M,4)*r - 8*(-3 + 76*pow(E,2))*pow(M,3)*pow(r,2) - 2*(415 + 716*pow(E,2))*pow(M,2)*pow(r,3) + 2*(143 - 8*pow(E,2))*M*pow(r,4) + (-17 + 160*pow(E,2))*pow(r,5))) - (2*M - r)*pow(r,8)*(2*pow(E,2)*M*(12*pow(L,6) + 30*pow(L,4)*pow(r,2) + 17*pow(L,2)*pow(r,4) + 5*pow(r,6)) + (2*M - r)*r*(2*pow(L,4)*(8*M - 3*r) + 4*M*pow(r,4) + pow(L,2)*pow(r,2)*(2*M + 3*r))) - pow(a,10)*M*(2*pow(L,2)*M*(40*pow(M,2) + (29 + 32*pow(E,2))*M*r + 3*(-1 + 4*pow(E,2))*pow(r,2)) + r*(r*(50*pow(M,3) + 11*pow(M,2)*r - 19*M*pow(r,2) + 6*pow(r,3)) + pow(E,2)*(-32*pow(M,4) + 336*pow(M,3)*r + 504*pow(M,2)*pow(r,2) + 156*M*pow(r,3) - 6*pow(r,4)))) - pow(a,8)*M*r*(64*pow(L,4)*pow(M,2) + 2*pow(L,2)*(16*(-1 + 36*pow(E,2))*pow(M,4) + 4*(27 + 166*pow(E,2))*pow(M,3)*r + (141 + 364*pow(E,2))*pow(M,2)*pow(r,2) + 3*(-13 + 6*pow(E,2))*M*pow(r,3) + 3*(1 - 7*pow(E,2))*pow(r,4)) + r*(r*(-92*pow(M,4) + 20*pow(M,3)*r + 65*pow(M,2)*pow(r,2) - 63*M*pow(r,3) + 28*pow(r,4)) - 2*pow(E,2)*(128*pow(M,5) + 368*pow(M,4)*r + 120*pow(M,3)*pow(r,2) - 424*pow(M,2)*pow(r,3) - 202*M*pow(r,4) + 17*pow(r,5))))))/(3.*pow(r,10)*pow(pow(a,2) + r*(-2*M + r),2));
	A124 = 8*(pow(L,2) + pow(r,2) + (pow(a,2)*(2*M + r))/r)*(-(pow(a,2)*r) - pow(L,2)*r - 3*pow(r,3) - (4*L*(pow(a,3)*E*M - pow(a,2)*L*M + 3*a*E*M*pow(r,2) + L*pow(r,2)*(-2*M + r)))/(pow(a,2) + r*(-2*M + r)));
	A126 = (-4*(2*pow(a,11)*E*L*pow(M,2)*(-16*(-16 + 9*pow(E,2))*pow(M,2) + (293 - 144*pow(E,2))*M*r + 6*(13 - 6*pow(E,2))*pow(r,2)) + 2*a*E*L*M*pow(r,5)*(2*pow(L,4)*pow(r,2)*(328*pow(M,2) + 2*(-79 + 21*pow(E,2))*M*r - 3*pow(r,2)) + 42*pow(L,6)*(4*pow(M,2) - pow(r,2)) + 3*pow(L,2)*pow(r,4)*(128*pow(M,2) + 2*(-37 + 20*pow(E,2))*M*r + 5*pow(r,2)) - 2*pow(r,6)*(118*pow(M,2) - 125*M*r + 33*pow(r,2))) - (2*M - r)*pow(r,7)*(2*pow(r,6)*(16*pow(M,2) + 2*(-7 + 10*pow(E,2))*M*r + 3*pow(r,2)) + 4*pow(L,6)*(16*pow(M,2) + (-20 + 21*pow(E,2))*M*r + 6*pow(r,2)) + 2*pow(L,2)*pow(r,4)*(10*pow(M,2) + (-17 + 66*pow(E,2))*M*r + 6*pow(r,2)) + pow(L,4)*pow(r,2)*(88*pow(M,2) + 2*(-61 + 94*pow(E,2))*M*r + 39*pow(r,2))) + pow(a,12)*M*(3*r*(21*pow(M,2) + 23*M*r + 6*pow(r,2)) - 2*pow(E,2)*(128*pow(M,3) + 172*pow(M,2)*r + 72*M*pow(r,2) + 9*pow(r,3))) + 4*pow(a,3)*E*L*M*pow(r,3)*(21*pow(L,6)*(12*pow(M,2) - 2*M*r - pow(r,2)) + pow(r,5)*(-1104*pow(M,3) + 4*(82 + 27*pow(E,2))*pow(M,2)*r + (311 - 18*pow(E,2))*M*pow(r,2) - 103*pow(r,3)) + 2*pow(L,4)*r*(256*pow(M,3) + 2*(127 + 21*pow(E,2))*pow(M,2)*r + 6*(-24 + 7*pow(E,2))*M*pow(r,2) - 21*pow(r,3)) - 2*pow(L,2)*pow(r,3)*(102*pow(M,3) - 2*(92 + 87*pow(E,2))*pow(M,2)*r + 5*(17 - 15*pow(E,2))*M*pow(r,2) + 11*pow(r,3))) + 2*pow(a,5)*E*L*M*pow(r,2)*(2*pow(L,4)*(1224*pow(M,3) + 28*(-1 + 3*pow(E,2))*pow(M,2)*r + 2*(-95 + 21*pow(E,2))*M*pow(r,2) - 39*pow(r,3)) + pow(r,3)*(-3300*pow(M,4) + 4*(-517 + 216*pow(E,2))*pow(M,3)*r + (2435 + 504*pow(E,2))*pow(M,2)*pow(r,2) + 12*(59 - 12*pow(E,2))*M*pow(r,3) - 250*pow(r,4)) + pow(L,2)*r*(-192*pow(M,4) + 8*(169 + 114*pow(E,2))*pow(M,3)*r + 72*(4 + 17*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-131 + 102*pow(E,2))*M*pow(r,3) - 131*pow(r,4))) - 2*pow(a,9)*E*L*M*(2*pow(L,2)*M*(216*pow(M,2) + 4*(2 + 9*pow(E,2))*M*r + 3*(-17 + 6*pow(E,2))*pow(r,2)) + r*(-288*(1 + 3*pow(E,2))*pow(M,4) - 4*(275 + 144*pow(E,2))*pow(M,3)*r + 3*(-403 + 72*pow(E,2))*pow(M,2)*pow(r,2) + 6*(-43 + 24*pow(E,2))*M*pow(r,3) + 36*pow(r,4))) + pow(a,2)*pow(r,3)*(-168*pow(L,8)*pow(M,2)*(2*M - r) - 2*pow(L,6)*r*(56*pow(M,4) + 12*(23 + 28*pow(E,2))*pow(M,3)*r + 14*(-13 + 12*pow(E,2))*pow(M,2)*pow(r,2) + 3*(13 - 28*pow(E,2))*M*pow(r,3) - 12*pow(r,4)) + 2*pow(r,7)*(132*pow(M,4) + 4*(-39 + 64*pow(E,2))*pow(M,3)*r + (113 - 246*pow(E,2))*pow(M,2)*pow(r,2) + (-58 + 79*pow(E,2))*M*pow(r,3) + 12*pow(r,4)) + 2*pow(L,4)*pow(r,3)*(592*pow(M,4) - 4*(139 + 329*pow(E,2))*pow(M,3)*r + 2*(91 - 22*pow(E,2))*pow(M,2)*pow(r,2) + 5*(-22 + 53*pow(E,2))*M*pow(r,3) + 42*pow(r,4)) + pow(L,2)*pow(r,5)*(2112*pow(M,4) - 4*(419 + 324*pow(E,2))*pow(M,3)*r - 4*(-73 + 4*pow(E,2))*pow(M,2)*pow(r,2) + (-93 + 464*pow(E,2))*M*pow(r,3) + 51*pow(r,4))) + 2*pow(a,7)*E*L*M*r*(-24*pow(L,4)*M*(9*M - r) + 2*pow(L,2)*(1296*pow(M,4) + 8*(35 + 51*pow(E,2))*pow(M,3)*r + 4*(-52 + 57*pow(E,2))*pow(M,2)*pow(r,2) - (19 + 6*pow(E,2))*M*pow(r,3) - 36*pow(r,4)) + r*(-1600*pow(M,5) + 4*(-515 + 216*pow(E,2))*pow(M,4)*r + 72*(13 + 22*pow(E,2))*pow(M,3)*pow(r,2) + 3*(857 + 72*pow(E,2))*pow(M,2)*pow(r,3) + 12*(43 - 18*pow(E,2))*M*pow(r,4) - 146*pow(r,5))) + pow(a,4)*pow(r,2)*(-4*pow(L,6)*M*(408*pow(M,3) + 4*(-32 + 63*pow(E,2))*pow(M,2)*r + 2*(-26 + 21*pow(E,2))*M*pow(r,2) - (5 + 21*pow(E,2))*pow(r,3)) + 2*pow(r,5)*(192*pow(M,5) + 24*(-1 + 53*pow(E,2))*pow(M,4)*r + 12*(-18 + 11*pow(E,2))*pow(M,3)*pow(r,2) + 4*(27 - 152*pow(E,2))*pow(M,2)*pow(r,3) + (-47 + 107*pow(E,2))*M*pow(r,4) + 18*pow(r,5)) + pow(L,4)*r*(1056*pow(M,5) - 80*(12 + 47*pow(E,2))*pow(M,4)*r - 8*(60 + 491*pow(E,2))*pow(M,3)*pow(r,2) + 12*(12 + 65*pow(E,2))*pow(M,2)*pow(r,3) + 2*(38 + 239*pow(E,2))*M*pow(r,4) + 45*pow(r,5)) + pow(L,2)*pow(r,3)*(3096*pow(M,5) + 12*(55 - 192*pow(E,2))*pow(M,4)*r - 2*(1367 + 788*pow(E,2))*pow(M,3)*pow(r,2) + (275 + 1116*pow(E,2))*pow(M,2)*pow(r,3) + 6*(21 + 94*pow(E,2))*M*pow(r,4) + 66*pow(r,5))) + pow(a,6)*r*(6*pow(L,6)*M*(24*pow(M,2) - 10*M*r + 3*pow(r,2)) - 2*pow(L,4)*M*(864*pow(M,4) + 8*(-1 + 306*pow(E,2))*pow(M,3)*r + 4*(-29 + 164*pow(E,2))*pow(M,2)*pow(r,2) - 32*(-2 + 9*pow(E,2))*M*pow(r,3) - (69 + 59*pow(E,2))*pow(r,4)) + pow(r,3)*(-152*pow(M,6) + 4*(23 + 800*pow(E,2))*pow(M,5)*r + 2*(-31 + 1672*pow(E,2))*pow(M,4)*pow(r,2) - (331 + 1792*pow(E,2))*pow(M,3)*pow(r,3) + (109 - 1704*pow(E,2))*pow(M,2)*pow(r,4) + 4*(8 + 19*pow(E,2))*M*pow(r,5) + 24*pow(r,6)) + pow(L,2)*r*(1600*pow(M,6) - 8*(-167 + 300*pow(E,2))*pow(M,5)*r - 4*(243 + 1228*pow(E,2))*pow(M,4)*pow(r,2) + 6*(-333 + 40*pow(E,2))*pow(M,3)*pow(r,3) + (59 + 1196*pow(E,2))*pow(M,2)*pow(r,4) + (293 + 228*pow(E,2))*M*pow(r,5) + 27*pow(r,6))) - pow(a,10)*M*(pow(L,2)*(-32*(-8 + 27*pow(E,2))*pow(M,3) + (242 - 448*pow(E,2))*pow(M,2)*r + 3*(-19 + 20*pow(E,2))*M*pow(r,2) + 18*(-3 + 2*pow(E,2))*pow(r,3)) + r*(r*(290*pow(M,3) + 165*pow(M,2)*r - 141*M*pow(r,2) - 68*pow(r,3)) + 2*pow(E,2)*(144*pow(M,4) + 680*pow(M,3)*r + 864*pow(M,2)*pow(r,2) + 358*M*pow(r,3) + 37*pow(r,4)))) + pow(a,8)*(2*pow(L,4)*M*(144*pow(M,3) + 8*(-8 + 27*pow(E,2))*pow(M,2)*r + 6*(-6 + 7*pow(E,2))*M*pow(r,2) - 9*(-3 + pow(E,2))*pow(r,3)) - pow(L,2)*M*r*(288*(1 + 18*pow(E,2))*pow(M,4) + 8*(105 + 286*pow(E,2))*pow(M,3)*r - 2*(-589 + 516*pow(E,2))*pow(M,2)*pow(r,2) + (15 - 268*pow(E,2))*M*pow(r,3) + 2*(-93 + 20*pow(E,2))*pow(r,4)) + pow(r,2)*(4*pow(E,2)*M*(400*pow(M,5) + 696*pow(M,4)*r - 12*pow(M,3)*pow(r,2) - 732*pow(M,2)*pow(r,3) - 368*M*pow(r,4) - 19*pow(r,5)) + r*(404*pow(M,5) + 8*pow(M,4)*r - 375*pow(M,3)*pow(r,2) + 103*pow(M,2)*pow(r,3) + 100*M*pow(r,4) + 6*pow(r,5))))))/(3.*pow(r,7)*pow(pow(a,2) + r*(-2*M + r),2));
	A142 = pow(r,2)*(2*pow(a,2)*M - 2*pow(r,3) - (4*L*(pow(a,3)*E*M - pow(a,2)*L*M + 3*a*E*M*pow(r,2) + L*pow(r,2)*(-2*M + r)))/(pow(a,2) + r*(-2*M + r)) - 4*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2))));
	A144 = (2*pow(a,9)*E*L*M*(16*(-16 + 9*pow(E,2))*pow(M,2) + 12*(-13 + 6*pow(E,2))*M*r - 3*pow(r,2)) + (2*M - r)*pow(r,4)*(16*pow(L,6)*M*(2*M - r) + 2*pow(L,2)*pow(r,4)*(-6*pow(M,2) + (-5 + 66*pow(E,2))*M*r + 4*pow(r,2)) + 2*pow(L,4)*pow(r,2)*(34*pow(M,2) + (-33 + 47*pow(E,2))*M*r + 8*pow(r,2)) + pow(r,6)*(48*pow(M,2) + 10*(-5 + 6*pow(E,2))*M*r + 13*pow(r,2))) - 2*a*E*L*M*pow(r,4)*(2*pow(L,4)*(200*pow(M,2) - 94*M*r - 3*pow(r,2)) + 3*pow(L,2)*pow(r,2)*(188*pow(M,2) + 4*(-27 + 5*pow(E,2))*M*r + 7*pow(r,2)) - 3*pow(r,4)*(100*pow(M,2) - 116*M*r + 33*pow(r,2))) + pow(a,10)*M*(-9*r*(13*M + 7*r) + pow(E,2)*(256*pow(M,2) + 264*M*r + 66*pow(r,2))) - 2*pow(a,3)*E*L*M*pow(r,2)*(4*pow(L,4)*(180*pow(M,2) - 92*M*r - 3*pow(r,2)) + pow(L,2)*r*(768*pow(M,3) + 4*(67 + 84*pow(E,2))*pow(M,2)*r + 4*(-149 + 9*pow(E,2))*M*pow(r,2) + 15*pow(r,3)) - 4*pow(r,3)*(426*pow(M,3) - (85 + 108*pow(E,2))*pow(M,2)*r + 12*(-13 + 3*pow(E,2))*M*pow(r,2) + 52*pow(r,3))) - 2*pow(a,5)*E*L*M*r*(-6*pow(L,4)*(18*M + r) + pow(L,2)*(1728*pow(M,3) + 24*(-33 + 10*pow(E,2))*pow(M,2)*r - 4*(83 + 15*pow(E,2))*M*pow(r,2) - 9*pow(r,3)) + 4*r*(-336*pow(M,4) + 2*(-107 + 108*pow(E,2))*pow(M,3)*r + 3*(79 + 36*pow(E,2))*pow(M,2)*pow(r,2) + (161 - 90*pow(E,2))*M*pow(r,3) - 34*pow(r,4))) - pow(a,2)*pow(r,2)*(-2*pow(L,6)*M*(240*pow(M,2) - 154*M*r + 17*pow(r,2)) + 2*pow(L,4)*r*(48*pow(M,4) - 4*(35 + 176*pow(E,2))*pow(M,3)*r + 8*(8 + 15*pow(E,2))*pow(M,2)*pow(r,2) + (-35 + 73*pow(E,2))*M*pow(r,3) + 16*pow(r,4)) + pow(r,5)*(264*pow(M,4) + 4*(-11 + 144*pow(E,2))*pow(M,3)*r + 2*(11 - 222*pow(E,2))*pow(M,2)*pow(r,2) + (-125 + 138*pow(E,2))*M*pow(r,3) + 46*pow(r,4)) + pow(L,2)*pow(r,3)*(1656*pow(M,4) - 12*(65 + 166*pow(E,2))*pow(M,3)*r + 2*(-161 + 374*pow(E,2))*pow(M,2)*pow(r,2) + (47 + 256*pow(E,2))*M*pow(r,3) + 51*pow(r,4))) - pow(a,4)*r*(6*pow(L,6)*M*(12*M + r) + 2*pow(L,4)*(-576*pow(M,4) - 36*(-9 + 20*pow(E,2))*pow(M,3)*r + 2*(-19 + 137*pow(E,2))*pow(M,2)*pow(r,2) + (47 + 29*pow(E,2))*M*pow(r,3) + 8*pow(r,4)) + pow(r,3)*(192*pow(M,5) + 4*(85 + 468*pow(E,2))*pow(M,4)*r + 228*(-1 + 2*pow(E,2))*pow(M,3)*pow(r,2) - 3*(149 + 268*pow(E,2))*pow(M,2)*pow(r,3) - 12*(-9 + pow(E,2))*M*pow(r,4) + 60*pow(r,5)) + pow(L,2)*r*(1344*pow(M,5) + 32*(13 - 105*pow(E,2))*pow(M,4)*r - 4*(243 + 280*pow(E,2))*pow(M,3)*pow(r,2) + 28*(-29 + 75*pow(E,2))*pow(M,2)*pow(r,3) + (391 + 56*pow(E,2))*M*pow(r,4) + 66*pow(r,5))) + 2*pow(a,7)*E*L*M*(4*r*(-8*(5 + 18*pow(E,2))*pow(M,3) + 3*(-59 + 12*pow(E,2))*pow(M,2)*r + (-121 + 72*pow(E,2))*M*pow(r,2) + 6*pow(r,3)) + 3*pow(L,2)*(144*pow(M,2) + pow(r,2) + 4*M*(r + 3*pow(E,2)*r))) + pow(a,8)*(pow(L,2)*M*(-32*(-8 + 27*pow(E,2))*pow(M,2) + 12*(4 - 19*pow(E,2))*M*r + 3*(-41 + 20*pow(E,2))*pow(r,2)) + r*(-(r*(-420*pow(M,3) + 115*pow(M,2)*r + 256*M*pow(r,2) + 7*pow(r,3))) + 4*pow(E,2)*M*(40*pow(M,3) + 178*pow(M,2)*r + 209*M*pow(r,2) + 60*pow(r,3)))) - pow(a,6)*(6*pow(L,4)*M*(48*pow(M,2) + 2*(-5 + 18*pow(E,2))*M*r + (11 + pow(E,2))*pow(r,2)) + pow(L,2)*r*(-32*(5 + 108*pow(E,2))*pow(M,4) + 8*(-88 + 153*pow(E,2))*pow(M,3)*r + 2*(-247 + 658*pow(E,2))*pow(M,2)*pow(r,2) + (417 - 128*pow(E,2))*M*pow(r,3) + 23*pow(r,4)) + pow(r,2)*(4*pow(E,2)*M*(336*pow(M,4) + 324*pow(M,3)*r - 128*pow(M,2)*pow(r,2) - 263*M*pow(r,3) - 66*pow(r,4)) + r*(276*pow(M,4) - 652*pow(M,3)*r - 323*pow(M,2)*pow(r,2) + 350*M*pow(r,3) + 34*pow(r,4)))))/(3.*pow(r,4)*pow(pow(a,2) + r*(-2*M + r),2));
	A160 = -pow(r,5)/2.;
	A162 = (6*pow(a,7)*E*L*M*(12*(-4 + pow(E,2))*M - 11*r) + 3*pow(a,8)*(3*r*(-3*M + r) + 2*pow(E,2)*M*(24*M + 13*r)) + 2*a*E*L*M*(2*M - r)*pow(r,3)*((82*M - 51*r)*pow(r,2) + pow(L,2)*(-124*M + 30*r)) + (2*M - r)*pow(r,3)*(8*pow(L,4)*M*(2*M - r) + pow(L,2)*pow(r,2)*(-36*pow(M,2) + 4*(8 + 11*pow(E,2))*M*r - 7*pow(r,2)) + 8*pow(r,4)*(4*pow(M,2) + (-4 + 5*pow(E,2))*M*r + pow(r,2))) - 2*pow(a,5)*E*L*M*(-12*pow(L,2)*(9*M - 2*r) + r*(8*(-11 + 15*pow(E,2))*pow(M,2) - 24*(-9 + 5*pow(E,2))*M*r + 63*pow(r,2))) - 2*pow(a,3)*E*L*M*r*(pow(L,2)*(360*pow(M,2) - 352*M*r + 54*pow(r,2)) - r*(400*pow(M,3) + 4*(19 - 54*pow(E,2))*pow(M,2)*r + 4*(-40 + 21*pow(E,2))*M*pow(r,2) + 21*pow(r,3))) - pow(a,2)*r*(-8*pow(L,4)*M*(30*pow(M,2) - 29*M*r + 7*pow(r,2)) + pow(r,3)*(88*pow(M,4) + 4*(27 + 64*pow(E,2))*pow(M,3)*r - 2*(57 + 50*pow(E,2))*pow(M,2)*pow(r,2) + (-11 + 26*pow(E,2))*M*pow(r,3) + 15*pow(r,4)) + pow(L,2)*r*(400*pow(M,4) + 8*(-26 + 63*pow(E,2))*pow(M,2)*pow(r,2) + 4*(17 + 4*pow(E,2))*M*pow(r,3) + pow(r,4) - 64*pow(M,3)*(r + 14*pow(E,2)*r))) + pow(a,6)*(-3*pow(L,2)*(24*(-2 + 3*pow(E,2))*pow(M,2) + (4 - 8*pow(E,2))*M*r - 3*pow(r,2)) + r*(r*(86*pow(M,2) - 139*M*r + 19*pow(r,2)) + 2*pow(E,2)*M*(-44*pow(M,2) + 98*M*r + 105*pow(r,2)))) + pow(a,4)*(24*pow(L,4)*M*(-3*M + r) + pow(L,2)*r*(8*(-11 + 90*pow(E,2))*pow(M,3) + 4*(59 - 178*pow(E,2))*pow(M,2)*r + 2*(-41 + 26*pow(E,2))*M*pow(r,2) + pow(r,3)) + pow(r,2)*(r*(-20*pow(M,3) + 296*pow(M,2)*r - 149*M*pow(r,2) + 3*pow(r,3)) + 2*pow(E,2)*M*(-200*pow(M,3) - 108*pow(M,2)*r + 36*M*pow(r,2) + 73*pow(r,3)))))/(12.*r*pow(pow(a,2) + r*(-2*M + r),2));
	A180 = -(pow(r,2)*(60*pow(a,3)*E*L*M + 4*a*E*L*M*r*(8*M + 3*r) - 3*pow(a,4)*(10*pow(E,2)*M + 3*r) + r*(8*pow(L,2)*M*(-2*M + r) + pow(r,2)*(8*pow(M,2) + 2*(-3 + 5*pow(E,2))*M*r + pow(r,2))) - 2*pow(a,2)*(15*pow(L,2)*M + r*(r*(-7*M + 4*r) + 2*pow(E,2)*M*(4*M + 5*r)))))/(48.*(pow(a,2) + r*(-2*M + r)));
	A204 = (48*pow(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)),2))/(pow(a,2) + r*(-2*M + r));
	A206 = (8*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)))*(-4*pow(a,7)*E*L*M*(12*pow(M,2) + 23*M*r + 3*pow(r,2)) - 4*a*E*L*M*pow(r,5)*(pow(L,2)*(67*M - 36*r) + (17*M - 15*r)*pow(r,2)) - 4*pow(a,5)*E*L*M*r*(-64*pow(M,3) - 104*pow(M,2)*r + 55*M*pow(r,2) + 7*pow(r,3) + 3*pow(L,2)*(8*M + r)) + pow(a,8)*M*(-(r*(M + 6*r)) + 12*pow(E,2)*(2*pow(M,2) + 3*M*r + pow(r,2))) + 4*pow(a,3)*E*L*M*pow(r,2)*(pow(L,2)*(16*pow(M,2) - 73*M*r + 3*pow(r,2)) + pow(r,2)*(196*pow(M,2) - 17*M*r + 11*pow(r,2))) + pow(r,5)*(8*pow(L,4)*(10*pow(M,2) - 11*M*r + 3*pow(r,2)) + pow(r,4)*(40*pow(M,2) - 2*(19 + 2*pow(E,2))*M*r + 9*pow(r,2)) + 2*pow(L,2)*pow(r,2)*(54*pow(M,2) + (-57 + 2*pow(E,2))*M*r + 15*pow(r,2))) - 2*pow(a,2)*pow(r,2)*(pow(r,4)*(58*pow(M,3) + 3*(-23 + 26*pow(E,2))*pow(M,2)*r + 33*M*pow(r,2) - 6*pow(r,3)) + 2*pow(L,4)*(8*pow(M,3) - 33*pow(M,2)*r + 10*M*pow(r,2) + 3*pow(r,3)) - pow(L,2)*pow(r,2)*(-194*pow(M,3) + (171 + 104*pow(E,2))*pow(M,2)*r + (-57 + 10*pow(E,2))*M*pow(r,2) + 6*pow(r,3))) + pow(a,4)*r*(48*pow(L,4)*pow(M,2) + pow(r,2)*(4*pow(M,4) + 20*(1 - 18*pow(E,2))*pow(M,3)*r - 3*(-53 + 76*pow(E,2))*pow(M,2)*pow(r,2) + 12*(-7 + 2*pow(E,2))*M*pow(r,3) - 3*pow(r,4)) - 2*pow(L,2)*(64*pow(M,4) + 8*(11 + 2*pow(E,2))*pow(M,3)*r - (109 + 80*pow(E,2))*pow(M,2)*pow(r,2) + (51 - 14*pow(E,2))*M*pow(r,3) + 9*pow(r,4))) + 2*pow(a,6)*(pow(L,2)*M*(12*pow(M,2) + 4*(7 + 6*pow(E,2))*M*r + 3*(-1 + 2*pow(E,2))*pow(r,2)) - r*(pow(r,3)*(31*M + 3*r) + 2*pow(E,2)*M*(32*pow(M,3) + 60*pow(M,2)*r + 9*M*pow(r,2) - 8*pow(r,3))))))/(3.*pow(r,6)*pow(pow(a,2) + r*(-2*M + r),2));
	A222 = (24*pow(r,3)*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2))))/(pow(a,2) + r*(-2*M + r));
	A224 = (2*(4*pow(a,7)*E*L*M*(24*pow(M,2) + 19*M*r + 9*pow(r,2)) - 4*a*E*L*M*pow(r,3)*(4*pow(L,4)*M + 9*pow(L,2)*(12*M - 7*r)*pow(r,2) + 6*(9*M - 7*r)*pow(r,4)) + 4*pow(a,5)*E*L*M*r*(48*pow(M,3) + 142*pow(M,2)*r + 76*M*pow(r,2) + 50*pow(r,3) + pow(L,2)*(-6*M + 9*r)) - pow(a,8)*(24*pow(E,2)*pow(M,2)*(2*M + r) + r*(49*pow(M,2) + 42*M*r + 6*pow(r,2))) + pow(r,3)*(4*pow(L,6)*M*(2*M - r) + 12*pow(L,2)*pow(r,4)*(28*pow(M,2) - 32*M*r + 9*pow(r,2)) + 3*pow(r,6)*(40*pow(M,2) - 2*(23 + 2*pow(E,2))*M*r + 13*pow(r,2)) + 2*pow(L,4)*pow(r,2)*(118*pow(M,2) + (-131 + 2*pow(E,2))*M*r + 36*pow(r,2))) + 4*pow(a,3)*E*L*M*pow(r,3)*(-6*pow(L,2)*(9*M - 8*r) + r*(114*pow(M,2) + 23*M*r + 83*pow(r,2))) + 2*pow(a,2)*pow(r,3)*(pow(L,4)*(4*(18 + pow(E,2))*pow(M,2) + (-71 + 2*pow(E,2))*M*r + 15*pow(r,2)) - 3*pow(r,3)*(12*pow(M,3) + 2*(-7 + 24*pow(E,2))*pow(M,2)*r + (31 + 6*pow(E,2))*M*pow(r,2) - 15*pow(r,3)) + 3*pow(L,2)*r*(-32*pow(M,3) + 2*(19 + 9*pow(E,2))*pow(M,2)*r - 61*M*pow(r,2) + 25*pow(r,3))) + pow(a,4)*r*(6*pow(L,4)*(2*pow(M,2) - 6*M*r - pow(r,2)) - 3*pow(r,2)*(76*pow(M,4) + 12*(-1 + 16*pow(E,2))*pow(M,3)*r + 3*(-19 + 56*pow(E,2))*pow(M,2)*pow(r,2) + 4*(7 + 3*pow(E,2))*M*pow(r,3) - 19*pow(r,4)) - 2*pow(L,2)*(48*pow(M,4) + 68*pow(M,3)*r - 4*(13 + 9*pow(E,2))*pow(M,2)*pow(r,2) + 84*M*pow(r,3) - 15*pow(r,4))) - 2*pow(a,6)*(pow(L,2)*(24*pow(M,3) + (26 - 6*pow(E,2))*pow(M,2)*r + 39*M*pow(r,2) + 6*pow(r,3)) + M*r*(r*(-106*pow(M,2) - 37*M*r + 39*pow(r,2)) + 6*pow(E,2)*(8*pow(M,3) + 36*pow(M,2)*r + 20*M*pow(r,2) + pow(r,3))))))/(3.*pow(r,3)*pow(pow(a,2) + r*(-2*M + r),2));
	A240 = (3*pow(r,6))/(pow(a,2) + r*(-2*M + r));
	A242 = (24*pow(a,5)*E*L*M*(7*M + 4*r) + 12*pow(a,3)*E*L*M*r*(4*pow(L,2) - 16*pow(M,2) + 9*M*r + 23*pow(r,2)) - 12*a*E*L*M*pow(r,2)*(pow(L,2)*(8*M - 4*r) + (19*M - 13*r)*pow(r,2)) - 3*pow(a,6)*(4*pow(E,2)*M*(7*M + 3*r) + r*(10*M + 3*r)) + pow(r,2)*(24*pow(L,4)*M*(2*M - r) + pow(r,4)*(120*pow(M,2) - 2*(71 + 6*pow(E,2))*M*r + 41*pow(r,2)) + pow(L,2)*pow(r,2)*(228*pow(M,2) - 232*M*r + 59*pow(r,2))) + pow(a,2)*r*(-24*pow(L,4)*M + pow(r,2)*(-36*pow(M,3) + (22 - 84*pow(E,2))*pow(M,2)*r - 6*(19 + 10*pow(E,2))*M*pow(r,2) + 67*pow(r,3)) + 2*pow(L,2)*(48*pow(M,3) - 6*(13 + 2*pow(E,2))*M*pow(r,2) + 25*pow(r,3) + 6*pow(M,2)*(r + 4*pow(E,2)*r))) - pow(a,4)*(3*pow(L,2)*(28*pow(M,2) + 4*(5 + 2*pow(E,2))*M*r + 3*pow(r,2)) + r*(r*(-78*pow(M,2) + 14*M*r - 17*pow(r,2)) + 12*pow(E,2)*M*(-8*pow(M,2) + 10*M*r + 7*pow(r,2)))))/(6.*pow(pow(a,2) + r*(-2*M + r),2));
	A260 = (pow(r,3)*(48*pow(a,3)*E*L*M - 3*pow(a,4)*(8*pow(E,2)*M + r) + 16*a*E*L*M*r*(-5*M + 3*r) + r*(20*pow(L,2)*M*(2*M - r) + pow(r,2)*(40*pow(M,2) - 2*(21 + 2*pow(E,2))*M*r + 11*pow(r,2))) - 2*pow(a,2)*(12*pow(L,2)*M + r*((7*M - 4*r)*r + 2*pow(E,2)*M*(-10*M + 7*r)))))/(24.*pow(pow(a,2) + r*(-2*M + r),2));
	A304 = (-8*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2)))*(4*pow(a,3)*E*L*M + pow(a,4)*r + 12*a*E*L*M*pow(r,2) + pow(r,2)*(pow(r,2)*(-5*M + 2*r) + pow(L,2)*(-9*M + 4*r)) + pow(a,2)*(pow(L,2)*(-4*M + r) + r*(2*pow(M,2) - 3*M*r + 3*pow(r,2)))))/(r*pow(pow(a,2) + r*(-2*M + r),2));
	A306 = (4*(-4*pow(a,11)*E*L*pow(M,2)*(4*(35 + 9*pow(E,2))*pow(M,2) + (127 + 54*pow(E,2))*M*r + 3*(5 + 6*pow(E,2))*pow(r,2)) - 2*a*E*L*M*pow(r,5)*(pow(r,6)*(-184*pow(M,2) + 12*(17 + 9*pow(E,2))*M*r - 57*pow(r,2)) + pow(L,2)*pow(r,4)*(212*pow(M,2) + 4*(-23 + 57*pow(E,2))*M*r - 3*pow(r,2)) + 42*pow(L,6)*(4*pow(M,2) - pow(r,2)) + 2*pow(L,4)*pow(r,2)*(486*pow(M,2) + 2*(-170 + 21*pow(E,2))*M*r + 51*pow(r,2))) - 2*pow(a,3)*E*L*M*pow(r,3)*(4*pow(L,4)*pow(r,2)*((829 + 42*pow(E,2))*pow(M,2) + 2*(-107 + 21*pow(E,2))*M*r - 84*pow(r,2)) + 42*pow(L,6)*(12*pow(M,2) - 2*M*r - pow(r,2)) + pow(r,5)*(-3444*pow(M,3) + 2*(961 - 270*pow(E,2))*pow(M,2)*r + 4*(106 + 135*pow(E,2))*M*pow(r,2) - 257*pow(r,3)) + 2*pow(L,2)*pow(r,3)*(-3230*pow(M,3) + 3*(895 + 224*pow(E,2))*pow(M,2)*r + (-109 + 384*pow(E,2))*M*pow(r,2) - 214*pow(r,3))) + pow(r,7)*(2*pow(L,6)*(2*M - r)*(54*pow(M,2) + (-71 + 42*pow(E,2))*M*r + 24*pow(r,2)) + pow(L,4)*pow(r,2)*(272*pow(M,3) + 2*(-247 + 222*pow(E,2))*pow(M,2)*r + (311 - 226*pow(E,2))*M*pow(r,2) - 66*pow(r,3)) + pow(L,2)*pow(r,4)*(220*pow(M,3) + 340*(-1 + pow(E,2))*pow(M,2)*r + (193 - 184*pow(E,2))*M*pow(r,2) - 39*pow(r,3)) + pow(r,6)*(116*pow(M,3) + 2*(-83 + 44*pow(E,2))*pow(M,2)*r + 6*(14 - 9*pow(E,2))*M*pow(r,2) - 15*pow(r,3))) + pow(a,12)*M*(3*r*(-8*pow(M,2) - 3*M*r + 2*pow(r,2)) + pow(E,2)*(280*pow(M,3) + 308*pow(M,2)*r + 72*M*pow(r,2) - 6*pow(r,3))) - 4*pow(a,9)*E*L*M*(pow(L,2)*M*(108*pow(M,2) + 2*(137 + 24*pow(E,2))*M*r + 3*(13 + 6*pow(E,2))*pow(r,2)) + r*(24*(-1 + 12*pow(E,2))*pow(M,4) + 2*(215 + 234*pow(E,2))*pow(M,3)*r + 12*(49 + 39*pow(E,2))*pow(M,2)*pow(r,2) + 9*(-5 + 18*pow(E,2))*M*pow(r,3) - 30*pow(r,4))) - pow(a,2)*pow(r,3)*(-168*pow(L,8)*pow(M,2)*(2*M - r) + pow(L,4)*pow(r,3)*(4168*pow(M,4) - 4*(1441 + 1086*pow(E,2))*pow(M,3)*r + 4*(583 + 6*pow(E,2))*pow(M,2)*pow(r,2) + (-85 + 732*pow(E,2))*M*pow(r,3) - 81*pow(r,4)) + 2*pow(L,6)*r*(240*pow(M,4) - 4*(233 + 84*pow(E,2))*pow(M,3)*r - 6*(-83 + 28*pow(E,2))*pow(M,2)*pow(r,2) + (23 + 84*pow(E,2))*M*pow(r,3) - 36*pow(r,4)) + pow(L,2)*pow(r,5)*(3640*pow(M,4) + 12*(-413 + 58*pow(E,2))*pow(M,3)*r + 2*(1091 - 1004*pow(E,2))*pow(M,2)*pow(r,2) + (-243 + 800*pow(E,2))*M*pow(r,3) - 27*pow(r,4)) + pow(r,7)*(696*pow(M,4) + 4*(-315 + 274*pow(E,2))*pow(M,3)*r + 4*(204 - 241*pow(E,2))*pow(M,2)*pow(r,2) + (-205 + 268*pow(E,2))*M*pow(r,3) + 15*pow(r,4))) - 2*pow(a,5)*E*L*M*pow(r,2)*(2*pow(L,4)*(1008*pow(M,3) + 4*(257 + 21*pow(E,2))*pow(M,2)*r + 2*(-104 + 21*pow(E,2))*M*pow(r,2) - 87*pow(r,3)) + pow(L,2)*r*(-1888*pow(M,4) + 20*(-161 + 24*pow(E,2))*pow(M,3)*r + 288*(21 + 8*pow(E,2))*pow(M,2)*pow(r,2) + 12*(-89 + 74*pow(E,2))*M*pow(r,3) - 617*pow(r,4)) - pow(r,3)*(768*pow(M,4) + 4*(823 + 378*pow(E,2))*pow(M,3)*r + 4*(-1091 + 18*pow(E,2))*pow(M,2)*pow(r,2) - 18*(-1 + 56*pow(E,2))*M*pow(r,3) + 403*pow(r,4))) + pow(a,4)*pow(r,2)*(2*pow(L,6)*M*(672*pow(M,3) + 4*(101 + 126*pow(E,2))*pow(M,2)*r + 28*(-11 + 3*pow(E,2))*M*pow(r,2) - (31 + 42*pow(E,2))*pow(r,3)) + pow(r,5)*(180*pow(M,5) + 12*(1 - 238*pow(E,2))*pow(M,4)*r - 9*(-161 + 60*pow(E,2))*pow(M,3)*pow(r,2) + 8*(-219 + 287*pow(E,2))*pow(M,2)*pow(r,3) + (421 - 538*pow(E,2))*M*pow(r,4) + 45*pow(r,5)) + pow(L,4)*r*(-1888*pow(M,5) + 8*(-239 + 180*pow(E,2))*pow(M,4)*r + 32*(163 + 293*pow(E,2))*pow(M,3)*pow(r,2) + 12*(-223 + 58*pow(E,2))*pow(M,2)*pow(r,3) + 2*(79 - 396*pow(E,2))*M*pow(r,4) + 147*pow(r,5)) + pow(L,2)*pow(r,3)*(-364*pow(M,5) - 4*(763 + 2976*pow(E,2))*pow(M,4)*r + (6129 + 4072*pow(E,2))*pow(M,3)*pow(r,2) + 4*(-965 + 831*pow(E,2))*pow(M,2)*pow(r,3) + 3*(181 - 436*pow(E,2))*M*pow(r,4) + 171*pow(r,5))) - 2*pow(a,7)*E*L*M*r*(48*pow(L,4)*M*(6*M + r) + 2*pow(L,2)*(864*pow(M,4) + 4*(163 + 84*pow(E,2))*pow(M,3)*r + (1271 + 624*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-101 + 96*pow(E,2))*M*pow(r,3) - 96*pow(r,4)) - r*(896*pow(M,5) + 1840*pow(M,4)*r + 36*(17 + 2*pow(E,2))*pow(M,3)*pow(r,2) - 4*(769 + 324*pow(E,2))*pow(M,2)*pow(r,3) + 2*(179 - 432*pow(E,2))*M*pow(r,4) + 263*pow(r,5))) + pow(a,6)*r*(6*pow(L,6)*M*(32*pow(M,2) + 2*M*r + pow(r,2)) + 4*pow(L,4)*M*(288*pow(M,4) + 4*(23 + 252*pow(E,2))*pow(M,3)*r + (389 + 1450*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-125 + 69*pow(E,2))*M*pow(r,3) - (6 + 73*pow(E,2))*pow(r,4)) + pow(r,3)*(184*pow(M,6) - 4*(11 + 268*pow(E,2))*pow(M,5)*r - 6*(37 + 632*pow(E,2))*pow(M,4)*pow(r,2) + (77 + 2080*pow(E,2))*pow(M,3)*pow(r,3) + (-1231 + 2360*pow(E,2))*pow(M,2)*pow(r,4) + 3*(203 - 184*pow(E,2))*M*pow(r,5) + 75*pow(r,6)) + pow(L,2)*r*(-896*pow(M,6) - 32*(38 + 59*pow(E,2))*pow(M,5)*r - 16*(21 + 292*pow(E,2))*pow(M,4)*pow(r,2) + 2*(1643 + 4488*pow(E,2))*pow(M,3)*pow(r,3) + (-1991 + 2372*pow(E,2))*pow(M,2)*pow(r,4) + (571 - 964*pow(E,2))*M*pow(r,5) + 105*pow(r,6))) + pow(a,10)*M*(pow(L,2)*(8*(35 + 54*pow(E,2))*pow(M,3) + 8*(25 + 109*pow(E,2))*pow(M,2)*r + 3*(-7 + 76*pow(E,2))*M*pow(r,2) + 6*(3 - 2*pow(E,2))*pow(r,3)) + r*(r*(134*pow(M,3) - 21*pow(M,2)*r - 69*M*pow(r,2) + 58*pow(r,3)) + pow(E,2)*(-48*pow(M,4) + 1072*pow(M,3)*r + 1704*pow(M,2)*pow(r,2) + 452*M*pow(r,3) - 76*pow(r,4)))) + pow(a,8)*(2*pow(L,4)*M*(72*pow(M,3) + 4*(55 + 72*pow(E,2))*pow(M,2)*r + 78*pow(E,2)*M*pow(r,2) - 3*(-3 + pow(E,2))*pow(r,3)) + pow(L,2)*M*r*(48*(-1 + 72*pow(E,2))*pow(M,4) + 8*(81 + 514*pow(E,2))*pow(M,3)*r + (869 + 5400*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-285 + 472*pow(E,2))*M*pow(r,3) + 4*(24 - 71*pow(E,2))*pow(r,4)) + pow(r,2)*(r*(-264*pow(M,5) + 92*pow(M,4)*r + 135*pow(M,3)*pow(r,2) - 189*pow(M,2)*pow(r,3) + 361*M*pow(r,4) + 30*pow(r,5)) - 2*pow(E,2)*M*(448*pow(M,5) + 1232*pow(M,4)*r + 680*pow(M,3)*pow(r,2) - 1460*pow(M,2)*pow(r,3) - 660*M*pow(r,4) + 149*pow(r,5))))))/(3.*pow(r,7)*pow(pow(a,2) + r*(-2*M + r),3));
	A322 = (-4*pow(r,2)*(2*pow(a,3)*E*L*M + 6*a*E*L*M*pow(r,2) + pow(a,4)*(3*M + 2*r) + pow(r,2)*(pow(r,2)*(-5*M + 2*r) + pow(L,2)*(-7*M + 3*r)) - 2*pow(a,2)*(pow(L,2)*(M - r) + r*(2*pow(M,2) + M*r - 2*pow(r,2)))))/pow(pow(a,2) + r*(-2*M + r),2);
	A324 = (4*pow(a,9)*E*L*M*(2*(-113 + 9*pow(E,2))*pow(M,2) - 57*M*r + 3*pow(r,2)) - 4*a*E*L*M*pow(r,4)*(pow(L,4)*(240*pow(M,2) - 22*M*r - 51*pow(r,2)) + pow(L,2)*pow(r,2)*(314*pow(M,2) + (-83 + 114*pow(E,2))*M*r - 42*pow(r,2)) + 3*pow(r,4)*(-84*pow(M,2) + 4*(26 + 9*pow(E,2))*M*r - 33*pow(r,2))) + pow(a,10)*M*(-21*r*(7*M + 2*r) + 4*pow(E,2)*(113*pow(M,2) + 81*M*r + 9*pow(r,2))) - 4*pow(a,7)*E*L*M*(3*pow(L,2)*(-18*pow(M,2) + 8*M*r + pow(r,2)) + r*(4*(31 + 108*pow(E,2))*pow(M,3) + 9*(109 + 30*pow(E,2))*pow(M,2)*r + (37 + 108*pow(E,2))*M*pow(r,2) - 69*pow(r,3))) - 4*pow(a,3)*E*L*M*pow(r,2)*(pow(L,4)*(612*pow(M,2) - 92*M*r - 45*pow(r,2)) + pow(r,3)*(-2250*pow(M,3) + (955 - 54*pow(E,2))*pow(M,2)*r + 9*(49 + 36*pow(E,2))*M*pow(r,2) - 211*pow(r,3)) + pow(L,2)*r*(16*pow(M,3) + 2*(409 + 114*pow(E,2))*pow(M,2)*r + 2*(-185 + 108*pow(E,2))*M*pow(r,2) - 45*pow(r,3))) + pow(r,4)*(2*pow(L,6)*M*(12*pow(M,2) - 16*M*r + 5*pow(r,2)) + pow(L,2)*pow(r,4)*(304*pow(M,3) + 2*(-243 + 340*pow(E,2))*pow(M,2)*r + (257 - 368*pow(E,2))*M*pow(r,2) - 45*pow(r,3)) + 2*pow(L,4)*pow(r,2)*(128*pow(M,3) + 2*(-113 + 111*pow(E,2))*pow(M,2)*r + (129 - 113*pow(E,2))*M*pow(r,2) - 24*pow(r,3)) + 3*pow(r,6)*(116*pow(M,3) + 2*(-85 + 44*pow(E,2))*pow(M,2)*r + 2*(41 - 27*pow(E,2))*M*pow(r,2) - 13*pow(r,3))) - 4*pow(a,5)*E*L*M*r*(6*pow(L,4)*r + pow(L,2)*(1296*pow(M,3) + 6*(47 + 34*pow(E,2))*pow(M,2)*r + (-371 + 102*pow(E,2))*M*pow(r,2) - 60*pow(r,3)) + r*(-1008*pow(M,4) + 2*(-641 + 216*pow(E,2))*pow(M,3)*r + 6*(335 + 63*pow(E,2))*pow(M,2)*pow(r,2) + (209 + 324*pow(E,2))*M*pow(r,3) - 178*pow(r,4))) + pow(a,2)*pow(r,2)*(2*pow(L,6)*M*(408*pow(M,2) - 194*M*r - 5*pow(r,2)) - 8*pow(L,4)*M*r*(104*pow(M,3) - 3*(73 + 77*pow(E,2))*pow(M,2)*r + (115 - 101*pow(E,2))*M*pow(r,2) + 17*(-1 + 3*pow(E,2))*pow(r,3)) - 3*pow(r,5)*(464*pow(M,4) + 4*(-187 + 194*pow(E,2))*pow(M,3)*r - 8*(-66 + 83*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-83 + 98*pow(E,2))*M*pow(r,3) + 15*pow(r,4)) + pow(L,2)*pow(r,3)*(-4696*pow(M,4) + 8*(672 + 127*pow(E,2))*pow(M,3)*r + 6*(-347 + 324*pow(E,2))*pow(M,2)*pow(r,2) - 3*(-87 + 352*pow(E,2))*M*pow(r,3) + 15*pow(r,4))) + pow(a,4)*r*(12*pow(L,6)*M*r + 3*pow(r,3)*(60*pow(M,5) - 4*(23 + 326*pow(E,2))*pow(M,4)*r + (751 - 436*pow(E,2))*pow(M,3)*pow(r,2) + 2*(-387 + 586*pow(E,2))*pow(M,2)*pow(r,3) - 4*(-50 + 63*pow(E,2))*M*pow(r,4) + 19*pow(r,5)) + 2*pow(L,2)*r*(-1008*pow(M,5) + 4*(-227 + 328*pow(E,2))*pow(M,4)*r + 2*(1325 + 758*pow(E,2))*pow(M,3)*pow(r,2) + 14*(-95 + 6*pow(E,2))*pow(M,2)*pow(r,3) + 3*(39 - 160*pow(E,2))*M*pow(r,4) + 69*pow(r,5)) + 2*pow(L,4)*(864*pow(M,4) + 2*(-193 + 107*pow(E,2))*pow(M,2)*pow(r,2) + 17*(1 - 5*pow(E,2))*M*pow(r,3) + 24*pow(r,4) + 12*pow(M,3)*(r + 102*pow(E,2)*r))) + pow(a,8)*(2*pow(L,2)*M*((226 - 108*pow(E,2))*pow(M,2) + 24*(-2 + pow(E,2))*M*r + 3*(-15 + 8*pow(E,2))*pow(r,2)) + r*(2*pow(E,2)*M*(124*pow(M,3) + 982*pow(M,2)*r + 686*M*pow(r,2) + 3*pow(r,3)) + 3*r*(179*pow(M,3) - 70*pow(M,2)*r + 30*M*pow(r,2) + 10*pow(r,3)))) + pow(a,6)*(-12*pow(L,4)*M*(6*pow(M,2) - 4*M*r - (-3 + pow(E,2))*pow(r,2)) + 2*pow(L,2)*r*(4*(31 + 648*pow(E,2))*pow(M,4) + 28*(35 + 39*pow(E,2))*pow(M,3)*r - (361 + 140*pow(E,2))*pow(M,2)*pow(r,2) + (43 - 112*pow(E,2))*M*pow(r,3) + 39*pow(r,4)) + pow(r,2)*(3*r*(-192*pow(M,4) + 210*pow(M,3)*r - 381*pow(M,2)*pow(r,2) + 160*M*pow(r,3) + 31*pow(r,4)) - 4*pow(E,2)*M*(504*pow(M,4) + 828*pow(M,3)*r - 433*pow(M,2)*pow(r,2) - 709*M*pow(r,3) + 90*pow(r,4)))))/(3.*pow(r,4)*pow(pow(a,2) + r*(-2*M + r),3));
	A340 = -(pow(r,5)*(3*pow(a,2) + r*(-5*M + 2*r)))/(2.*pow(pow(a,2) + r*(-2*M + r),2));
	A342 = (24*pow(a,7)*E*L*M*((-25 + 3*pow(E,2))*M - 7*r) + 3*pow(a,8)*(r*(8*M + 9*r) + 2*pow(E,2)*M*(50*M + 19*r)) - 4*a*E*L*M*pow(r,3)*(3*pow(r,2)*(-76*pow(M,2) + (103 + 18*pow(E,2))*M*r - 36*pow(r,2)) + pow(L,2)*(208*pow(M,2) - 64*M*r - 27*pow(r,2))) - 12*pow(a,5)*E*L*M*(pow(L,2)*(-18*M + 8*r) + r*((4 + 48*pow(E,2))*pow(M,2) + (71 - 6*pow(E,2))*M*r + 20*pow(r,2))) - 12*pow(a,3)*E*L*M*r*(pow(L,2)*(144*pow(M,2) - 64*M*r - pow(r,2)) + 2*r*(-88*pow(M,3) + 2*(-5 + 18*pow(E,2))*pow(M,2)*r + 9*(4 + pow(E,2))*M*pow(r,2) - 13*pow(r,3))) + pow(r,3)*(-8*pow(L,4)*M*(2*pow(M,2) - 3*M*r + pow(r,2)) + pow(r,4)*(348*pow(M,3) + 8*(-64 + 33*pow(E,2))*pow(M,2)*r + (245 - 162*pow(E,2))*M*pow(r,2) - 38*pow(r,3)) + pow(L,2)*pow(r,2)*(-52*pow(M,3) + 4*(19 + 85*pow(E,2))*pow(M,2)*r - (33 + 184*pow(E,2))*M*pow(r,2) + 4*pow(r,3))) - pow(a,2)*r*(-12*pow(L,4)*M*(48*pow(M,2) - 29*M*r + pow(r,2)) + pow(L,2)*r*(1056*pow(M,4) - 16*(33 + 107*pow(E,2))*pow(M,3)*r + 4*(-35 + 16*pow(E,2))*pow(M,2)*pow(r,2) + 4*(9 + 64*pow(E,2))*M*pow(r,3) + 25*pow(r,4)) + pow(r,3)*(696*pow(M,4) + 8*(-103 + 171*pow(E,2))*pow(M,3)*r + 4*(166 - 273*pow(E,2))*pow(M,2)*pow(r,2) + 3*(-111 + 124*pow(E,2))*M*pow(r,3) + 57*pow(r,4))) + pow(a,6)*(pow(L,2)*((300 - 216*pow(E,2))*pow(M,2) + 6*(9 + 8*pow(E,2))*M*r + 27*pow(r,2)) + r*(12*pow(E,2)*M*(2*pow(M,2) + 61*M*r + 15*pow(r,2)) + r*(-312*pow(M,2) - 25*M*r + 61*pow(r,2)))) + pow(a,4)*(24*pow(L,4)*M*(-3*M + 2*r) + pow(L,2)*r*(24*(1 + 72*pow(E,2))*pow(M,3) + 12*(10 - 41*pow(E,2))*pow(M,2)*r - 3*(23 + 8*pow(E,2))*M*pow(r,2) - 2*pow(r,3)) + pow(r,2)*(-12*pow(E,2)*M*(88*pow(M,3) + 64*pow(M,2)*r - 97*M*pow(r,2) + 12*pow(r,3)) + r*(876*pow(M,3) - 512*pow(M,2)*r + 63*M*pow(r,2) + 15*pow(r,3)))))/(12.*r*pow(pow(a,2) + r*(-2*M + r),3));
	A360 = (pow(r,2)*(-144*pow(a,5)*E*L*M + 28*pow(a,3)*E*L*M*(4*M - 3*r)*r + 3*pow(a,6)*(24*pow(E,2)*M + 7*r) + 4*a*E*L*M*pow(r,2)*(68*pow(M,2) - 54*M*r + 15*pow(r,2)) + pow(r,2)*(-2*pow(L,2)*M*(68*pow(M,2) - 76*M*r + 21*pow(r,2)) + pow(r,2)*(116*pow(M,3) + 8*(-21 + 11*pow(E,2))*pow(M,2)*r + (83 - 54*pow(E,2))*M*pow(r,2) - 14*pow(r,3))) + pow(a,4)*(72*pow(L,2)*M + r*(r*(-53*M + 28*r) + pow(E,2)*(-56*pow(M,2) + 74*M*r))) - pow(a,2)*r*(2*pow(L,2)*M*(28*M - 5*r) + r*(r*(36*pow(M,2) - 30*M*r + 7*pow(r,2)) + 4*pow(E,2)*M*(34*pow(M,2) - 16*M*r + 13*pow(r,2))))))/(48.*pow(pow(a,2) + r*(-2*M + r),3));
	A402 = (12*pow(r,3)*(pow(a,2)*(2*M + r) + r*(pow(L,2) + pow(r,2))))/pow(pow(a,2) + r*(-2*M + r),2);
	A404 = (2*(-4*pow(a,7)*E*L*M*(48*pow(M,2) + 40*M*r + 3*pow(r,2)) - 4*pow(a,5)*E*L*M*r*(-96*pow(M,3) - 110*pow(M,2)*r + 73*M*pow(r,2) + 8*pow(r,3) + 9*pow(L,2)*(6*M + r)) + pow(a,8)*(12*pow(E,2)*M*(8*pow(M,2) + 10*M*r + 3*pow(r,2)) - r*(62*pow(M,2) + 42*M*r + 3*pow(r,2))) + pow(r,3)*(4*pow(L,6)*M*(-2*M + r) + 3*pow(r,6)*(23*pow(M,2) + 4*(-5 + pow(E,2))*M*r + 3*pow(r,2)) + 3*pow(L,2)*pow(r,4)*(42*pow(M,2) + 4*(-9 + 2*pow(E,2))*M*r + 5*pow(r,2)) + pow(L,4)*pow(r,2)*(73*pow(M,2) + 4*(-17 + 2*pow(E,2))*M*r + 12*pow(r,2))) - 4*pow(a,3)*E*L*M*pow(r,2)*(6*pow(L,4) + pow(L,2)*(-48*pow(M,2) + 81*M*r + 9*pow(r,2)) + pow(r,2)*(-234*pow(M,2) + 50*M*r + 17*pow(r,2))) + 4*a*E*L*M*pow(r,3)*(pow(L,4)*(4*M - 6*r) + 3*(M - 4*r)*pow(r,4) + pow(L,2)*(-45*M*pow(r,2) + 6*pow(r,3))) + pow(a,2)*pow(r,2)*(12*pow(L,6)*M + 6*pow(r,4)*(-16*pow(M,3) + (26 - 36*pow(E,2))*pow(M,2)*r + (-25 + 12*pow(E,2))*M*pow(r,2) + 6*pow(r,3)) + 3*pow(L,2)*pow(r,2)*(-152*pow(M,3) + 4*(35 + 9*pow(E,2))*pow(M,2)*r + 2*(-33 + 16*pow(E,2))*M*pow(r,2) + 19*pow(r,3)) + pow(L,4)*(-96*pow(M,3) + 4*(39 - 2*pow(E,2))*pow(M,2)*r + 10*(-5 + 2*pow(E,2))*M*pow(r,2) + 24*pow(r,3))) + pow(a,4)*r*(3*pow(L,4)*(36*pow(M,2) + 4*(-1 + pow(E,2))*M*r - pow(r,2)) - 2*pow(L,2)*(96*pow(M,4) + 4*(19 + 12*pow(E,2))*pow(M,3)*r - (143 + 84*pow(E,2))*pow(M,2)*pow(r,2) + 15*(3 - 4*pow(E,2))*M*pow(r,3) - 18*pow(r,4)) - 3*pow(r,2)*(36*pow(M,4) + 4*(-5 + 48*pow(E,2))*pow(M,3)*r + 9*(-5 + 8*pow(E,2))*pow(M,2)*pow(r,2) + 2*(19 - 24*pow(E,2))*M*pow(r,3) - 14*pow(r,4))) + 2*pow(a,6)*(pow(L,2)*(48*pow(M,3) + 2*(10 + 27*pow(E,2))*pow(M,2)*r + 3*(-11 + 8*pow(E,2))*M*pow(r,2) - 3*pow(r,3)) + r*(-12*pow(E,2)*M*(8*pow(M,3) + 12*pow(M,2)*r - 5*M*pow(r,2) - 5*pow(r,3)) + r*(74*pow(M,3) + 17*pow(M,2)*r - 33*M*pow(r,2) + 6*pow(r,3))))))/(3.*pow(r,3)*pow(pow(a,2) + r*(-2*M + r),3));
	A420 = (3*pow(r,6))/pow(pow(a,2) + r*(-2*M + r),2);
	A422 = (2*pow(a,5)*E*L*M*(-8*M + r) - 2*pow(a,3)*E*L*M*r*(4*pow(L,2) - 8*pow(M,2) + M*r - 7*pow(r,2)) - 2*a*E*L*M*pow(r,3)*(4*pow(L,2) + r*(3*M + 4*r)) + pow(a,6)*M*(-7*r + pow(E,2)*(8*M + 6*r)) + pow(r,4)*(pow(r,2)*(23*pow(M,2) + 2*(-11 + 2*pow(E,2))*M*r + 4*pow(r,2)) + pow(L,2)*(25*pow(M,2) + (-23 + 4*pow(E,2))*M*r + 4*pow(r,2))) + pow(a,2)*r*(4*pow(L,4)*M + pow(r,2)*(-16*pow(M,3) - 4*(-5 + 8*pow(E,2))*pow(M,2)*r + 7*(-5 + 2*pow(E,2))*M*pow(r,2) + 14*pow(r,3)) + pow(L,2)*(-8*pow(M,3) + 10*pow(M,2)*r + (-23 + 8*pow(E,2))*M*pow(r,2) + 15*pow(r,3))) + pow(a,4)*(4*pow(L,2)*M*(2*M + (-2 + pow(E,2))*r) + r*(-8*pow(E,2)*M*(pow(M,2) + M*r - 2*pow(r,2)) + r*(17*pow(M,2) - 8*M*r + 10*pow(r,2)))))/pow(pow(a,2) + r*(-2*M + r),3);
	A440 = (pow(r,4)*(9*pow(a,4) - 48*a*E*L*pow(M,2) + 12*pow(L,2)*M*(2*M - r) + pow(r,2)*(69*pow(M,2) + 2*(-31 + 6*pow(E,2))*M*r + 10*pow(r,2)) + 2*pow(a,2)*(6*pow(E,2)*M*(2*M + r) + r*(-30*M + 17*r))))/(24.*pow(pow(a,2) + r*(-2*M + r),3));
	A502 = (-2*pow(r,2)*(2*pow(a,3)*E*L*M + 6*a*E*L*M*pow(r,2) + pow(a,4)*(3*M + 2*r) + pow(r,2)*(pow(r,2)*(-4*M + r) + pow(L,2)*(-6*M + 2*r)) + pow(a,2)*(-2*pow(L,2)*(M - r) + r*(-2*pow(M,2) - 3*M*r + 3*pow(r,2)))))/pow(pow(a,2) + r*(-2*M + r),3);
	A504 = (-2*pow(a,9)*E*L*M*(4*(49 + 27*pow(E,2))*pow(M,2) + 6*(-7 + 12*pow(E,2))*M*r - 9*pow(r,2)) + pow(a,10)*M*(3*r*(-10*M + 7*r) + 2*pow(E,2)*(98*pow(M,2) + 30*M*r - 15*pow(r,2))) - 2*a*E*L*M*pow(r,4)*(pow(L,2)*pow(r,2)*(-725*pow(M,2) + 8*(172 + 21*pow(E,2))*M*r - 570*pow(r,2)) + 3*pow(r,4)*(-135*pow(M,2) + 4*(49 + 18*pow(E,2))*M*r - 82*pow(r,2)) + 16*pow(L,4)*(2*pow(M,2) + 15*M*r - 9*pow(r,2))) - 2*pow(a,7)*E*L*M*(pow(L,2)*(324*pow(M,2) + 12*(5 + 3*pow(E,2))*M*r + 9*pow(r,2)) + 2*r*(16*(2 + 9*pow(E,2))*pow(M,3) + 6*(94 + 57*pow(E,2))*pow(M,2)*r + 2*(-65 + 126*pow(E,2))*M*pow(r,2) - 39*pow(r,3))) + pow(r,4)*(-2*pow(L,6)*M*(44*pow(M,2) - 52*M*r + 15*pow(r,2)) + 3*M*pow(r,6)*(52*pow(M,2) + (-65 + 44*pow(E,2))*M*r + (23 - 30*pow(E,2))*pow(r,2)) + 2*pow(L,4)*pow(r,2)*(-153*pow(M,3) + (271 + 140*pow(E,2))*pow(M,2)*r - 2*(82 + 39*pow(E,2))*M*pow(r,2) + 36*pow(r,3)) + pow(L,2)*pow(r,4)*(34*pow(M,3) + (51 + 428*pow(E,2))*pow(M,2)*r - (109 + 248*pow(E,2))*M*pow(r,2) + 48*pow(r,3))) - 2*pow(a,3)*E*L*M*pow(r,2)*(2*pow(L,4)*(252*pow(M,2) + 116*M*r - 63*pow(r,2)) + pow(L,2)*r*(-784*pow(M,3) + 3*(219 + 40*pow(E,2))*pow(M,2)*r + 18*(39 + 22*pow(E,2))*M*pow(r,2) - 168*pow(r,3)) - 2*pow(r,3)*(765*pow(M,3) + 2*(8 + 135*pow(E,2))*pow(M,2)*r - 3*(89 + 132*pow(E,2))*M*pow(r,2) + 119*pow(r,3))) - 2*pow(a,5)*E*L*M*r*(18*pow(L,4)*(6*M + r) + pow(L,2)*(864*pow(M,3) + 84*(15 + 2*pow(E,2))*pow(M,2)*r + 2*(-145 + 132*pow(E,2))*M*pow(r,2) - 123*pow(r,3)) + r*(-576*pow(M,4) - 1390*pow(M,3)*r + 3*(667 + 108*pow(E,2))*pow(M,2)*pow(r,2) + 4*(89 + 252*pow(E,2))*M*pow(r,3) - 61*pow(r,4))) - pow(a,2)*pow(r,2)*(4*pow(L,6)*M*(-84*pow(M,2) + 8*M*r + 23*pow(r,2)) + 3*pow(r,5)*(210*pow(M,4) + (-355 + 484*pow(E,2))*pow(M,3)*r + 12*(23 - 33*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-62 + 65*pow(E,2))*M*pow(r,3) + 21*pow(r,4)) + 2*pow(L,4)*r*(392*pow(M,4) - (409 + 196*pow(E,2))*pow(M,3)*r - (39 + 584*pow(E,2))*pow(M,2)*pow(r,2) + (-41 + 167*pow(E,2))*M*pow(r,3) + 84*pow(r,4)) + pow(L,2)*pow(r,3)*(1618*pow(M,4) + (-1559 + 2272*pow(E,2))*pow(M,3)*r + (683 - 4036*pow(E,2))*pow(M,2)*pow(r,2) + (-413 + 848*pow(E,2))*M*pow(r,3) + 165*pow(r,4))) + pow(a,4)*r*(18*pow(L,6)*M*(4*M + r) + 8*pow(L,4)*(72*pow(M,4) + 18*(4 + 7*pow(E,2))*pow(M,3)*r + (-91 + 128*pow(E,2))*pow(M,2)*pow(r,2) - 4*(-4 + 5*pow(E,2))*M*pow(r,3) + 9*pow(r,4)) + 3*pow(r,3)*(44*pow(M,5) - 56*(-1 + 8*pow(E,2))*pow(M,4)*r + 2*(233 - 192*pow(E,2))*pow(M,3)*pow(r,2) + (-499 + 736*pow(E,2))*pow(M,2)*pow(r,3) + (151 - 220*pow(E,2))*M*pow(r,4) - 27*pow(r,5)) - pow(L,2)*r*(576*pow(M,5) + 4*(263 + 196*pow(E,2))*pow(M,4)*r - 13*(215 + 88*pow(E,2))*pow(M,3)*pow(r,2) - 6*(-313 + 524*pow(E,2))*pow(M,2)*pow(r,3) + 4*(-82 + 241*pow(E,2))*M*pow(r,4) + 96*pow(r,5))) + pow(a,8)*(pow(L,2)*M*(4*(49 + 162*pow(E,2))*pow(M,2) + 12*(-12 + 23*pow(E,2))*M*r + 3*(11 - 4*pow(E,2))*pow(r,2)) + r*(2*pow(E,2)*M*(32*pow(M,3) + 608*pow(M,2)*r + 286*M*pow(r,2) - 105*pow(r,3)) + 3*r*(24*pow(M,3) + 2*pow(M,2)*r + 130*M*pow(r,2) + 15*pow(r,3)))) + pow(a,6)*(6*pow(L,4)*M*(36*pow(M,2) + 2*(-1 + 18*pow(E,2))*M*r + (5 + 3*pow(E,2))*pow(r,2)) + pow(L,2)*r*(64*(1 + 27*pow(E,2))*pow(M,4) + 16*(65 + 207*pow(E,2))*pow(M,3)*r + 68*(-14 + 17*pow(E,2))*pow(M,2)*pow(r,2) + (595 - 376*pow(E,2))*M*pow(r,3) + 117*pow(r,4)) + pow(r,2)*(3*r*(-58*pow(M,4) - 45*pow(M,3)*r - 420*pow(M,2)*pow(r,2) + 173*M*pow(r,3) + 9*pow(r,4)) - 4*pow(E,2)*M*(144*pow(M,4) + 432*pow(M,3)*r - 218*pow(M,2)*pow(r,2) - 416*M*pow(r,3) + 135*pow(r,4)))))/(3.*pow(r,4)*pow(pow(a,2) + r*(-2*M + r),4));
	A520 = -(pow(r,5)*(3*pow(a,2) + r*(-4*M + r)))/(2.*pow(pow(a,2) + r*(-2*M + r),3));
	A522 = (-6*pow(a,7)*E*L*M*(4*(14 + 3*pow(E,2))*M + 23*r) + 3*pow(a,8)*(2*pow(E,2)*M*(28*M - r) + r*(43*M + 9*r)) - 2*a*E*L*M*pow(r,3)*(pow(L,2)*(64*pow(M,2) + 344*M*r - 222*pow(r,2)) + 3*pow(r,2)*(-199*pow(M,2) + 4*(79 + 18*pow(E,2))*M*r - 134*pow(r,2))) - 6*pow(a,5)*E*L*M*(4*pow(L,2)*(9*M + 2*r) + r*(12*(7 + 6*pow(E,2))*pow(M,2) + 6*(17 + 16*pow(E,2))*M*r + 7*pow(r,2))) - 6*pow(a,3)*E*L*M*r*(2*pow(L,2)*(108*pow(M,2) + 52*M*r - 33*pow(r,2)) + r*(-272*pow(M,3) + (-73 + 72*pow(E,2))*pow(M,2)*r + 2*(103 + 78*pow(E,2))*M*pow(r,2) - 116*pow(r,3))) + pow(r,3)*(-4*pow(L,4)*M*(38*pow(M,2) - 45*M*r + 13*pow(r,2)) + 6*pow(r,4)*(52*pow(M,3) + (-63 + 44*pow(E,2))*pow(M,2)*r - 6*(-3 + 5*pow(E,2))*M*pow(r,2) + 2*pow(r,3)) + pow(L,2)*pow(r,2)*(-158*pow(M,3) + (369 + 428*pow(E,2))*pow(M,2)*r - 4*(73 + 62*pow(E,2))*M*pow(r,2) + 84*pow(r,3))) + pow(a,2)*r*(24*pow(L,4)*M*(18*pow(M,2) + M*r - 7*pow(r,2)) + pow(L,2)*r*(-816*pow(M,4) + 2*(213 + 356*pow(E,2))*pow(M,3)*r + (103 + 1444*pow(E,2))*pow(M,2)*pow(r,2) + 4*(32 - 125*pow(E,2))*M*pow(r,3) - 156*pow(r,4)) - 3*pow(r,3)*(210*pow(M,4) + (-295 + 548*pow(E,2))*pow(M,3)*r + (313 - 488*pow(E,2))*pow(M,2)*pow(r,2) + 2*(-97 + 95*pow(E,2))*M*pow(r,3) + 40*pow(r,4))) + 3*pow(a,6)*(pow(L,2)*(8*(7 + 9*pow(E,2))*pow(M,2) + 8*(6 + pow(E,2))*M*r + 9*pow(r,2)) + r*(2*pow(E,2)*M*(42*pow(M,2) + 148*M*r - 37*pow(r,2)) + r*(-240*pow(M,2) + 87*M*r + 11*pow(r,2)))) + 3*pow(a,4)*(8*pow(L,4)*M*(3*M + r) + pow(L,2)*r*(12*(7 + 36*pow(E,2))*pow(M,3) + 4*(-23 + 98*pow(E,2))*pow(M,2)*r + 4*(3 - 19*pow(E,2))*M*pow(r,2) - 13*pow(r,3)) - pow(r,2)*(r*(-385*pow(M,3) + 317*pow(M,2)*r - 120*M*pow(r,2) + 42*pow(r,3)) + 2*pow(E,2)*M*(136*pow(M,3) + 144*pow(M,2)*r - 312*M*pow(r,2) + 101*pow(r,3)))))/(12.*r*pow(pow(a,2) + r*(-2*M + r),4));
	A540 = (pow(r,2)*(-72*pow(a,5)*E*L*M + 9*pow(a,6)*(4*pow(E,2)*M + r) + 12*pow(a,3)*E*L*M*r*(-12*M + 11*r) + 12*a*E*L*M*pow(r,2)*(32*pow(M,2) - 42*M*r + 17*pow(r,2)) + pow(a,4)*(36*pow(L,2)*M + r*(18*pow(E,2)*M*(4*M - 3*r) + (24*M - 13*r)*r)) + pow(r,2)*(-6*pow(L,2)*M*(32*pow(M,2) - 38*M*r + 11*pow(r,2)) + pow(r,2)*(156*pow(M,3) + (-193 + 132*pow(E,2))*pow(M,2)*r + 2*(32 - 45*pow(E,2))*M*pow(r,2) + 2*pow(r,3))) - pow(a,2)*r*(pow(L,2)*(-72*pow(M,2) + 78*M*r) + r*(12*pow(E,2)*M*(16*pow(M,2) - 23*M*r + 15*pow(r,2)) + r*(183*pow(M,2) - 196*M*r + 62*pow(r,2))))))/(48.*pow(pow(a,2) + r*(-2*M + r),4));
	A600 = pow(r,6)/pow(pow(a,2) + r*(-2*M + r),3);
	A602 = (-24*pow(a,5)*E*L*M*(5*M + r) - 16*pow(a,3)*E*L*M*r*(3*pow(L,2) - 8*pow(M,2) + 2*M*r + 3*pow(r,2)) + 8*a*E*L*M*pow(r,2)*(pow(L,2)*(4*M - 6*r) + (11*M - 15*r)*pow(r,2)) + 3*pow(a,6)*(r*(-6*M + r) + 4*pow(E,2)*M*(5*M + 3*r)) + pow(r,2)*(8*pow(L,4)*M*(-2*M + r) + M*pow(r,4)*(43*M + 4*(-8 + 5*pow(E,2))*r) + pow(L,2)*pow(r,2)*(3*pow(M,2) + 4*(3 + 4*pow(E,2))*M*r - 12*pow(r,2))) + pow(a,2)*r*(24*pow(L,4)*M + pow(r,2)*(-34*pow(M,3) + (39 - 100*pow(E,2))*pow(M,2)*r + 2*(-47 + 38*pow(E,2))*M*pow(r,2) + 36*pow(r,3)) + pow(L,2)*(-64*pow(M,3) - 8*(-3 + 2*pow(E,2))*pow(M,2)*r + 10*(-3 + 4*pow(E,2))*M*pow(r,2) + 48*pow(r,3))) + pow(a,4)*(3*pow(L,2)*(20*pow(M,2) + 4*(-1 + 2*pow(E,2))*M*r + pow(r,2)) + r*(r*(32*pow(M,2) - 14*M*r + 39*pow(r,2)) + pow(E,2)*(-64*pow(M,3) + 8*pow(M,2)*r + 92*M*pow(r,2)))))/(6.*pow(pow(a,2) + r*(-2*M + r),4));
	A620 = (pow(r,3)*(-48*pow(a,3)*E*L*M + 16*a*E*L*M*(M - 3*r)*r + 3*pow(a,4)*(8*pow(E,2)*M + 5*r) + M*r*(pow(L,2)*(-8*M + 4*r) + pow(r,2)*(43*M + 4*(-8 + 5*pow(E,2))*r)) + pow(a,2)*(24*pow(L,2)*M + 2*r*(r*(-31*M + 18*r) + pow(E,2)*(-4*pow(M,2) + 22*M*r)))))/(24.*pow(pow(a,2) + r*(-2*M + r),4));
	A700 = (pow(r,5)*(-pow(a,2) + M*r))/(2.*pow(pow(a,2) + r*(-2*M + r),4));
	A702 = (-12*pow(a,7)*E*L*M*((2 + 6*pow(E,2))*M + 3*r) + 3*pow(a,8)*(2*pow(E,2)*M*(2*M - 7*r) + r*(26*M + 3*r)) - 4*a*E*L*M*pow(r,3)*(pow(L,2)*(-52*pow(M,2) + 144*M*r - 69*pow(r,2)) + pow(r,2)*(-85*pow(M,2) + 2*(61 + 27*pow(E,2))*M*r - 51*pow(r,2))) - 8*pow(a,5)*E*L*M*(27*pow(L,2)*M + r*((35 + 12*pow(E,2))*pow(M,2) + 3*(8 + 17*pow(E,2))*M*r - 9*pow(r,2))) + M*pow(r,3)*(-4*pow(L,4)*(26*pow(M,2) - 31*M*r + 9*pow(r,2)) + pow(L,2)*pow(r,2)*(-25*pow(M,2) + 4*(3 + 44*pow(E,2))*M*r + 2*(5 - 54*pow(E,2))*pow(r,2)) + pow(r,4)*(91*pow(M,2) + 2*(-59 + 40*pow(E,2))*M*r + 2*(23 - 29*pow(E,2))*pow(r,2))) - 4*pow(a,3)*E*L*M*r*(pow(L,2)*(72*pow(M,2) + 172*M*r - 69*pow(r,2)) + r*(-80*pow(M,3) - 65*pow(M,2)*r + 2*(64 + 69*pow(E,2))*M*pow(r,2) - 84*pow(r,3))) + pow(a,2)*r*(4*pow(L,4)*M*(24*pow(M,2) + 35*M*r - 31*pow(r,2)) + pow(L,2)*r*(-160*pow(M,4) - 52*(-1 + 2*pow(E,2))*pow(M,3)*r + (27 + 1004*pow(E,2))*pow(M,2)*pow(r,2) + 2*(51 - 130*pow(E,2))*M*pow(r,3) - 84*pow(r,4)) - pow(r,3)*(166*pow(M,4) + (-289 + 532*pow(E,2))*pow(M,3)*r + (371 - 472*pow(E,2))*pow(M,2)*pow(r,2) + 16*(-17 + 14*pow(E,2))*M*pow(r,3) + 72*pow(r,4))) + pow(a,6)*(3*pow(L,2)*((4 + 72*pow(E,2))*pow(M,2) + 26*M*r + 3*pow(r,2)) + r*(4*pow(E,2)*M*(35*pow(M,2) + 88*M*r - 48*pow(r,2)) + r*(-322*pow(M,2) + 147*M*r - 9*pow(r,2)))) + pow(a,4)*(72*pow(L,4)*pow(M,2) + pow(L,2)*r*(4*(35 + 72*pow(E,2))*pow(M,3) + 4*(-40 + 239*pow(E,2))*pow(M,2)*r + (23 - 152*pow(E,2))*M*pow(r,2) - 36*pow(r,3)) - pow(r,2)*(4*pow(E,2)*M*(40*pow(M,3) + 78*pow(M,2)*r - 195*M*pow(r,2) + 79*pow(r,3)) + r*(-376*pow(M,3) + 365*pow(M,2)*r - 205*M*pow(r,2) + 90*pow(r,3)))))/(12.*r*pow(pow(a,2) + r*(-2*M + r),5));
	A720 = -(pow(r,2)*(-48*pow(a,5)*E*L*M + 4*pow(a,3)*E*L*M*(68*M - 57*r)*r + 3*pow(a,6)*(8*pow(E,2)*M + 3*r) - 4*a*E*L*M*pow(r,2)*(52*pow(M,2) - 98*M*r + 45*pow(r,2)) + M*pow(r,2)*(2*pow(L,2)*(52*pow(M,2) - 64*M*r + 19*pow(r,2)) + pow(r,2)*(-91*pow(M,2) + 2*(59 - 40*pow(E,2))*M*r + 2*(-23 + 29*pow(E,2))*pow(r,2))) + pow(a,4)*(24*pow(L,2)*M + r*(-2*pow(E,2)*M*(68*M - 65*r) + r*(-79*M + 42*r))) + pow(a,2)*r*(2*pow(L,2)*M*(-68*M + 49*r) + r*(4*pow(E,2)*M*(26*pow(M,2) - 66*M*r + 41*pow(r,2)) + r*(187*pow(M,2) - 212*M*r + 72*pow(r,2))))))/(48.*pow(pow(a,2) + r*(-2*M + r),5));
	A800 = (pow(r,3)*(-12*pow(a,3)*E*L*M + 4*a*E*L*M*(2*M - 3*r)*r + 3*pow(a,4)*(2*pow(E,2)*M + r) + M*r*(pow(L,2)*(-4*M + 2*r) + pow(r,2)*(5*M + 4*(-1 + pow(E,2))*r)) + 2*pow(a,2)*(3*pow(L,2)*M + r*(r*(-5*M + 3*r) + pow(E,2)*M*(-2*M + 5*r)))))/(12.*pow(pow(a,2) + r*(-2*M + r),5));
	A900 = -(pow(r,2)*(-18*pow(a,5)*E*L*M + 2*pow(a,3)*E*L*M*(26*M - 21*r)*r + 3*pow(a,6)*(3*pow(E,2)*M + r) - 8*a*E*L*M*pow(r,2)*(2*pow(M,2) - 6*M*r + 3*pow(r,2)) + M*pow(r,2)*(pow(L,2)*(8*pow(M,2) - 10*M*r + 3*pow(r,2)) + 2*pow(r,2)*(-5*pow(M,2) + (7 - 4*pow(E,2))*M*r + 3*(-1 + pow(E,2))*pow(r,2))) + pow(a,4)*(9*pow(L,2)*M + r*(-26*pow(E,2)*M*(M - r) + r*(-17*M + 9*r))) + pow(a,2)*r*(2*pow(L,2)*M*(-13*M + 8*r) + r*(r*(30*pow(M,2) - 35*M*r + 12*pow(r,2)) + pow(E,2)*M*(8*pow(M,2) - 38*M*r + 23*pow(r,2))))))/(24.*pow(pow(a,2) + r*(-2*M + r),6));
  }

  /* Compute aLpha, beta coefficients */
  alpha20 = r*r/(a*a+r*(r-2*M));
  alpha02 = r*r;
  beta = 4.0*(L*L + r*r + a*a*(r+2*M)/r);
}
