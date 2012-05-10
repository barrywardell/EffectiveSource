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
 * Version 3.0 alpha 2 - 10 May 2012
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
#include <gsl/gsl_sf_ellint.h>

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

/* Compute the singular field at the point x for the particle at xp */
void effsource_phis_m(int m, struct coordinate * x, double * phis)
{
  /* For odd m, the circular orbit singular field is exactly 0 */
  if(m % 2 == 1)
  {
    *phis = 0;
    return;
  }

  double A, A0, A2, A4, A6, A8, alpha, ellE, ellK;

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

  A0 = dr6*(A600 + A700*dr) + dr8*(A800 + A900*dr) + dr4*(A420 + A520*dr + dr2*(A620 + A720*dr))*dtheta2 + (A080 + A180*dr)*dtheta8 + dtheta4*(dr2*(A240 + A340*dr) + dr4*(A440 + A540*dr) + (A060 + A160*dr + dr2*(A260 + A360*dr))*dtheta2);
  A2 = dr4*(A402 + A502*dr + dr2*(A602 + A702*dr)) + (dr2*(A222 + A322*dr) + dr4*(A422 + A522*dr))*dtheta2 + dtheta4*(A042 + A142*dr + dr2*(A242 + A342*dr) + (A062 + A162*dr)*dtheta2);
  A4 = dr2*(A204 + A304*dr) + dr4*(A404 + A504*dr) + (A024 + A124*dr + dr2*(A224 + A324*dr))*dtheta2 + (A044 + A144*dr)*dtheta4;
  A6 = A006 + A106*dr + dr2*(A206 + A306*dr) + (A026 + A126*dr)*dtheta2;
  A8 = A008 + A108*dr;

  alpha = alpha20*dr2 + alpha02*dtheta2;

  const double C = alpha / beta;

  ellE = gsl_sf_ellint_Kcomp(sqrt(1.0/(1.0+C)), GSL_PREC_DOUBLE);
  ellK = gsl_sf_ellint_Ecomp(sqrt(1.0/(1.0+C)), GSL_PREC_DOUBLE);

  const double C2  = C*C;
  const double C3  = C2*C;
  const double C4  = C2*C2;
  const double C5  = C3*C2;
  const double C6  = C3*C3;
  const double C7  = C4*C3;
  const double C8  = C4*C4;
  const double C9  = C5*C4;
  const double C10 = C5*C5;
  const double C11 = C6*C5;
  const double C12 = C6*C6;
  const double C13 = C7*C6;
  const double C14 = C7*C7;
  const double C15 = C8*C7;
  const double C16 = C8*C8;

  switch(m) {
  case 0:
    A = (A0*(8 + 23*C + 23*C2) + C*(A2*(2 + 7*C - 3*C2) + C*(A4*(3 - 7*C - 2*C2) + C*(-(A6*(23 + 23*C + 8*C2)) + A8*(15 + 103*C + 128*C2 + 48*C3)))))*ellE - C*(A0*(4 + 8*C) + C*(A2*(1 - 3*C) + C*(-2*A4*(3 + C) - A6*(15 + 19*C + 8*C2) + 4*A8*C*(15 + 26*C + 12*C2))))*ellK;
    A = A / 15.;
    break;
  case 2:
    A = (A0*(8 + 19*C + 9*C2 + 6*C3) + C*(A2*(2 + C + 11*C2 + 4*C3) + C*(A4*(3 + 39*C + 44*C2 + 16*C3) + C*(-(A6*(53 + 229*C + 264*C2 + 96*C3)) + A8*(-5 + 153*C + 614*C2 + 704*C3 + 256*C4)))))*ellE - C*(A0*(4 + 6*C + 6*C2) + C*(A2*(1 + 9*C + 4*C2) + C*(4*A4*(6 + 9*C + 4*C2) - A6*(15 + 139*C + 216*C2 + 96*C3) + 2*A8*C*(25 + 187*C + 288*C2 + 128*C3))))*ellK;
    A = A / 15.;
    break;
  case 4:
    A = (A0*(8 + 7*C - 9*C2 - 32*C3 - 16*C4) + C*(-(A2*(-2 + 17*C + 131*C2 + 168*C3 + 64*C4)) + C*(A4*(3 + 297*C + 1006*C2 + 1088*C3 + 384*C4) + C*(-(A6*(63 + 1047*C + 2976*C2 + 3008*C3 + 1024*C4)) + A8*(-1 + 199*C + 2400*C2 + 6304*C3 + 6144*C4 + 2048*C5)))))*ellE - C*(-4*A0*(-1 + 6*C2 + 4*C3) + C*(-(A2*(-1 + 75*C + 136*C2 + 64*C3)) + C*(2*A4*(57 + 315*C + 448*C2 + 192*C3) + 4*A8*C*(13 + 276*C + 1032*C2 + 1280*C3 + 512*C4) - A6*(15 + 459*C + 1920*C2 + 2496*C3 + 1024*C4))))*ellK;
    A = A / 15.;
    break;
  case 6:
    A = (7*A0*(8 - 13*C + 41*C2 + 454*C3 + 640*C4 + 256*C5) + C*(-7*A2*(-2 + 47*C + 1461*C2 + 4364*C3 + 4480*C4 + 1536*C5) + C*(7*A4*(3 + 807*C + 6156*C2 + 14064*C3 + 12800*C4 + 4096*C5) + C*(-7*A6*(69 + 2661*C + 15288*C2 + 31136*C3 + 26624*C4 + 8192*C5) + A8*(-3 + 1551*C + 41226*C2 + 208320*C3 + 398080*C4 + 327680*C5 + 98304*C6)))))*ellE - C*(14*A0*(2 - 5*C + 123*C2 + 256*C3 + 128*C4) + C*(-7*A2*(-1 + 615*C + 2796*C2 + 3712*C3 + 1536*C4) + C*(28*A4*(66 + 765*C + 2364*C2 + 2688*C3 + 1024*C4) - 7*A6*(15 + 987*C + 8040*C2 + 21408*C3 + 22528*C4 + 8192*C5) + 2*A8*C*(183 + 8037*C + 56352*C2 + 138624*C3 + 139264*C4 + 49152*C5))))*ellK;
    A = A / 105.;
    break;
  case 8:
    A = (21*A0*(8 - 41*C + 279*C2 + 6784*C3 + 18752*C4 + 18432*C5 + 6144*C6) + C*(-21*A2*(-2 + 89*C + 6019*C2 + 33184*C3 + 65152*C4 + 54272*C5 + 16384*C6) + C*(21*A4*(3 + 1593*C + 21950*C2 + 88192*C3 + 149760*C4 + 114688*C5 + 32768*C6) + C*(-3*A6*(513 + 36129*C + 358104*C2 + 1248128*C3 + 1957888*C4 + 1425408*C5 + 393216*C6) + A8*(-5 + 4979*C + 237312*C2 + 2018352*C3 + 6488064*C4 + 9682944*C5 + 6815744*C6 + 1835008*C7)))))*ellE + C*(-84*A0 + 21*(24*A0 - A2)*C - 63*(992*A0 - 737*A2 + 158*A4 - 5*A6)*C2 + (-256704*A0 + 371616*A2 - 198198*A4 + 36207*A6 - 1100*A8)*C3 - 24*(13440*A0 - 39536*A2 + 43904*A4 - 20427*A6 + 3463*A8)*C4 - 48*(2688*A0 - 20160*A2 + 46704*A4 - 45816*A6 + 19801*A8)*C5 + 24576*(14*A2 - 84*A4 + 173*A6 - 158*A8)*C6 - 49152*(14*A4 - 75*A6 + 144*A8)*C7 + 1179648*(A6 - 5*A8)*C8 - 1835008*A8*C9)*ellK;
    A = A / 315.;
    break;
  case 10:
    A = (231*A0*(8 - 77*C + 873*C2 + 37190*C3 + 169600*C4 + 297216*C5 + 229376*C6 + 65536*C7) + C*(-231*A2*(-2 + 143*C + 16949*C2 + 149804*C3 + 478592*C4 + 706048*C5 + 491520*C6 + 131072*C7) + C*(33*A4*(21 + 18705*C + 408724*C2 + 2568656*C3 + 6999552*C4 + 9408512*C5 + 6160384*C6 + 1572864*C7) + C*(-11*A6*(1609 + 180809*C + 2776824*C2 + 14908320*C3 + 37140480*C4 + 47112192*C5 + 29622272*C6 + 7340032*C7) + 7*A8*(-5 + 8217*C + 619622*C2 + 8058432*C3 + 39494400*C4 + 92897280*C5 + 113278976*C6 + 69206016*C7 + 16777216*C8)))))*ellE - C*(924*A0 - 231*(42*A0 - A2)*C + 117440512*A8*C10 + 693*(4850*A0 - 1901*A2 + 248*A4 - 5*A6)*C2 + (21880320*A0 - 16250388*A2 + 5246604*A4 - 616957*A6 + 12110*A8)*C3 + 2*(24393600*A0 - 33190080*A2 + 21283944*A4 - 6292572*A6 + 703381*A8)*C4 + 96*(473088*A0 - 1245552*A2 + 1495824*A4 - 888085*A6 + 249774*A8)*C5 + 768*(19712*A0 - 128128*A2 + 301488*A4 - 337480*A6 + 191135*A8)*C6 - 24576*(1232*A2 - 7216*A4 + 15895*A6 - 17010*A8)*C7 + 98304*(528*A4 - 2904*A6 + 6125*A8)*C8 - 7340032*(11*A6 - 58*A8)*C9)*ellK;
    A = A / 3465.;
    break;
  case 12:
    A = (3003*A0*(8 - 121*C + 2039*C2 + 135392*C3 + 919664*C4 + 2490368*C5 + 3276800*C6 + 2097152*C7 + 524288*C8) + C*(-429*A2*(-14 + 1463*C + 269973*C2 + 3499864*C3 + 16518848*C4 + 37502976*C5 + 44400640*C6 + 26476544*C7 + 6291456*C8) + C*(143*A4*(63 + 85213*C + 2714630*C2 + 24667200*C3 + 98485632*C4 + 202604544*C5 + 224722944*C6 + 127926272*C7 + 29360128*C8) + C*(-13*A6*(18329 + 3015617*C + 66622112*C2 + 512430912*C3 + 1858323456*C4 + 3590553600*C5 + 3810525184*C6 + 2099249152*C7 + 469762048*C8) + 7*A8*(-45 + 110971*C + 12194528*C2 + 226385568*C3 + 1579419648*C4 + 5381167104*C5 + 9956229120*C6 + 10235150336*C7 + 5502926848*C8 + 1207959552*C9)))))*ellE - C*(12012*A0 - 3003*(64*A0 - A2)*C - 469762048*(13*A6 - 73*A8)*C10 + 8455716864*A8*C11 + 9009*(16312*A0 - 4057*A2 + 358*A4 - 5*A6)*C2 + (1377199824*A0 - 642192408*A2 + 140151154*A4 - 11493001*A6 + 157500*A8)*C3 + 16*(292131840*A0 - 237732924*A2 + 100657128*A4 - 20641244*A6 + 1631651*A8)*C4 + 96*(76876800*A0 - 108725760*A2 + 81732508*A4 - 32802614*A6 + 6486949*A8)*C5 + 3072*(1793792*A0 - 4736160*A2 + 6250816*A4 - 4480463*A6 + 1738387*A8)*C6 + 6144*(256256*A0 - 1629056*A2 + 4040608*A4 - 5105360*A6 + 3545479*A8)*C7 - 1572864*(1716*A2 - 10296*A4 + 24518*A6 - 30065*A8)*C8 + 7340032*(572*A4 - 3302*A6 + 7641*A8)*C9)*ellK;
    A = A / 45045.;
    break;
  case 14:
    A = (429*A0*(56 - 1211*C + 28287*C2 + 2715050*C3 + 25687552*C4 + 98860032*C5 + 194281472*C6 + 206503936*C7 + 113246208*C8 + 25165824*C9) + C*(-143*A2*(-42 + 6027*C + 1604761*C2 + 28711820*C3 + 187731456*C4 + 604133376*C5 + 1066500096*C6 + 1054998528*C7 + 549453824*C8 + 117440512*C9) + C*(13*A4*(693 + 1330737*C + 58296436*C2 + 724401200*C3 + 3985422336*C4 + 11580112896*C5 + 19102826496*C6 + 18004574208*C7 + 9042919424*C8 + 1879048192*C9) + C*(7*A8*(-33 + 114469*C + 6442450944*C10 + 17303918*C2 + 435865536*C3 + 4116646656*C4 + 19180781568*C5 + 49832067072*C6 + 75839307776*C7 + 67259858944*C8 + 32212254720*C9) - A6*(245207 + 55657543*C + 1676535864*C2 + 17530109856*C3 + 87229071360*C4 + 237290029056*C5 + 373612085248*C6 + 340115062784*C7 + 166295764992*C8 + 33822867456*C9)))))*ellE - C*(12012*A0 - 3003*(90*A0 - A2)*C + 117440512*(208*A4 - 1272*A6 + 3217*A8)*C10 - 33822867456*(A6 - 6*A8)*C11 + 45097156608*A8*C12 + 9009*(44002*A0 - 7645*A2 + 488*A4 - 5*A6)*C2 + (5038725120*A0 - 1633734388*A2 + 258068460*A4 - 15595561*A6 + 157542*A8)*C3 + 2*(11964446208*A0 - 6593283840*A2 + 1994342584*A4 - 300669300*A6 + 17670793*A8)*C4 + 96*(582506496*A0 - 528189376*A2 + 275229760*A4 - 80007919*A6 + 11732938*A8)*C5 + 768*(89872640*A0 - 136236672*A2 + 117774592*A4 - 59270656*A6 + 16793791*A8)*C6 + 98304*(439296*A0 - 1209780*A2 + 1757964*A4 - 1472881*A6 + 723898*A8)*C7 + 196608*(54912*A0 - 356928*A2 + 945880*A4 - 1335492*A6 + 1094485*A8)*C8 - 7340032*(2288*A2 - 14352*A4 + 37025*A6 - 51198*A8)*C9)*ellK;
    A = A / 45045.;
    break;
  case 16:
    A = (17017*A0*(24 - 699*C + 67108864*C10 + 21573*C2 + 2838016*C3 + 35665664*C4 + 184541184*C5 + 501915648*C6 + 780140544*C7 + 698351616*C8 + 335544320*C9) + C*(-221*A2*(-462 + 87087*C + 7516192768*C10 + 31718581*C2 + 749350272*C3 + 6482574848*C4 + 27995025408*C5 + 68318232576*C6 + 98747547648*C7 + 83898662912*C8 + 38755368960*C9) + C*(17*A4*(9009 + 23394531*C + 135291469824*C10 + 1350145418*C2 + 22013165056*C3 + 159657974784*C4 + 621026770944*C5 + 1413698224128*C6 + 1944273813504*C7 + 1591377657856*C8 + 714038312960*C9) + C*(-17*A6*(251213 + 75317197*C + 180388626432*C10 + 2971493976*C2 + 40602852864*C3 + 265566044160*C4 + 964964745216*C5 + 2092985221120*C6 + 2776440504320*C7 + 2209760673792*C8 + 969588867072*C9) + 7*A8*(-429 + 1997387*C + 3092376453120*C10 + 566935683072*C11 + 398053504*C2 + 13089192048*C3 + 161090961408*C4 + 984295931904*C5 + 3409374609408*C6 + 7136879312896*C7 + 9211631108096*C8 + 7173669126144*C9)))))*ellE - C*(204204*A0 - 51051*(120*A0 - A2)*C - 7516192768*(221*A2 - 1462*A4 + 4080*A6 - 6284*A8)*C10 + 67645734912*(34*A4 - 221*A6 + 608*A8)*C11 - 180388626432*(17*A6 - 109*A8)*C12 + 3968549781504*A8*C13 + 153153*(102272*A0 - 13185*A2 + 638*A4 - 5*A6)*C2 + (258998195456*A0 - 62178353536*A2 + 7452256170*A4 - 345580913*A6 + 2678676*A8)*C3 + 8*(203851135488*A0 - 81863817152*A2 + 18640947712*A4 - 2154290235*A6 + 97737311*A8)*C4 + 816*(6401787392*A0 - 4097144064*A2 + 1577517376*A4 - 348289952*A6 + 39381419*A8)*C5 + 49152*(190590400*A0 - 192479066*A2 + 118901332*A4 - 44570447*A6 + 9643907*A8)*C6 + 98304*(96928832*A0 - 158943200*A2 + 155633538*A4 - 93878573*A6 + 34545959*A8)*C7 + 786432*(6534528*A0 - 19055504*A2 + 5*(6078656*A4 - 5824863*A6 + 3454087*A8))*C8 + 3670016*(311168*A0 - 2107456*A2 + 5991888*A4 - 9368360*A6 + 8835285*A8)*C9)*ellK;
    A = A / 765765.;
    break;
  case 18:
    A = (29393*A0*(264 - 9933*C + 23622320128*C10 + 4294967296*C11 + 391017*C2 + 67761734*C3 + 1091746304*C4 + 7298180096*C5 + 26123403264*C6 + 55228563456*C7 + 71257030656*C8 + 55205429248*C9) + C*(-2261*A2*(-858 + 205491*C + 435939180544*C10 + 77309411328*C11 + 98397201*C2 + 2968271212*C3 + 32822001152*C4 + 182873069568*C5 + 586838802432*C6 + 1153052049408*C7 + 1411412656128*C8 + 1050908033024*C9) + C*(323*A4*(9009 + 30494333*C + 4148938407936*C10 + 721554505728*C11 + 2243705956*C2 + 46495090416*C3 + 429916870656*C4 + 2153776250880*C5 + 6438918684672*C6 + 12026429046784*C7 + 14170959380480*C8 + 10242691694592*C9) + C*(-19*A6*(4360711 + 1670533303*C + 92719753986048*C10 + 15874199126016*C11 + 83675183608*C2 + 1448982507168*C3 + 12056898207744*C4 + 56335156936704*C5 + 160272581132288*C6 + 288449641840640*C7 + 330230505930752*C8 + 233253768265728*C9) + 7*A8*(-6435 + 38810927*C + 822565694078976*C10 + 322019467984896*C11 + 54425825574912*C12 + 9869213386*C2 + 411059370816*C3 + 6398235161856*C4 + 49652914618368*C5 + 220854870540288*C6 + 605729541062656*C7 + 1059716826398720*C8 + 1186158298005504*C9)))))*ellE - C*(3879876*A0 - 969969*(154*A0 - A2)*C + 117440512*(1074944*A0 + 5*(-1529728*A2 + 4666704*A4 - 8015416*A6 + 8543857*A8))*C10 - 33822867456*(5168*A2 - 36176*A4 + 108889*A6 - 184950*A8)*C11 + 45097156608*(5168*A4 - 35720*A6 + 106383*A8)*C12 - 15874199126016*(19*A6 - 130*A8)*C13 + 380980779024384*A8*C14 + 2909907*(213282*A0 - 21277*A2 + 808*A4 - 5*A6)*C2 + (12939978082304*A0 - 2401452695964*A2 + 226093771268*A4 - 8298287867*A6 + 50900850*A8)*C3 + 2*(52007456883200*A0 - 15980517631232*A2 + 2844155629032*A4 - 260048074364*A6 + 9377358011*A8)*C4 + 96*(4526382794752*A0 - 2171848589120*A2 + 645634208960*A4 - 112095493333*A6 + 10069259438*A8)*C5 + 14592*(72142377216*A0 - 52760038016*A2 + 24627236864*A4 - 7166256256*A6 + 1224045179*A8)*C6 + 98304*(15621086208*A0 - 17477231548*A2 + 12481023620*A4 - 5717138095*A6 + 1639010198*A8)*C7 + 196608*(6768384896*A0 - 12020705984*A2 + 13168051080*A4 - 9245744204*A6 + 4177657295*A8)*C8 + 7340032*(85995520*A0 - 266994384*A2 + 465207856*A4 - 501780595*A6 + 347885946*A8)*C9)*ellK;
    A = A / 14549535.;
    break;
  case 20:
    A = ellE*(7759752*A0 - 969969*(377*A0 - 2*A2)*C
          + 969969*(18423*A0 - 593*A2 + 3*A4)*C2
          + (3955749733600*A0 - 351164499539*A2 + 12474319641*A4 - 84385039*A6 - 36465*A8)*C3
          + (79492913586800*A0 - 13174868723688*A2 + 1140297608590*A4 - 40239516839*A6 + 277060439*A8)*C4
          + 32*(20818955759616*A0 - 5664762945002*A2 + 915152814058*A4 - 78006324927*A6 + 2737771187*A8)*C5
          + 672*(4504930881536*A0 - 1882222478848*A2 + 500415837996*A4 - 79548030030*A6 + 6713391837*A8)*C6
          + 3072*(2702408554496*A0 - 1676268138720*A2 + 685895857696*A4 - 179536095105*A6 + 28225526310*A8)*C7
          + 6144*(2342338118400*A0 - 2135179582592*A2 + 1299542459488*A4 - 523980136752*A6 + 135630820905*A8)*C8
          + 1048576*(15164979200*A0 - 20569099306*A2 + 18426693810*A4 - 11060962898*A6 + 4411424763*A8)*C9
          + 1048576*(10367752192*A0 - 21837570048*A2 + 29146764332*A4 - 25772784478*A6 + 15308490613*A8)*C10
          + 67108864*(62512128*A0 - 226100000*A2 + 469143136*A4 - 618503655*A6 + 541396502*A8)*C11
          + 402653184*(1736448*A0 - 14180992*A2 + 50574048*A4 - 103719280*A6 + 135416585*A8)*C12
          - 1443109011456*(646*A2 - 5206*A4 + 18361*A6 - 37305*A8)*C13
          + 2267742732288*(532*A4 - 4242*A6 + 14827*A8)*C14
          - 72567767433216*(21*A6 - 166*A8)*C15 + 1886761953263616*A8*C16)
        - ellK*C*(3879876*A0 - 969969*(192*A0 - A2)*C
          + 2909907*(409400*A0 - 32601*A2 + 998*A4 - 5*A6)*C2
          + (30564801234800*A0 - 4525891759816*A2 + 343787826774*A4 - 10234141787*A6 + 50905140*A8)*C3
          + 16*(19055537249280*A0 - 4639500959380*A2 + 663973941176*A4 - 49234940364*A6 + 1444459537*A8)*C4
          + 96*(16706271282176*A0 - 6264545247744*A2 + 1484881335028*A4 - 208219197170*A6 + 15215085719*A8)*C5
          + 3072*(1629853824256*A0 - 910624225504*A2 + 333971614144*A4 - 77786876945*A6 + 10761139805*A8)*C6
          + 6144*(1588951212800*A0 - 1309738188416*A2 + 717976998368*A4 - 259305985520*A6 + 59657231945*A8)*C7
          + 1572864*(7640371200*A0 - 9408174748*A2 + 7626148712*A4 - 4123583706*A6 + 1472639635*A8)*C8
          + 7340032*(1237012480*A0 - 2374799360*A2 + 2880615156*A4 - 2306105554*A6 + 1234117467*A8)*C9
          + 67108864*(57302784*A0 - 189634592*A2 + 359095744*A4 - 430622535*A6 + 341442827*A8)*C10
          + 1207959552*(578816*A0 - 4341120*A2 + 14185248*A4 - 26577040*A6 + 31585235*A8)*C11
          - 103079215104*(9044*A2 - 67032*A4 + 216860*A6 - 402865*A8)*C12
          + 2267742732288*(532*A4 - 3906*A6 + 12535*A8)*C13
          - 217703302299648*(7*A6 - 51*A8)*C14 + 1886761953263616*A8*C15);
    A = A / 14549535.;
    break;
  default:
    A=0;
  }

  *phis = 4.0/(beta*C3*pow(alpha+beta, 2.5));
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
