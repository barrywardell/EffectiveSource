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
 ******************************************************************************/

struct coordinate {
  double r;
  double theta;
  double phi;
  double t;
};

void effsource_init(double M, double a);
void effsource_set_particle(struct coordinate * x_p, struct coordinate * u_p);
void effsource_set_particle_el(struct coordinate * x_p, double e, double l, double ur_p);
void effsource_phis(struct coordinate * x, double * phis);
void effsource_calc(struct coordinate * x, double *phis, double *dphis_dr, double *dphis_dth,
                    double *dphis_dph, double *dphis_dt, double *box_phis);
void effsource_phis_m(int m, struct coordinate * x, double * phis_re, double * phis_im);
void effsource_calc_m(int m, struct coordinate * x, double *phis, double *dphis_dr, double *dphis_dth,
                    double *dphis_dph, double *dphis_dt, double *box_phis);
double phis_calc(struct coordinate * x);
