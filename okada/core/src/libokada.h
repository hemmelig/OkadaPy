/*
 * =============================================================================
 *
 *       Filename:  libokada.h
 *
 *        Purpose:  Header file for analytical deformation equations, after
 *                  Okada, 1992.
 *
 *      Copyright:  Conor A. Bacon, 2024
 *        License:  GNU General Public License, Version 3
 *                  (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * =============================================================================
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef _OPENMP
    /* Generate a compiler error to stop the build */
    mustLinkOpenMP
#endif

#define PI 3.14159265358979323846

 typedef struct
 {
     double c;
     double s;
     double ss;
     double cc;
     double sc;
     double s2;
     double c2;
 } trigonometric_constants;

 typedef struct
{
    double a1;
    double a2;
    double a3;
    double a4;
    double a5;
} alpha_constants;

typedef struct
{
    double x;
    double y;
    double z;
    double d;
    double p;
    double q;
    double s;
    double t;
    double xy;
    double x2;
    double y2;
    double d2;
    double r2;
    double r;
    double r3;
    double r5;
    double r7;
    double a3;
    double a5;
    double b3;
    double c3;
    double qr;
    double qrx;
    double uy;
    double uz;
    double vy;
    double vz;
    double wy;
    double wz;
} geometric_constants_ps;

typedef struct
{
    double xi;
    double xi2;
    double eta;
    double eta2;
    double q;
    double q2;
    double r;
    double r2;
    double r3;
    double r5;
    double y;
    double d;
    double tt;
    double alx;
    double x11;
    double x32;
    double rxi;
    double y11;
    double y32;
    double ale;
    double reta;
    double ey;
    double ez;
    double fy;
    double fz;
    double gy;
    double gz;
    double hy;
    double hz;
} geometric_constants_frs;

trigonometric_constants compute_trig_functions(double delta);

alpha_constants compute_alpha_constants(double alpha);

void convert2okada_coordinates(
    double x_grid,
    double y_grid,
    double x_min,
    double y_min,
    double x_max,
    double y_max,
    double z_max,
    double z_min,
    double dip,
    double *result
);

void transform_tensor(
    double sin_beta,
    double cos_beta,
    double *s_in,
    double *s_out
);

void printer(double *arr, int n);

geometric_constants_ps compute_geometric_constants_ps(
    double x,
    double y,
    double z,
    double d,
    trigonometric_constants trig_constants
);

/*
 * Function: compute_ps_deformation
 * --------------------------------
 * Compute the internal displacement field due to a point source in a
 * half-space, after the analytical equations presented in Okada, 1992.
 * 
 * alpha: medium constant, calculated from the Lamé's constants.
 * x: x-coordinate of point in space for which to compute displacement.
 * y: y-coordinate of point in space for which to compute displacement.
 * z: z-coordinate of point in space for which to compute displacement.
 * depth: depth of point source.
 * dip: dip of point source.
 * strike_potency: Moment divided by Lamé's mu.
 * dip_potency: Moment divided by Lamé's mu.
 * tensile_potency: Moment divided by Lamé's mu.
 * inflation_potency: Moment divided by Lamé's mu.
 *
 * returns: void - the computed displacement for a given point in space is
 *                 written to pre-specified memory.
 */

void compute_ps_deformation(
    double alpha,
    double x,
    double y,
    double z,
    double depth,
    double dip,
    double strike_potency,
    double dip_potency,
    double tensile_potency,
    double inflation_potency,
    double *u
);

geometric_constants_frs compute_geometric_constants_frs(
    double xi,
    double eta,
    double q,
    trigonometric_constants t
);

void compute_frs_deformation(
    double alpha,
    double x,
    double y,
    double z,
    double depth,
    double dip,
    double fault_half_length,
    double fault_half_width,
    double strike_dislocation,
    double dip_dislocation,
    double tensile_dislocation,
    double *u
);
