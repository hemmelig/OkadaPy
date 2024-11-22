/*
 * =============================================================================
 *
 *       Filename:  finite_rectangular_source.c
 *
 *        Purpose:  Compute strain/stress solution to analytical equations
 *                  for a finite rectangular source embedded in an infinite
 *                  elastic halfspace.
 *
 *      Copyright:  Conor A. Bacon, 2024
 *        License:  GNU General Public License, Version 3
 *                  (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * =============================================================================
 */

#include "libokada.h"

/*
 * Function: compute_ua
 * --------------------
 * Compute the Somigliani tensor, which represents the infinite-medium terms
 * of the internal displacement field due to a finite rectangular source in a
 * half-space, after the analytical equations presented in Okada, 1992.
 * 
 * g: Struct containing pre-computed, position-dependent terms.
 * t: Struct containing pre-computed, dip-dependent terms.
 * a: Struct containing pre-computed, material-dependent terms.
 * strike_potency: Moment divided by Lamé's mu.
 * dip_potency: Moment divided by Lamé's mu.
 * tensile_potency: Moment divided by Lamé's mu.
 * inflation_potency: Moment divided by Lamé's mu.
 * u: Pointer to an array in which to write the computed results.
 *
 * returns: void - the computed displacement for a given point in space is
 *                 written to pre-specified memory.
 */
void compute_ua(
    geometric_constants_frs g,
    trigonometric_constants t,
    alpha_constants a,
    double strike_dislocation,
    double dip_dislocation,
    double tensile_dislocation,
    double *u
)
{
    double du[12] = {0};

    double xy = g.xi * g.y11;
    double qx = g.q * g.x11;
    double qy = g.q * g.y11;

    /* ================================ */
    /* === Strike-slip contribution === */
    /* ================================ */
    if (strike_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Ua components)
        du[0] = g.tt / 2. + a.a2 * g.xi * qy;
        du[1] = a.a2 * g.q / g.r;
        du[2] = a.a1 * g.ale - a.a2 * g.q * qy;

        // x-derivatives - Table 7, Okada 1992 (Ua components)
        du[3] = -a.a1 * qy - a.a2 * g.xi2 * g.q * g.y32;
        du[4] = -a.a2 * g.xi * g.q / g.r3;
        du[5] = a.a1 * xy + a.a2 * g.xi * g.q2 * g.y32;

        // y-derivatives - Table 8, Okada 1992 (Ua components)
        du[6] = a.a1 * xy * t.s + a.a2 * g.xi * g.fy + g.d / 2. * g.x11;
        du[7] = a.a2 * g.ey;
        du[8] = a.a1 * (t.c / g.r + qy * t.s) - a.a2 * g.q * g.fy;

        // z-derivatives - Table 9, Okada 1992 (Ua components)
        du[9]  = a.a1 * xy * t.c + a.a2 * g.xi * g.fz + g.y / 2. * g.x11;
        du[10] = a.a2 * g.ez;
        du[11] = -a.a1 * (t.s / g.r - qy * t.c) - a.a2 * g.q * g.fz;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * strike_dislocation / (2.* PI);
        }
    }
    /* ================================ */

    /* ============================= */
    /* === Dip-slip contribution === */
    /* ============================= */
    if (dip_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Ua components)
        du[0] = a.a2 * g.q / g.r;
        du[1] = g.tt / 2. + a.a2 * g.eta * qx;
        du[2] = a.a1 * g.alx - a.a2 * g.q * qx;

        // x-derivatives - Table 7, Okada 1992 (Ua components)
        du[3] = -a.a2 * g.xi * g.q / g.r3;
        du[4] = -qy / 2. - a.a2 * g.eta * g.q / g.r3;
        du[5] = a.a1 / g.r + a.a2 * g.q2 / g.r3;

        // y-derivatives - Table 8, Okada 1992 (Ua components)
        du[6] = a.a2 * g.ey;
        du[7] = a.a1 * g.d * g.x11 + xy / 2. * t.s + a.a2 * g.eta * g.gy;
        du[8] = a.a1 * g.y * g.x11 - a.a2 * g.q * g.gy;

        // z-derivatives - Table 9, Okada 1992 (Ua components)
        du[9]  = a.a2 * g.ez;
        du[10] = a.a1 * g.y * g.x11 + xy / 2. * t.c + a.a2 * g.eta * g.gz;
        du[11] = -a.a1 * g.d * g.x11 - a.a2 * g.q * g.gz;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * dip_dislocation / (2.* PI);
        }
    }
    /* ============================= */

    /* ================================== */
    /* === Tensile fault contribution === */
    /* ================================== */
    if (tensile_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Ua components)
        du[0] = -a.a1 * g.ale - a.a2 * g.q * qy;
        du[1] = -a.a1 * g.alx - a.a2 * g.q * qx;
        du[2] = g.tt / 2. - a.a2 * (g.eta * qx + g.xi * qy);

        // x-derivatives - Table 7, Okada 1992 (Ua components)
        du[3] = -a.a1 * xy + a.a2 * g.xi * g.q2 * g.y32;
        du[4] = -a.a1 / g.r + a.a2 * g.q2 / g.r3;
        du[5] = -a.a1 * qy - a.a2 * g.q * g.q2 * g.y32;

        // y-derivatives - Table 8, Okada 1992 (Ua components)
        du[6] = -a.a1 * (t.c / g.r + qy * t.s) - a.a2 * g.q * g.fy;
        du[7] = -a.a1 * g.y * g.x11 - a.a2 * g.q * g.gy;
        du[8] =  a.a1 * (g.d * g.x11 + xy * t.s) + a.a2 * g.q * g.hy;

        // z-derivatives - Table 9, Okada 1992 (Ua components)
        du[9]  = a.a1 * (t.s / g.r - qy * t.c) - a.a2 * g.q * g.fz;
        du[10] = a.a1 * g.d * g.x11 - a.a2 * g.q * g.gz;
        du[11] = a.a1 * (g.y * g.x11 + xy * t.c) + a.a2 * g.q * g.hz;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * tensile_dislocation / (2.* PI);
        }
    }
    /* ================================== */

    return;
}

/*
 * Function: compute_ub
 * --------------------
 * Compute the surface-deformation related term of the internal displacement
 * field due to a finite rectangular source in a half-space, after the
 * analytical equations presented in Okada, 1992.
 * 
 * g: Struct containing pre-computed, position-dependent terms.
 * t: Struct containing pre-computed, dip-dependent terms.
 * a: Struct containing pre-computed, material-dependent terms.
 * strike_potency: Moment divided by Lamé's mu.
 * dip_potency: Moment divided by Lamé's mu.
 * tensile_potency: Moment divided by Lamé's mu.
 * inflation_potency: Moment divided by Lamé's mu.
 * u: Pointer to an array in which to write the computed results.
 *
 * returns: void - the computed displacement for a given point in space is
 *                 written to pre-specified memory.
 */
void compute_ub(
    geometric_constants_frs g,
    trigonometric_constants t,
    alpha_constants a,
    double strike_dislocation,
    double dip_dislocation,
    double tensile_dislocation,
    double *u
)
{
    double du[12] = {0};
    double c, cc;

    double rd = g.r + g.d;
    double d11 = 1. / (g.r * rd);
    double aj2 = g.xi * g.y / rd * d11;
    double aj5 = -(g.d + g.y * g.y / rd) * d11;

    double x = sqrt(g.xi2 + g.q2);
    double rd2 = rd * rd;

    int c1 = !(t.c <= 1e-12);
    int c2 = (t.c <= 1e-12);
    if (c1) {
        c = t.c;
        cc = t.cc; 
    } else {
        c = 1e-12;
        cc = 1e-12; 
    }

    int s1 = (fabs(g.xi) <= 1e-12);
    int s2 = !(fabs(g.xi) <= 1e-12);

    double ai4 = c1 * (s1 * 0. + s2 * (1. / cc * (g.xi / rd * t.sc 
            + 2. * atan((g.eta * (x + g.q * c)
            + x * (g.r + x) * t.s) / (g.xi * (g.r + x) * c)))))
        + c2 * (g.xi * g.y / rd2 / 2.);
    double ai3 = c1 * ((g.y * c / rd - g.ale + t.s * log(rd)) / cc)
        + c2 * ((g.eta / rd + g.y * g.q / rd2 - g.ale) / 2.);
    double ak1 = c1 * (g.xi * (d11 - g.y11 * t.s) / c) + c2 * (g.xi * g.q / rd * d11);
    double ak3 = c1 * ((g.q * g.y11 - g.y * d11) / c) + c2 * (t.s / rd * (g.xi2 * d11 - 1.));
    double aj3 = c1 * ((ak1 - aj2 * t.s) / c) + c2 * (-g.xi / rd2 * (g.q2 * d11 - 0.5));
    double aj6 = c1 * ((ak3 - aj5 * t.s) / c) + c2 * (-g.y / rd2 * (g.xi2 * d11 - 0.5));

    double xy = g.xi * g.y11;
    double ai1 = -g.xi / rd * t.c - ai4 * t.s;
    double ai2 = log(rd) + ai3 * t.s;
    double ak2 = 1. / g.r + ak3 * t.s;
    double ak4 = xy * t.c - ak1 * t.s;
    double aj1 = aj5 * t.c - aj6 * t.s;
    double aj4 = -xy - aj2 * t.c + aj3 * t.s;

    double qx = g.q * g.x11;
    double qy = g.q * g.y11;

    /* ================================ */
    /* === Strike-slip contribution === */
    /* ================================ */
    if (strike_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Ub components)
        du[0] = -g.xi * qy - g.tt - a.a3 * ai1 * t.s;
        du[1] = -g.q / g.r + a.a3 * g.y / rd * t.s;
        du[2] =  g.q * qy - a.a3 * ai2 * t.s;

        // x-derivatives - Table 7, Okada 1992 (Ub components)
        du[3] = g.xi2 * g.q * g.y32 - a.a3 * aj1 * t.s;
        du[4] = g.xi * g.q / g.r3 - a.a3 * aj2 * t.s;
        du[5] = -g.xi * g.q2 * g.y32 - a.a3 * aj3 * t.s;

        // y-derivatives - Table 8, Okada 1992 (Ub components)
        du[6] = -g.xi * g.fy - g.d * g.x11 + a.a3 * (xy + aj4) * t.s;
        du[7] = -g.ey + a.a3 * (1. / g.r + aj5) * t.s;
        du[8] = g.q * g.fy - a.a3 * (qy - aj6) * t.s;

        // z-derivatives - Table 9, Okada 1992 (Ub components)
        du[9] = -g.xi * g.fz - g.y * g.x11 + a.a3 * ak1 * t.s;
        du[10] = -g.ez + a.a3 * g.y * d11 * t.s;
        du[11] = g.q * g.fz + a.a3 * ak2 * t.s;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * strike_dislocation / (2.* PI);
        }
    }
    /* ================================ */

    /* ============================= */
    /* === Dip-slip contribution === */
    /* ============================= */
    if (dip_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Ub components)
        du[0] = -g.q / g.r + a.a3 * ai3 * t.sc;
        du[1] = -g.eta * qx - g.tt - a.a3 * g.xi / rd * t.sc;
        du[2] =  g.q * qx + a.a3 * ai4 * t.sc;

        // x-derivatives - Table 7, Okada 1992 (Ub components)
        du[3] = g.xi * g.q / g.r3 + a.a3 * aj4 * t.sc;
        du[4] = g.eta * g.q / g.r3 + qy + a.a3 * aj5 * t.sc;
        du[5] = -g.q2 / g.r3 + a.a3 * aj6 * t.sc;

        // y-derivatives - Table 8, Okada 1992 (Ub components)
        du[6] = -g.ey + a.a3 * aj1 * t.sc;
        du[7] = -g.eta * g.gy - xy * t.s + a.a3 * aj2 * t.sc;
        du[8] = g.q * g.gy + a.a3 * aj3 * t.sc;

        // z-derivatives - Table 9, Okada 1992 (Ub components)
        du[9] = -g.ez - a.a3 * ak3 * t.sc;
        du[10] = -g.eta * g.gz - xy * t.c - a.a3 * g.xi * d11 * t.sc;
        du[11] = g.q * g.gz - a.a3 * ak4 * t.sc;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * dip_dislocation / (2.* PI);
        }
    }
    /* ============================= */

    /* ================================== */
    /* === Tensile fault contribution === */
    /* ================================== */
    if (tensile_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Ub components)
        du[0] = g.q * qy - a.a3 * ai3 * t.ss;
        du[1] = g.q * qx + a.a3 * g.xi / rd * t.ss;
        du[2] = g.eta * qx + g.xi * qy - g.tt - a.a3 * ai4 * t.ss;

        // x-derivatives - Table 7, Okada 1992 (Ub components)
        du[3] = -g.xi * g.q2 * g.y32 - a.a3 * aj4 * t.ss;
        du[4] = -g.q2 / g.r3 - a.a3 * aj5 * t.ss;
        du[5] =  g.q * g.q2 * g.y32 - a.a3 * aj6 * t.ss;

        // y-derivatives - Table 8, Okada 1992 (Ub components)
        du[6] =  g.q * g.fy - a.a3 * aj1 * t.ss;
        du[7] =  g.q * g.gy - a.a3 * aj2 * t.ss;
        du[8] = -g.q * g.hy - a.a3 * aj3 * t.ss;

        // z-derivatives - Table 9, Okada 1992 (Ub components)
        du[9] =   g.q * g.fz + a.a3 * ak3 * t.ss;
        du[10] =  g.q * g.gz + a.a3 * g.xi * d11 * t.ss;
        du[11] = -g.q * g.hz + a.a3 * ak4 * t.ss;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * tensile_dislocation / (2.* PI);
        }
    }
    /* ================================== */

    return;
}

/*
 * Function: compute_uc
 * --------------------
 * Compute the depth-multipled term of the internal displacement field due to a
 * finite rectangular source in a half-space, after the analytical equations
 * presented in Okada, 1992.
 * 
 * g: Struct containing pre-computed, position-dependent terms.
 * t: Struct containing pre-computed, dip-dependent terms.
 * a: Struct containing pre-computed, material-dependent terms.
 * strike_potency: Moment divided by Lamé's mu.
 * dip_potency: Moment divided by Lamé's mu.
 * tensile_potency: Moment divided by Lamé's mu.
 * inflation_potency: Moment divided by Lamé's mu.
 * u: Pointer to an array in which to write the computed results.
 *
 * returns: void - the computed displacement for a given point in space is
 *                 written to pre-specified memory.
 */
void compute_uc(
    geometric_constants_frs g,
    trigonometric_constants t,
    alpha_constants a,
    double z,
    double strike_dislocation,
    double dip_dislocation,
    double tensile_dislocation,
    double *u
)
{
    double du[12] = {0};

    double c = g.d + z;

    double x53 = (8.0 * g.r2 + 9.0 * g.r * g.xi + 3. * g.xi2) * g.x11 * g.x11 * g.x11 / g.r2;
    double y53 = (8.0 * g.r2 + 9.0 * g.r * g.eta + 3. * g.eta2) * g.y11 * g.y11 * g.y11 / g.r2;

    double h = g.q * t.c - z;
    double z32 = t.s / g.r3 - h * g.y32;
    double z53 = 3. * t.s / g.r5 - h * y53;

    double y0 = g.y11 - g.xi2 * g.y32;
    double z0 = z32 - g.xi2 * z53;

    double ppy = t.c / g.r3 + g.q * g.y32 * t.s;
    double ppz = t.s / g.r3 - g.q * g.y32 * t.c;

    double qq = z * g.y32 + z32 + z0;
    double qqy = 3. * c * g.d / g.r5 - qq * t.s;
    double qqz = 3. * c * g.y / g.r5 - qq * t.c + g.q * g.y32;
    
    double xy = g.xi * g.y11;
    double qx = g.q * g.x11;
    double qy = g.q * g.y11;
    double qr = 3. * g.q / g.r5;
    
    double cqx = c * g.q * x53;
    double cdr = (c + g.d) / g.r3;

    double yy0 = g.y / g.r3 - y0 * t.c;

    /* ================================ */
    /* === Strike-slip contribution === */
    /* ================================ */
    if (strike_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Uc components)
        du[0] = a.a4 * xy * t.c - a.a5 * g.xi * g.q * z32;
        du[1] = a.a4 * (t.c / g.r + 2. * qy * t.s) - a.a5 * c * g.q / g.r3;
        du[2] = a.a4 * qy * t.c - a.a5 * (c * g.eta / g.r3 - z * g.y11 + g.xi2 * z32);

        // x-derivatives - Table 7, Okada 1992 (Uc components)
        du[3] =  a.a4 * y0 * t.c - a.a5 * g.q * z0;
        du[4] = -a.a4 * g.xi * (t.c / g.r3 + 2. * g.q * g.y32 * t.s) + a.a5 * c * g.xi * qr;
        du[5] = -a.a4 * g.xi * g.q * g.y32 * t.c + a.a5 * g.xi * (3. * c * g.eta / g.r5 - qq);

        // y-derivatives - Table 8, Okada 1992 (Uc components)
        du[6] = -a.a4 * g.xi * ppy * t.c - a.a5 * g.xi * qqy;
        du[7] =  a.a4 * 2. * (g.d / g.r3 - y0 * t.s) * t.s - g.y / g.r3 * t.c - a.a5 * (cdr * t.s - g.eta / g.r3 - c * g.y * qr);
        du[8] = -a.a4 * g.q / g.r3 + yy0 * t.s + a.a5 * (cdr * t.c + c * g.d * qr - (y0 * t.c + g.q * z0) * t.s);

        // z-derivatives - Table 9, Okada 1992 (Uc components)
        du[9] = a.a4 * g.xi * ppz * t.c - a.a5 * g.xi * qqz;
        du[10] = a.a4 * 2. * (g.y / g.r3 - y0 * t.c) * t.s + g.d / g.r3 * t.c - a.a5 * (cdr * t.c + c * g.d * qr);
        du[11] = yy0 * t.c - a.a5 * (cdr * t.s - c * g.y * qr - y0 * t.ss + g.q * z0 * t.c);

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * strike_dislocation / (2.* PI);
        }
    }
    /* ================================ */

    /* ============================= */
    /* === Dip-slip contribution === */
    /* ============================= */
    if (dip_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Uc components)
        du[0] = a.a4 * t.c / g.r - qy * t.s - a.a5 * c * g.q / g.r3;
        du[1] = a.a4 * g.y * g.x11 - a.a5 * c * g.eta * g.q * g.x32;
        du[2] = -g.d * g.x11 - xy * t.s - a.a5 * c * (g.x11 - g.q2 * g.x32);

        // x-derivatives - Table 7, Okada 1992 (Uc components)
        du[3] = -a.a4 * g.xi / g.r3 * t.c + a.a5 * c * g.xi * qr + g.xi * g.q * g.y32 * t.s;
        du[4] = -a.a4 * g.y / g.r3 + a.a5 * c * g.eta * qr;
        du[5] = g.d / g.r3 - y0 * t.s + a.a5 * c / g.r3 * (1. - 3. * g.q2 / g.r2);

        // y-derivatives - Table 8, Okada 1992 (Uc components)
        du[6] = -a.a4 * g.eta / g.r3 + y0 * t.ss - a.a5 * (cdr * t.s - c * g.y * qr);
        du[7] =  a.a4 * (g.x11 - g.y * g.y * g.x32) - a.a5 * c * ((g.d + 2. * g.q * t.c) * g.x32 - g.y * g.eta * g.q * x53);
        du[8] = g.xi * ppy * t.s + g.y * g.d * g.x32 + a.a5 * c *((g.y + 2. * g.q * t.s) * g.x32 - g.y * g.q2 * x53);

        // z-derivatives - Table 9, Okada 1992 (Uc components)
        du[9] = -g.q / g.r3 + y0 * t.sc - a.a5 * (cdr * t.c + c * g.d * qr);
        du[10] = a.a4 * g.y * g.d * g.x32 - a.a5 * c * ((g.y - 2. * g.q * t.s) * g.x32 + g.d * g.eta * g.q * x53);
        du[11] = -g.xi * ppz * t.s + g.x11 - g.d * g.d * g.x32 - a.a5 * c * ((g.d - 2. * g.q * t.c) * g.x32 - g.d * g.q2 * x53);

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * dip_dislocation / (2.* PI);
        }
    }
    /* ============================= */

    /* ================================== */
    /* === Tensile fault contribution === */
    /* ================================== */
    if (tensile_dislocation != 0.)
    {
        // Displacements - Table 6, Okada 1992 (Uc components)
        du[0] = -a.a4 * (t.s / g.r + qy * t.c) - a.a5 * (z * g.y11 - g.q2 * z32);
        du[1] =  a.a4 * 2. * xy * t.s + g.d * g.x11 - a.a5 * c * (g.x11 - g.q2 * g.x32);
        du[2] =  a.a4 * (g.y * g.x11 + xy * t.c) + a.a5 * g.q * (c * g.eta * g.x32 + g.xi * z32);

        // x-derivatives - Table 7, Okada 1992 (Uc components)
        du[3] =  a.a4 * g.xi / g.r3 * t.s + g.xi * g.q * g.y32 * t.c + a.a5 * g.xi * (3. * c * g.eta / g.r5 - 2. * z32 - z0);
        du[4] =  a.a4 * 2. * y0 * t.s - g.d / g.r3 + a.a5 * c / g.r3 * (1. - 3. * g.q2 / g.r2);
        du[5] = -a.a4 * yy0 - a.a5 * (c * g.eta * qr - g.q * z0);

        // y-derivatives - Table 8, Okada 1992 (Uc components)
        du[6] =  a.a4 * (g.q / g.r3 + y0 * t.sc) + a.a5 * (z / g.r3 * t.c + c * g.d * qr - g.q * z0 * t.s);
        du[7] = -a.a4 * 2. * g.xi * ppy * t.s - g.y * g.d * g.x32 + a.a5 * c * ((g.y + 2. * g.q * t.s) * g.x32 - g.y * g.q2 * x53);
        du[8] = -a.a4 * (g.xi * ppy * t.c - g.x11 + g.y * g.y * g.x32) + a.a5 * (c * ((g.d + 2. * g.q * t.c) * g.x32 - g.y * g.eta * g.q * x53) + g.xi * qqy);

        // z-derivatives - Table 9, Okada 1992 (Uc components)
        du[9] = -g.eta / g.r3 + y0 * t.cc - a.a5 * (z / g.r3 * t.s - c * g.y * qr - y0 * t.ss + g.q * z0 * t.c);
        du[10] = a.a4 * 2. * g.xi * ppz * t.s - g.x11 + g.d * g.d * g.x32 - a.a5 * c * ((g.d - 2. * g.q * t.c) * g.x32 - g.d * g.q2 * x53);
        du[11] = a.a4 * (g.xi * ppz * t.c + g.y * g.d * g.x32) + a.a5 * (c * ((g.y - 2. * g.q * t.s) * g.x32 + g.d * g.eta * g.q * x53) + g.xi * qqz);

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * tensile_dislocation / (2.* PI);
        }
    }
    /* ================================== */

    return; 
}

/*
 * Function: compute_frs_deformation
 * --------------------------------
 * Compute the internal displacement field due to a finite rectangular source
 * in a half-space, after the analytical equations presented in Okada, 1992.
 * 
 * alpha: medium constant, calculated from the Lamé's constants.
 * x: x-coordinate of point in space for which to compute displacement.
 * y: y-coordinate of point in space for which to compute displacement.
 * z: z-coordinate of point in space for which to compute displacement.
 * depth: depth of centre of finite rectangular source.
 * dip: dip of finite rectangular source.
 * fault_half_length: Half-length of the fault.
 * fault_half_width: Half-length of the fault.
 * strike_dislocation: Amount of slip along strike.
 * dip_dislocation: Amount of slip down dip.
 * tensile_dislocation: Amount of tensile opening.
 * u: Pointer to an array in which to write the computed results.
 *
 * returns: void - the computed displacement for a given point in space is
 *                 written to pre-specified memory.
 */
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
)
{
    // === Variable declaration ===
    double ua[12] = {0}, ua_image[12] = {0}, ub[12] = {0}, uc[12] = {0};
    double du[12] = {0};
    geometric_constants_frs geometric_constants_frs;
    trigonometric_constants trigonometric_constants;
    alpha_constants alpha_constants;
    double xi, eta, effective_depth, p, q;

    // Compute basic constants
    trigonometric_constants = compute_trig_functions(dip);
    alpha_constants = compute_alpha_constants(alpha);
    effective_depth = depth + z;
    p = y * trigonometric_constants.c
        + effective_depth * trigonometric_constants.s;
    q = y * trigonometric_constants.s
        - effective_depth * trigonometric_constants.c;

    for (int k = 0; k < 2; k++)
    {
        if (k == 0)
        {
            eta = p + fault_half_width;
        }
        else
        {
            eta = p - fault_half_width;
        }
        for (int j = 0; j < 2; j++)
        {
            if (j == 0)
            {
                xi = x + fault_half_length;
            }
            else
            {
                xi = x - fault_half_length;
            }

            // Pre-compute some more useful terms
            geometric_constants_frs = compute_geometric_constants_frs(
                xi,
                eta,
                q,
                trigonometric_constants
            );

            // Detect singular point

            // === Compute real-source Ua terms ===
            compute_ua(
                geometric_constants_frs,
                trigonometric_constants,
                alpha_constants,
                strike_dislocation,
                dip_dislocation,
                tensile_dislocation,
                ua
            );

            for (int i = 0; i <= 9; i += 3)
            {
                du[i]   = -ua[i];
                du[i+1] = -ua[i+1] * trigonometric_constants.c
                    + ua[i+2] * trigonometric_constants.s;
                du[i+2] = -ua[i+1] * trigonometric_constants.s
                    - ua[i+2] * trigonometric_constants.c;
            }
            du[9]  = -du[9];
            du[10] = -du[10];
            du[11] = -du[11];

            for (int i = 0; i < 12; i++)
            {
                if (j + k != 1)
                {
                    u[i] += du[i];
                }
                else
                {
                    u[i] -= du[i];
                }
                ua[i] = 0.;
            }
        }
    }

    effective_depth = depth - z;
    p = y * trigonometric_constants.c
        + effective_depth * trigonometric_constants.s;
    q = y * trigonometric_constants.s
        - effective_depth * trigonometric_constants.c;

    for (int k = 0; k < 2; k++)
    {
        if (k == 0)
        {
            eta = p + fault_half_width;
        }
        else
        {
            eta = p - fault_half_width;
        }
        for (int j = 0; j < 2; j++)
        {
            if (j == 0)
            {
                xi = x + fault_half_length;
            }
            else
            {
                xi = x - fault_half_length;
            }

            // Pre-compute some more useful terms
            geometric_constants_frs = compute_geometric_constants_frs(
                xi,
                eta,
                q,
                trigonometric_constants
            );

            // Detect singular point

            // === Compute image-source Ua terms ===
            compute_ua(
                geometric_constants_frs,
                trigonometric_constants,
                alpha_constants,
                strike_dislocation,
                dip_dislocation,
                tensile_dislocation,
                ua_image
            );

            // === Compute real-source Ub terms ===
            compute_ub(
                geometric_constants_frs,
                trigonometric_constants,
                alpha_constants,
                strike_dislocation,
                dip_dislocation,
                tensile_dislocation,
                ub
            );

            // === Compute real-source Uc terms ===
            compute_uc(
                geometric_constants_frs,
                trigonometric_constants,
                alpha_constants,
                z,
                strike_dislocation,
                dip_dislocation,
                tensile_dislocation,
                uc
            );

            for (int i = 0; i <= 9; i += 3)
            {
                du[i] = ua_image[i] + ub[i] + z * uc[i];
                du[i+1] = (ua_image[i+1] + ub[i+1]
                    + z * uc[i+1]) * trigonometric_constants.c
                    - (ua_image[i+2] + ub[i+2]
                        + z * uc[i+2]) * trigonometric_constants.s;
                du[i+2] = (ua_image[i+1] + ub[i+1]
                    - z * uc[i+1]) * trigonometric_constants.s
                    + (ua_image[i+2] + ub[i+2]
                        - z * uc[i+2]) * trigonometric_constants.c;
            }
            du[9] += uc[0];
            du[10] += uc[1] * trigonometric_constants.c
                - uc[2] * trigonometric_constants.s;
            du[11] -= uc[1] * trigonometric_constants.s
                + uc[2] * trigonometric_constants.c;

            for (int i = 0; i < 12; i++)
            {
                if (j + k != 1)
                {
                    u[i] += du[i];
                }
                else
                {
                    u[i] -= du[i];
                }
                ua_image[i] = 0.;
                ub[i] = 0;
                uc[i] = 0;
            }
        }
    }
}
