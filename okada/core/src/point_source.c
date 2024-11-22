/*
 * =============================================================================
 *
 *       Filename:  point_source.c
 *
 *        Purpose:  Routines for computing the internal displacement field due
 *                  to a point source.
 *
 *      Copyright:  Conor A. Bacon, 2024
 *        License:  GNU General Public License, Version 3
 *                  (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * =============================================================================
 */

#include "libokada.h"

/*
 * Function: compute_ua0
 * ---------------------
 * Compute the Somigliani tensor, which represents the infinite-medium terms
 * of the internal displacement field due to a point source in a half-space,
 * after the analytical equations presented in Okada, 1992.
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
void compute_ua0(
    geometric_constants_ps g,
    trigonometric_constants t,
    alpha_constants a,
    double strike_potency,
    double dip_potency,
    double tensile_potency,
    double inflation_potency,
    double *u
)
{
    double du[12] = {0};

    /* ================================ */
    /* === Strike-slip contribution === */
    /* ================================ */
    if (strike_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Ua components)
        du[0] =  a.a1 * g.q / g.r3 + a.a2 * g.x2 * g.qr;
        du[1] =  a.a1 * g.x / g.r3 * t.s + a.a2 * g.xy * g.qr;
        du[2] = -a.a1 * g.x / g.r3 * t.c + a.a2 * g.x * g.d * g.qr;

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] =  g.x * g.qr * (-a.a1 + (a.a2 * (1 + g.a5)));
        du[4] =  a.a1 * g.a3 / g.r3 * t.s + a.a2 * g.y * g.qr * g.a5;
        du[5] = -a.a1 * g.a3 / g.r3 * t.c + a.a2 * g.d * g.qr * g.a5;

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] = a.a1 * (t.s / g.r3 - g.y * g.qr)
            + a.a2 * 3. * g.x2 / g.r5 * g.uy;
        du[7] = 3. * g.x / g.r5 * (
            -a.a1 * g.y * t.s + a.a2 * (g.y * g.uy + g.q)
        );
        du[8] = 3. * g.x / g.r5 * (
             a.a1 * g.y * t.c + a.a2 * g.d * g.uy
        );

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] = a.a1 * (t.c / g.r3 + g.d * g.qr)
            + a.a2 * 3. * g.x2 / g.r5 * g.uz;
        du[10] = 3. * g.x / g.r5 * (
             a.a1 * g.d * t.s + a.a2 * g.y * g.uz
        );
        du[11] = 3. * g.x / g.r5 * (
            -a.a1 * g.d * t.c + a.a2 * (g.d * g.uz - g.q)
        );

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * strike_potency / (2.* PI);
        }
    }
    /* ================================ */

    /* ============================= */
    /* === Dip-slip contribution === */
    /* ============================= */
    if (dip_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Ua components)
        du[0] =  a.a2 * g.x * g.p * g.qr;
        du[1] =  a.a1 * g.s / g.r3 + a.a2 * g.y * g.p * g.qr;
        du[2] = -a.a1 * g.t / g.r3 + a.a2 * g.d * g.p * g.qr; 

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] =  a.a2 * g.p * g.qr * g.a5;
        du[4] = -a.a1 * 3. * g.x * g.s / g.r5 - a.a2 * g.y * g.p * g.qrx;
        du[5] =  a.a1 * 3. * g.x * g.t / g.r5 - a.a2 * g.d * g.p * g.qrx;

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] =  a.a2 * 3. * g.x / g.r5 * g.vy;
        du[7] =  a.a1 * (t.s2 / g.r3 - 3. * g.y * g.s / g.r5)
            + a.a2 * (3. * g.y / g.r5 * g.vy + g.p * g.qr);
        du[8] = -a.a1 * (t.c2 / g.r3 - 3. * g.y * g.t / g.r5)
            + a.a2 * 3. * g.d / g.r5 * g.vy;

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] =  a.a2 * 3. * g.x / g.r5 * g.vz;
        du[10] = a.a1 * (t.c2 / g.r3 + 3. * g.d * g.s / g.r5)
            + a.a2 * 3. * g.y / g.r5 * g.vz;
        du[11] = a.a1 * (t.s2 / g.r3 - 3. * g.d * g.t / g.r5)
            + a.a2 * (3. * g.d / g.r5 * g.vz - g.p * g.qr);

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * dip_potency / (2.* PI);
        }
    }
    /* ============================= */

    /* ================================== */
    /* === Tensile fault contribution === */
    /* ================================== */
    if (tensile_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Ua components)
        du[0] = a.a1 * g.x / g.r3 - a.a2 * g.x * g.q * g.qr;
        du[1] = a.a1 * g.t / g.r3 - a.a2 * g.y * g.q * g.qr;
        du[2] = a.a1 * g.s / g.r3 - a.a2 * g.d * g.q * g.qr;

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] =  a.a1 * g.a3 /g.r3 - a.a2 * g.q * g.qr * g.a5;
        du[4] = -a.a1 * 3. * g.x * g.t / g.r5 + a.a2 * g.y * g.q * g.qrx;
        du[5] = -a.a1 * 3. * g.x * g.s / g.r5 + a.a2 * g.d * g.q * g.qrx; 

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] = -a.a1 * 3. * g.xy / g.r5 - a.a2 * g.x * g.qr * g.wy;
        du[7] =  a.a1 * (t.c2 / g.r3 - 3. * g.y * g.t / g.r5)
            - a.a2 * (g.y * g.wy + g.q) * g.qr;
        du[8] =  a.a2 * (t.s2 / g.r3 - 3. * g.y * g.s / g.r5)
            - a.a2 * g.d * g.qr * g.wy;

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] = a.a1 * 3 * g.x * g.d / g.r5 - a.a2 * g.x * g.qr * g.wz;
        du[10] = -a.a1 * (t.s2 / g.r3 - 3. * g.d * g.t / g.r5)
            - a.a2 * g.y * g.qr * g.wz;
        du[11] = -a.a1 * (t.c2 / g.r3 + 3. * g.d * g.s / g.r5)
            - a.a2 * (g.d * g.wz - g.q) * g.qr;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * tensile_potency / (2.* PI);
        }
    }
    /* ================================== */

    /* ===================================== */
    /* === Inflation source contribution === */
    /* ===================================== */
    if (inflation_potency != 0)
    {
        // Displacements - Table 2, Okada 1992 (Ua components)
        du[0] = -a.a1 * g.x / g.r3;
        du[1] = -a.a1 * g.y / g.r3;
        du[2] = -a.a1 * g.d / g.r3;

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] = -a.a1 * g.a3 / g.r3;
        du[4] =  a.a1 * 3. * g.xy / g.r5;
        du[5] =  a.a1 * 3. * g.x * g.d / g.r5;

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] = du[4];
        du[7] = -a.a1 * g.b3 / g.r3;
        du[8] =  a.a1 * 3. * g.y * g.d / g.r5;

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] = -du[5];
        du[10] = -du[8];
        du[11] = a.a1 * g.c3 / g.r3;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * inflation_potency / (2.* PI);
        }
    }
    /* ===================================== */

    return;
}

/*
 * Function: compute_ub0
 * ---------------------
 * Compute the surface-deformation related term of the internal displacement
 * field due to a point source in a half-space, after the analytical equations
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
void compute_ub0(
    geometric_constants_ps g,
    trigonometric_constants t,
    alpha_constants a,
    double strike_potency,
    double dip_potency,
    double tensile_potency,
    double inflation_potency,
    double *u
)
{
    double du[12] = {0};

    double c = g.d + g.z;
    double rd = g.r + g.d;
    double d12, d32, d33, d53, d54;
    double fi1, fi2, fi3, fi4, fi5, fj1, fj2, fj3, fj4, fk1, fk2, fk3;

    d12 = 1. / (g.r * rd * rd);
    d32 = d12 * (2. * g.r + g.d) / g.r2;
    d33 = d12 * (3. * g.r + g.d) / (g.r2 * rd);
    d53 = d12 * (8. * g.r2 + 9. * g.r * g.d + 3. * g.d2) / (g.r2 * g.r2 * rd);
    d54 = d12 * (5. * g.r2 + 4. * g.r * g.d + g.d2) / g.r3 * d12;

    fi1 = g.y * (d12 - g.x2 * d33);
    fi2 = g.x * (d12 - g.y2 * d33);
    fi3 = g.x / g.r3 - fi2;
    fi4 = -g.xy * d32;
    fi5 = 1. / (g.r * rd) - g.x2 * d32;
    fj1 = -3. * g.xy * (d33 - g.x2 * d54);
    fj2 = 1. / g.r3 - 3. * d12 + 3. * g.x2 * g.y2 * d54;
    fj3 = g.a3 / g.r3 - fj2;
    fj4 = -3. * g.xy / g.r5 - fj1;
    fk1 = -g.y * (d32 - g.x2 * d53);
    fk2 = -g.x * (d32 - g.y2 * d53);
    fk3 = -3. * g.x * g.d / g.r5 - fk2;

    /* ================================ */
    /* === Strike-slip contribution === */
    /* ================================ */
    if (strike_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Ub components)
        du[0] = -g.x2 * g.qr - a.a3 * fi1 * t.s;
        du[1] = -g.xy * g.qr - a.a3 * fi2 * t.s;
        du[2] = -c * g.x * g.qr - a.a3 * fi4 * t.s;

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] = -g.x * g.qr * (1. + g.a5) - a.a3 * fj1 * t.s;
        du[4] = -g.y * g.qr * g.a5 - a.a3 * fj2 * t.s;
        du[5] = -c * g.qr * g.a5 - a.a3 * fk1 * t.s;

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] = -3. * g.x2 / g.r5 * g.uy - a.a3 * fj2 * t.s;
        du[7] = -3. * g.xy / g.r5 * g.uy - g.x * g.qr - a.a3 * fj4 * t.s;
        du[8] = -3. * c * g.x / g.r5 * g.uy - a.a3 * fk2 * t.s;

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] = -3. * g.x2 / g.r5 * g.uz + a.a3 * fk1 * t.s;
        du[10] = -3. * g.xy * g.r5 * g.uz - a.a3 * fk2 * t.s;
        du[11] = 3. * g.x / g.r5 * (-c * g.uz + a.a3 * g.y * t.s);

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * strike_potency / (2.* PI);
        }
    }
    /* ================================ */

    /* ============================= */
    /* === Dip-slip contribution === */
    /* ============================= */
    if (dip_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Ub components)
        du[0] = -g.x * g.p * g.qr + a.a3 * fi3 * t.sc;
        du[1] = -g.y * g.p * g.qr + a.a3 * fi1 * t.sc;
        du[2] = -c * g.p * g.qr + a.a3 * fi5 * t.sc;

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] = -g.p * g.qr * g.a5 + a.a3 * fj3 * t.sc;
        du[4] =  g.y * g.p * g.qrx + a.a3 * fj1 * t.sc;
        du[5] =  c * g.p * g.qrx + a.a3 * fk3 * t.sc;

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] = -3. * g.x / g.r5 * g.vy + a.a3 * fj1 * t.sc;
        du[7] = -3. * g.y / g.r5 * g.vy - g.p * g.qr + a.a3 * fj2 * t.sc;
        du[8] = -3. * c / g.r5 * g.vy + a.a3 * fk1 * t.sc;

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] =  -3. * g.x / g.r5 * g.vz - a.a3 * fk3 * t.sc;
        du[10] = -3. * g.y / g.r5 * g.vz - a.a3 * fk1 * t.sc;
        du[11] = -3. * c / g.r5 * g.vz + a.a3 * g.a3 / g.r3 * t.sc;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * dip_potency / (2.* PI);
        }
    }
    /* ============================= */

    /* ================================== */
    /* === Tensile fault contribution === */
    /* ================================== */
    if (tensile_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Ub components)
        du[0] = g.x * g.q * g.qr - a.a3 * fi3 * t.ss;
        du[1] = g.y * g.q * g.qr - a.a3 * fi1 * t.ss;
        du[2] = c * g.q * g.qr - a.a3 * fi5 * t.ss;

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] =  g.q * g.qr * g.a5 - a.a3 * fj3 * t.ss;
        du[4] = -g.y * g.q * g.qrx - a.a3 * fj1 * t.ss;
        du[5] = -c * g.q * g.qrx - a.a3 * fk3 * t.ss;

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] = g.x * g.qr * g.wy - a.a3 * fj1 * t.ss;
        du[7] = g.qr * (g.y * g.wy + g.q) - a.a3 * fj2 * t.ss;
        du[8] = c * g.qr * g.wy - a.a3 * fk3 * t.ss;

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] =  g.x * g.qr * g.wz * a.a3 * fk3 * t.ss;
        du[10] = g.y * g.qr * g.wz + a.a3 * fk1 * t.ss;
        du[11] = c * g.qr * g.wz - a.a3 * g.a3 / g.r3 * t.ss;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * tensile_potency / (2.* PI);
        }
    }
    /* ================================== */

    /* ===================================== */
    /* === Inflation source contribution === */
    /* ===================================== */
    if (inflation_potency != 0)
    {
        // Displacements - Table 2, Okada 1992 (Ub components)
        du[0] = a.a3 * g.x / g.r3;
        du[1] = a.a3 * g.y / g.r3;
        du[2] = a.a3 * g.d / g.r3;

        // x-derivatives - Table 3, Okada 1992 (Ua components)
        du[3] =  a.a3 * g.a3 / g.r3;
        du[4] = -a.a3 * 3. * g.xy / g.r5;
        du[5] = -a.a3 * 3. * g.x * g.d / g.r5;

        // y-derivatives - Table 4, Okada 1992 (Ua components)
        du[6] = du[4];
        du[7] =  a.a3 * g.b3 / g.r3;
        du[8] = -a.a3 * 3. * g.y * g.d / g.r5;

        // z-derivatives - Table 5, Okada 1992 (Ua components)
        du[9] =  -du[5];
        du[10] = -du[8];
        du[11] = -a.a3 * g.c3 / g.r3;

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * inflation_potency / (2.* PI);
        }
    }
    /* ===================================== */

    return;
}

/*
 * Function: compute_uc0
 * ---------------------
 * Compute the depth-multipled term of the internal displacement field due to a
 * point source in a half-space, after the analytical equations presented in
 * Okada, 1992.
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
void compute_uc0(
    geometric_constants_ps g,
    trigonometric_constants t,
    alpha_constants a,
    double strike_potency,
    double dip_potency,
    double tensile_potency,
    double inflation_potency,
    double *u
)
{
    double du[12] = {0};

    double c = g.d + g.z;
    double q2 = g.q * g.q;
    double a7 = 1. - 7. * g.x2 / g.r2;
    double b5 = 1. - 5. * g.y2 / g.r2;
    double b7 = 1. - 7. * g.y2 / g.r2;
    double c5 = 1. - 5. * g.d2 / g.r2;
    double c7 = 1. - 7. * g.d2 / g.r2;
    double d7 = 2. - 7. * q2 / g.r2;
    double qr5 = 5. * g.q / g.r2;
    double qr7 = 7. * g.q / g.r2;
    double dr5 = 5. * g.d / g.r2;

    /* ================================ */
    /* === Strike-slip contribution === */
    /* ================================ */
    if (strike_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Uc components)
        du[0] = -a.a4 * g.a3 / g.r3 * t.c + a.a5 * c * g.qr * g.a5;
        du[1] = 3. * g.x / g.r5 * (
             a.a4 * g.y * t.c + a.a5 * c * (t.s - g.y * qr5)
        );
        du[2] = 3. * g.x / g.r5 * (
            -a.a4 * g.y * t.s + a.a5 * c * (t.c + g.d * qr5)
        );

        // x-derivatives - Table 3, Okada 1992 (Uc components)
        du[3] = a.a4 * 3. * g.x / g.r5 * (2. + g.a5) * t.c
            - a.a5 * c * g.qrx * (2. + a7);
        du[4] = 3. / g.r5 * ( a.a4 * g.y * g.a5 * t.c
            + a.a5 * c * (g.a5 * t.s - g.y * qr5 * a7));
        du[5] = 3. / g.r5 * (-a.a4 * g.y * g.a5 + t.s
            + a.a5 * c * (g.a5 * t.c + g.d * qr5 * a7));

        // y-derivatives - Table 4, Okada 1992 (Uc components)
        du[6] = du[4];
        du[7] = 3. * g.x / g.r5 * ( a.a4 * b5 * t.c
            - a.a5 * 5. * c / g.r2 * (2. * g.y * t.s + g.q * b7)
        );
        du[8] = 3. * g.x / g.r5 * (-a.a4 * b5 * t.s
            + a.a5 * 5. * c / g.r2 * (g.d * b7 * t.s - g.y * c7 * t.c)
        );

        // z-derivatives - Table 5, Okada 1992 (Uc components)
        du[9] = 3. / g.r5 * (-a.a4 * g.d * g.a5 * t.c
            + a.a5 * c * (g.a5 * t.c + g.d * qr5 * a7)
        );
        du[10] = 15. * g.x / g.r7 * ( a.a4 * g.y * g.d * t.c
            + a.a5 * c * (g.d * b7 * t.s - g.y * c7 * t.c)
        );
        du[11] = 15. * g.x / g.r7 * (-a.a4 * g.y * g.d * t.s
            + a.a5 * c * (2. * g.d * t.c - g.q * c7)
        );

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * strike_potency / (2.* PI);
        }
    }
    /* ================================ */

    /* ============================= */
    /* === Dip-slip contribution === */
    /* ============================= */
    if (dip_potency != 0.)
    {
        // Displacements - Table 2, Okada 1992 (Uc components)
        du[0] =  a.a4 * 3. * g.x * g.t / g.r5 - a.a5 * c * g.p * g.qrx;
        du[1] = -a.a4 / g.r3 * (t.c2 - 3. * g.y * g.t / g.r2)
            + a.a5 * 3. * c / g.r5 * (g.s - g.y * g.t / g.r2);
        du[2] = -a.a4 * g.a3 / g.r3 * t.sc
            + a.a5 * 3. * c / g.r5 * (g.t + g.d * g.p * qr5);

        // x-derivatives - Table 3, Okada 1992 (Uc components)
        du[3] = a.a4 * 3. * g.t / g.r5 * g.a5
            - a.a5 * 5. * c * g.p * g.qr / g.r2 * a7;
        du[4] = 3. * g.x / g.r5 * (a.a4 * (t.c2 - 5. * g.y * g.t / g.r2)
            - a.a5 * 5. * c / g.r2 * (g.s - g.y * g.p * qr7)
        );
        du[5] = 3. * g.x / g.r5 * (a.a4 * (2. * g.a5) * t.sc
            - a.a5 * 5. * c / g.r2 * (g.t + g.d * g.p * qr7)
        );

        // y-derivatives - Table 4, Okada 1992 (Uc components)
        du[6] = du[4];
        du[7] = 3. / g.r5 * (a.a4 * (2. * g.y * t.c2 + g.t * b5)
            + a.a5 * c * (t.s2 - 10. * g.y * g.s / g.r2 - g.p * qr5 * b7)
        );
        du[8] = 3. / g.r5 * (a.a4 * g.y * g.a5 * t.sc
            - a.a5 * c * ((3. * g.a5) * t.c2 + g.y * g.p * dr5 * qr7)
        );

        // z-derivatives - Table 5, Okada 1992 (Uc components)
        du[9] = 3. * g.x / g.r5 * (-a.a4 * (t.s2 - g.t * dr5)
            - a.a5 * 5. * c / g.r2 * (g.t + g.d * g.p * qr7)
        );
        du[10] = 3. / g.r5 * (-a.a4 * (g.d * b5 * t.c2 + g.y * c5 * t.s2)
            - a.a5 * c * ((3. + g.a5) * t.c2 + g.y * g.p * dr5 * qr7)
        );
        du[11] = 3. / g.r5 * (-a.a4 * g.d * g.a5 * t.sc
            - a.a5 * c * (t.s2 - 10. * g.d * g.t / g.r2 + g.p * qr5 * c7)
        );

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * dip_potency / (2.* PI);
        }
    }
    /* ============================= */

    /* ================================== */
    /* === Tensile fault contribution === */
    /* ================================== */
    if (tensile_potency != 0.)
    {
         // Displacements - Table 2, Okada 1992 (Uc components)
        du[0] = 3. * g.x / g.r5 * (-a.a4 * g.s + a.a5 * (c * g.q * qr5 - g.z));
        du[1] =  a.a4 / g.r3 * (t.s2 - 3. * g.y * g.s / g.r2)
            + a.a5 * 3. / g.r5 * (c * (g.t - g.y + g.y * g.q * qr5) - g.y * g.z);
        du[2] = -a.a4 / g.r3 * (1. - g.a3 * t.ss)
            - a.a5 * 3. / g.r5 * (c * (g.s - g.d + g.d * g.q * qr5) - g.d * g.z);

        // x-derivatives - Table 3, Okada 1992 (Uc components)
        du[3] = a.a4 * 3. * g.s / g.r5 * g.a5
            + a.a5 * (c * g.qr * qr5 * a7 - 3. * g.z / g.r5 * g.a5);
        du[4] = 3. * g.x / g.r5 * (-a.a4 * (t.s2 - 5. * g.y * g.s / g.r2)
            - a.a5 * 5. / g.r2 * (c * (g.t - g.y + g.y * g.q * qr7) - g.y * g.z)
        );
        du[5] = 3. * g.x / g.r5 * ( a.a4 * (1. - (2. + g.a5) * t.ss)
            + a.a5 * 5. / g.r2 * (c * (g.s - g.d + g.d * g.q * qr7) - g.d * g.z)
        );

        // y-derivatives - Table 4, Okada 1992 (Uc components)
        du[6] = du[4];
        du[7] = 3. / g.r5 * (-a.a4 * (2. * g.y * t.s2 + g.s * b5)
            - a.a5 * (c * (2. * t.ss + 10. * g.y * (g.t - g.y) / g.r2
                - g.q * qr5 * b7
            ) + g.z * b5)
        );
        du[8] = 3. / g.r5 * ( a.a4 * g.y * (1 / - g.a5 * t.ss)
            + a.a5 * (c * (3. + g.a5) * t.s2 - g.y * dr5 * (c * d7 + g.z))
        );

        // z-derivatives - Table 5, Okada 1992 (Uc components)
        du[9] = 3. * g.x / g.r5 * (-a.a4 * (t.c2 + g.s * dr5)
            + a.a5 * (5. * c / g.r2 * (g.s - g.d + g.d * g.q * qr7)
                - 1. - g.z * dr5
            )
        );
        du[10] = 3. / g.r5 * ( a.a4 * (g.d * b5 * t.s2 - g.y * c5 * t.c2)
            + a.a5 * (c * ((3. + g.a5) * t.s2
                - g.y * dr5 * d7
            ) - g.y * (1. + g.z * dr5))
        );
        du[11] = 3. / g.r5 * (-a.a4 * g.d * (1. - g.a5 * t.ss)
            - a.a5 * (c * (t.c2 + 10. * g.d * (g.s - g.d) / g.r2
                - g.q * qr5 * c7
            ) + g.z * (1. + c5))
        );

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * tensile_potency / (2.* PI);
        }
    }
    /* ================================== */

    /* ===================================== */
    /* === Inflation source contribution === */
    /* ===================================== */
    if (inflation_potency != 0)
    {
        // Displacements - Table 2, Okada 1992 (Uc components)
        du[0] = a.a4 * 3. * g.x * g.d / g.r5;
        du[1] = a.a4 * 3. * g.y * g.d / g.r5;
        du[2] = a.a4 * g.c3 / g.r3;

        // x-derivatives - Table 3, Okada 1992 (Uc components)
        du[3] =  a.a4 * 3. * g.d / g.r5 * g.a5;
        du[4] = -a.a4 * 15. * g.xy * g.d / g.r7;
        du[5] = -a.a4 * 3. * g.x / g.r5 * c5;

        // y-derivatives - Table 4, Okada 1992 (Uc components)
        du[6] =  du[4];
        du[7] =  a.a4 * 3. * g.d / g.r5 * b5;
        du[8] = -a.a4 * 3. * g.y / g.r5 * c5;

        // z-derivatives - Table 5, Okada 1992 (Uc components)
        du[9] =  du[5];
        du[10] = du[8];
        du[11] = a.a4 * 3. * g.d / g.r5 * (2. + c5);

        for (int i = 0; i < 12; i++)
        {
            u[i] += du[i] * inflation_potency / (2.* PI);
        }
    }
    /* ===================================== */

    return; 
}

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
)
{
    // === Variable declaration ===
    double ua0[12] = {0}, ua0_image[12] = {0}, ub0[12] = {0}, uc0[12] = {0};
    geometric_constants_ps geometric_constants_ps;
    trigonometric_constants trigonometric_constants;
    alpha_constants alpha_constants;
    double effective_depth;

    trigonometric_constants = compute_trig_functions(dip);
    alpha_constants = compute_alpha_constants(alpha);

    // === Compute infinite-medium terms (real and image) ===
    // Somigliani tensor for real source
    effective_depth = depth + z;
    geometric_constants_ps = compute_geometric_constants_ps(
        x,
        y,
        z,
        effective_depth,
        trigonometric_constants
    );
    compute_ua0(
        geometric_constants_ps,
        trigonometric_constants,
        alpha_constants,
        strike_potency,
        dip_potency,
        tensile_potency,
        inflation_potency,
        ua0
    );

    // Somigliani tensor for the image source
    effective_depth = depth - z;
    geometric_constants_ps = compute_geometric_constants_ps(
        x,
        y,
        z,
        effective_depth,
        trigonometric_constants
    );
    compute_ua0(
        geometric_constants_ps,
        trigonometric_constants,
        alpha_constants,
        strike_potency,
        dip_potency,
        tensile_potency,
        inflation_potency,
        ua0_image
    );

    // === Compute surface deformation related term ===
    compute_ub0(
        geometric_constants_ps,
        trigonometric_constants,
        alpha_constants,
        strike_potency,
        dip_potency,
        tensile_potency,
        inflation_potency,
        ub0
    );

    // === Compute depth-multiplied term ===
    compute_uc0(
        geometric_constants_ps,
        trigonometric_constants,
        alpha_constants,
        strike_potency,
        dip_potency,
        tensile_potency,
        inflation_potency,
        uc0
    );

    // === Assemble displacement field ===
    for (int i = 0; i < 12; i++)
    {
        u[i] += ua0_image[i] - ua0[i] + ub0[i] + z*uc0[i];
        if (i >= 9)
        {
            u[i] += 2*ua0[i] + uc0[i-9];
        }
        u[i] /= (1e6); // * (1. - alpha));  // THIS FACTOR DIFFERS FROM COULOMB
    }

    return;
}
