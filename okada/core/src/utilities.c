/*
 * =============================================================================
 *
 *       Filename:  utilities.c
 *
 *        Purpose:  Collection of utility equations for computing constants.
 *
 *      Copyright:  Conor A. Bacon, 2024
 *        License:  GNU General Public License, Version 3
 *                  (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * =============================================================================
 */

 #include "libokada.h"


/*
 * Function: compute_trig_functions
 * --------------------------------
 * Pre-compute a collection of trigonemetric constants.
 * 
 * delta: the angle of dip of the plane spanned by the deformation element.
 *
 * returns: trigonometric_constants - struct containing computed constants.
 */
trigonometric_constants compute_trig_functions(double delta)
{
	trigonometric_constants result;

    result.s = sin(delta*PI/180.);
    result.c = cos(delta*PI/180.);

    // Handle tiny cos(delta) elements
    // for (int i; i < len(delta); i++) {
    //     if cos_delta[i] < eps {
    //         cos_delta[i] = 0;
    //     }
    // }

    // Handle some shit with sin(delta) too (FINISH)
    // for (int i; i < len(delta); i++) {
    //     if sin_delta[i] 
    // }

    result.ss = result.s * result.s;
    result.cc = result.c * result.c;
    result.sc = result.s * result.c;
    result.s2 = 2 * result.s * result.c;
    result.c2 = result.cc - result.ss;

    return result;
}

/*
 * Function: compute_alpha_constants
 * ---------------------------------
 * Pre-compute a collection of constants containing the alpha term, which is
 * a function of the Lamé constants.
 * 
 * alpha: an elastic constant defined by the Lamé constants.
 *
 * returns: alpha_constants - struct containing computed constants.
 */
alpha_constants compute_alpha_constants(double alpha)
{
    alpha_constants result =
    {
        .a1 = (1. - alpha) / 2.,
        .a2 = alpha / 2.,
        .a3 = (1. - alpha) / alpha,
        .a4 = 1. - alpha,
        .a5 = alpha
    };

    return result;
}

/*
 * Function: compute_geometric_constants_ps
 * ----------------------------------------
 * Pre-compute a collection of geometric constants for the infinitesimal point
 * source equations.
 *
 * x: x-coordinate of point in space for which to compute displacement.
 * y: y-coordinate of point in space for which to compute displacement.
 * z: z-coordinate of point in space for which to compute displacement.
 * d: Effective depth of point in space for which to compute displacement.
 * t: Struct containing some pre-computed trigonometric constants.
 *
 * returns: geometric_constants_ps - struct containing geometric constants.
 */
geometric_constants_ps compute_geometric_constants_ps(
    double x, // Units: km
    double y, // Units: km
    double z, // Units: km
    double d, // Units: km
    trigonometric_constants t
)
{
    geometric_constants_ps result;

    result.x = x; // Units: km
    result.y = y; // Units: km
    result.d = d; // Units: km

    result.p = y * t.c + d * t.s; // Units: km
    result.q = y * t.s - d * t.c; // Units: km
    result.s = result.p * t.s + result.q * t.c; // Units: km
    result.t = result.p * t.c - result.q * t.s; // Units: km

    result.xy = x * y; // Units: km**2
    result.x2 = x * x; // Units: km**2
    result.y2 = y * y; // Units: km**2
    result.d2 = d * d; // Units: km**2

    result.r2 = result.x2 + result.y2 + result.d2; // Units: km**2
    result.r = sqrt(result.r2); // Units: km
    result.r3 = pow(result.r, 3); // Units: km**3
    result.r5 = pow(result.r, 5); // Units: km**5
    result.r7 = pow(result.r, 7); // Units: km**7

    result.a3 = 1 - 3 * result.x2 / result.r2; // Dimensionless
    result.a5 = 1 - 5 * result.x2 / result.r2; // Dimensionless
    result.b3 = 1 - 3 * result.y2 / result.r2; // Dimensionless
    result.c3 = 1 - 3 * result.d2 / result.r2; // Dimensionless

    result.qr = 3 * result.q / result.r5; // Units: km**-4
    result.qrx = 5 * result.qr * x / result.r2; // Units: km**-5

    result.uy = t.s - 5 * y * result.q / result.r2; // Dimensionless
    result.uz = t.c + 5 * d * result.q / result.r2; // Dimensionless
    result.vy = result.s - 5 * y * result.p * result.q / result.r2; // Units: km
    result.vz = result.t + 5 * d * result.p * result.q / result.r2; // Units: km
    result.wy = result.uy + t.s; // Dimensionless
    result.wz = result.uz + t.c; // Dimensionless

    return result;
}

/*
 * Function: compute_geometric_constants_frs
 * -----------------------------------------
 * Pre-compute a collection of geometric constants for the finite rectangular
 * source equations.
 * 
 * xi: Transformed coordinate in space for which to compute displacement.
 * eta: Transformed coordinate in space for which to compute displacement.
 * q: Transformed coordinate in space for which to compute displacement.
 * t: Struct containing some pre-computed trigonometric constants.
 *
 * returns: geometric_constants_frs - struct containing geometric constants.
 */
geometric_constants_frs compute_geometric_constants_frs(
    double xi,
    double eta,
    double q,
    trigonometric_constants t
)
{
    double eps = 0.000001;

    if (fabs(xi) < eps)
    {
        xi = 0.;
    }

    if (fabs(eta) < eps)
    {
        eta = 0.;
    }

    if (fabs(q) < eps)
    {
        q = 0.;
    }

    geometric_constants_frs result;
    
    result.xi = xi;
    result.xi2 = xi * xi;
    result.eta = eta;
    result.eta2 = eta * eta;
    result.q = q;
    result.q2 = q * q;
    result.r2 = result.xi2 + result.eta2 + result.q2;
    result.r = sqrt(result.r2);

    // Handle r = 0
    // if (result.r) == 0
    // {
    //     return;
    // }

    result.r3 = result.r2 * result.r;
    result.r5 = result.r3 * result.r2;
    result.y = eta * t.c + q * t.s;
    result.d = eta * t.s - q * t.c;

    result.tt = atan(xi * eta / (result.q * result.r));

    if (xi < 0 && result.q == 0 && eta == 0)
    {
        result.alx = -log(result.r - xi);
        result.x11 = 0;
        result.x32 = 0;
    }
    else
    {
        result.rxi = result.r + xi;
        result.alx = log(result.rxi);
        result.x11 = 1. / (result.r * result.rxi);
        result.x32 = (result.r + result.rxi) * result.x11 * result.x11 / result.r;
    }

    if (eta < 0 && result.q == 0 && xi == 0)
    {
        result.ale = -log(result.r - eta);
        result.y11 = 0;
        result.y32 = 0;
    }
    else
    {
        double reta = result.r + eta;
        result.ale = log(reta);
        result.y11 = 1. / (result.r * reta);
        result.y32 = (result.r + reta) * result.y11 * result.y11 / result.r;
    }

    result.ey = t.s / result.r - result.y * result.q / result.r3;
    result.ez = t.c / result.r + result.d * result.q / result.r3;
    result.fy = result.d / result.r3 + result.xi2 * result.y32 * t.s;
    result.fz = result.y / result.r3 + result.xi2 * result.y32 * t.c;
    result.gy = 2. * result.x11 * t.s - result.y * result.q * result.x32;
    result.gz = 2. * result.x11 * t.c + result.d * result.q * result.x32;
    result.hy = result.d * result.q * result.x32 + xi * result.q * result.y32 * t.s;
    result.hz = result.y * result.q * result.x32 + xi * result.q * result.y32 * t.c;

    return result;
}

/*
 * Function: convert2okada_coordinates
 * -----------------------------------
 * Convert from real-space Cartesian frame to Okada-defined coordinate frame.
 * 
 * x_coord: x-coordinate of midpoint of deformation element.
 * y_coord: y-coordinate of midpoint of deformation element.
 * x_min: minimum x-coordinate defining the extent of deformation element.
 * y_min: minimum y-coordinate defining the extent of deformation element.
 * x_max: maximum x-coordinate defining the extent of deformation element.
 * y_max: maximum y-coordinate defining the extent of deformation element.
 * z_max: maximum z-coordinate defining the extent of deformation element.
 * z_min: minimum z-coorindate defining the extent of deformation element.
 * dip: dip of the plane spanned by the deformation element.
 * result: pointer to an array to which transformed coordinates are written.
 *
 * returns: void
 */
void convert2okada_coordinates(
    double x_coord,
    double y_coord,
    double x_min,
    double y_min,
    double x_max,
    double y_max,
    double z_max,
    double z_min,
    double dip,
    double *result
)
{
    double x_midpoint = (x_min + x_max) / 2.;
    double y_midpoint = (y_min + y_max) / 2.;
    double z_midpoint = (z_min - z_max) / 2.;

    double k = tan(dip * PI / 180.);
    if (fabs(k) <= 1e-6)
    {
        k = 1e-6;
    }

    double d = z_midpoint / k;
    double b = atan((y_max - y_min) / (x_max - x_min));

    double ydipshift = fabs(d * cos(b));
    double xdipshift = fabs(d * sin(b));

    if (x_max > x_min)
    {
        if (y_max > y_min)
        {
            x_midpoint += xdipshift;
            y_midpoint -= ydipshift;
        }
        else
        {
            x_midpoint -= xdipshift;
            y_midpoint -= ydipshift;
        }
    }
    else
    {
        if (y_max > y_min)
        {
            x_midpoint += xdipshift;
            y_midpoint += ydipshift;
        }
        else
        {
            x_midpoint -= xdipshift;
            y_midpoint += ydipshift;
        }
    }
    
    // Converting from global coordinate to Okada-fault coordinate
    result[0] = (x_coord - x_midpoint) * cos(b)
        + (y_coord - y_midpoint) * sin(b);
    result[1] = -(x_coord - x_midpoint) * sin(b)
        + (y_coord - y_midpoint) * cos(b);

    if (x_max - x_min < 0.)
    {
        result[0] *= -1.;
        result[1] *= -1.;
    }

    result[2] = sqrt((x_max - x_min) * (x_max - x_min)
        + (y_max - y_min) * (y_max - y_min)) / 2.;
    result[3] = ((z_min - z_max) / 2.) / sin(dip * PI / 180.);

    return;
}

/*
 * Function: transform_tensor
 * --------------------------
 * Transform a tensor from Okada's frame to the real-space frame.
 *
 * sin_beta: sine of angle between frames.
 * cos_beta: cosine of angle between frames.
 * tensor_in: Tensor in Okada's coordinate frame.
 * tensor_out: Tensor in real-space coordinate frame.
 *
 * returns: void
 */
void transform_tensor(
    double sin_beta,
    double cos_beta,
    double *tensor_in,
    double *tensor_out
)
{
    double transform[36] = {0};
    double beta = asin(sin_beta);
    double x_beta, x_del, y_beta, y_del, z_beta, z_del;

    if (cos_beta > 0.)
    {
        x_beta = -beta;
        x_del = 0.;
        y_beta = -beta + PI / 2.;
        y_del = 0.;
        z_beta = -beta - PI / 2.;
        z_del = PI / 2.;
    }
    else
    {
        x_beta = beta - PI;
        x_del = 0.;
        y_beta = beta - PI / 2.;
        y_del = 0.;
        z_beta = beta - PI / 2.;
        z_del = PI / 2.;
    }

    double xl = cos(x_del) * cos(x_beta);
    double xm = cos(x_del) * sin(x_beta);
    double xn = sin(x_del);

    double yl = cos(y_del) * cos(y_beta);
    double ym = cos(y_del) * sin(y_beta);
    double yn = sin(y_del);

    double zl = cos(z_del) * cos(z_beta);
    double zm = cos(z_del) * sin(z_beta);
    double zn = sin(z_del);

    transform[0] = xl * xl;
    transform[1] = xm * xm;
    transform[2] = xn * xn;
    transform[3] = 2. * xm * xn;
    transform[4] = 2. * xn * xl;
    transform[5] = 2. * xl * xm;

    transform[6]  = yl * yl;
    transform[7]  = ym * ym;
    transform[8]  = yn * yn;
    transform[9]  = 2. * ym * yn;
    transform[10] = 2. * yn * yl;
    transform[11] = 2. * yl * ym;

    transform[12] = zl * zl;
    transform[13] = zm * zm;
    transform[14] = zn * zn;
    transform[15] = 2. * zm * zn;
    transform[16] = 2. * zn * zl;
    transform[17] = 2. * zl * zm;

    transform[18] = yl * zl;
    transform[19] = ym * zm;
    transform[20] = yn * zn;
    transform[21] = ym * zn + zm * yn;
    transform[22] = yn * zl + zn * yl;
    transform[23] = yl * zm + zl * ym;

    transform[24] = zl * xl;
    transform[25] = zm * xn;
    transform[26] = zn * xn;
    transform[27] = xm * zn + zm * xn;
    transform[28] = xn * zl + zn * xl;
    transform[29] = xl * zm + zl * xm;

    transform[30] = xl * yl;
    transform[31] = xm * ym;
    transform[32] = xn * yn;
    transform[33] = xm * yn + ym * xn;
    transform[34] = xn * yl + yn * xl;
    transform[35] = xl * ym + yl * xm;

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            tensor_out[i] += transform[i * 6 + j] * tensor_in[j];
        }
    }
}
