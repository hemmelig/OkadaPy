/*
 * =============================================================================
 *
 *       Filename:  okada.c
 *
 *        Purpose:  Implementation of the analytical deformation equations,
 *                  after Okada, 1992.
 *
 *      Copyright:  Conor A. Bacon, 2024
 *        License:  GNU General Public License, Version 3
 *                  (https://www.gnu.org/licenses/gpl-3.0.html)
 *
 * =============================================================================
 */

#include "libokada.h"

/*
 * Function: compute_okada_stress
 * ------------------------------
 * Compute the stress tensor using the equations set out in Okada, 1992.
 *
 * x_coords: x-component of coordinates at which to compute stress tensors.
 * y_coords: y-component of coordinates at which to compute stress tensors.
 * n_coords: The number of coordinates at which to compute stress tensors.
 * model_elements: Elements of the deformation model.
 * n_elements: Total number of elements in deformation model.
 * youngs_mod: The Young's modulus for the elastic half-space.
 * poisson: The Poisson's ratio for the elastic half-space.
 * calculation_depth: The depth at which to compute stress tensors.
 * resultant_stress_tensor: Pre-allocated array structure to write tensors.
 * threads: The number of computer threads to use for computation.
 *
 * returns: void - computed stress tensors are written to pre-allocated memory.
 */
void compute_okada_stress(
    double *x_coords,
    double *y_coords,
    int n_coords,
    double *model_elements,
    int n_elements,
    double youngs_mod,
    double poisson,
    double calculation_depth,
    double *resultant_stress_tensor,
    int threads
)
{
    int i;

    double alpha = 1. / (2. * (1. - poisson));
    calculation_depth *= -1.;

    #pragma omp parallel for num_threads(threads)
    for (i = 0; i < n_coords; i++)
    {
        // === Variable declarations ===
        double *model_element, *u, *du, *ui, *s, *si;
        double element_depth, sk, gk, vol;

        s = &resultant_stress_tensor[i * 6];
        for (int j = 0; j < n_elements; j++)
        {
            double ui[3] = {0};
            double si[6] = {0};
            model_element = &model_elements[j * 10];

            element_depth = (model_element[8] + model_element[9]) / 2.;

            // Convert from base coordinates to Okada fault coordinates
            double converted_coords[4] = {0};
            convert2okada_coordinates(
                x_coords[i],
                y_coords[i],
                model_element[0],
                model_element[1],
                model_element[2],
                model_element[3],
                model_element[8],
                model_element[9],
                model_element[7],
                converted_coords
            );

            double du[12] = {0};
            switch ((int) model_element[4])
            {
                case 100:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        -model_element[5],   // Strike-slip dislocation
                        model_element[6],    // Dip-slip dislocation
                        0.,                  // Tensile-fault dislocation
                        du
                    );
                    break;

                case 200:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        -model_element[5],   // Strike-slip dislocation
                        0.,                  // Dip-slip dislocation
                        model_element[6],    // Tensile-fault dislocation
                        du
                    );
                    break;

                case 300:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        0.,                  // Strike-slip dislocation
                        model_element[6],    // Dip-slip dislocation
                        model_element[5],    // Tensile-fault dislocation
                        du
                    );
                    break;

                case 400:
                    compute_ps_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        -model_element[5],   // Strike-slip potency
                        model_element[6],    // Dip-slip potency
                        0.,                  // Tensile-fault potency
                        0.,                  // Inflation potency
                        du
                    );
                    break;

                case 500:
                    compute_ps_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        0.,                  // Strike-slip potency
                        0.,                  // Dip-slip potency
                        model_element[5],    // Tensile-fault potency
                        model_element[6],    // Inflation potency
                        du
                    );
                    break;
            }

            // Convert displacement from Okada's field to given field
            double sw = sqrt(
                (model_element[3] - model_element[1])
                    * (model_element[3] - model_element[1])
                + (model_element[2] - model_element[0])
                    * (model_element[2] - model_element[0])
            );
            double sina = (model_element[3] - model_element[1]) / sw;
            double cosa = (model_element[2] - model_element[0]) / sw;

            // Convert displacement to strain
            sk = youngs_mod / (1. + poisson);
            gk = poisson / (1. - 2. * poisson);
            vol = du[3] + du[7] + du[11];

            si[0] = sk * (gk * vol + du[3]) * 0.001;  // Sxx
            si[1] = sk * (gk * vol + du[7]) * 0.001;  // Syy
            si[2] = sk * (gk * vol + du[11]) * 0.001;  // Szz
            si[5] = sk * (du[4] + du[6]) / 2. * 0.001;  // Syz
            si[4] = sk * (du[5] + du[9]) / 2. * 0.001;  // Sxz
            si[3] = sk * (du[8] + du[10]) / 2. * 0.001;  // Sxy

            double si_transformed[6] = {0};
            transform_tensor(sina, cosa, si, si_transformed);

            // Linearly stack with current stress tensor at this position
            for (int k = 0; k < 6; k++)
            {
                s[k] += si_transformed[k];
            }
        }
    }
}

/*
 * Function: compute_okada_strain
 * ------------------------------
 * Compute the strain tensor using the equations set out in Okada, 1992.
 *
 * x_coords: x-component of coordinates at which to compute strain tensors.
 * y_coords: y-component of coordinates at which to compute strain tensors.
 * n_coords: The number of coordinates at which to compute strain tensors.
 * model_elements: Elements of the deformation model.
 * n_elements: Total number of elements in deformation model.
 * youngs_mod: The Young's modulus for the elastic half-space.
 * poisson: The Poisson's ratio for the elastic half-space.
 * calculation_depth: The depth at which to compute strain tensors.
 * resultant_strain_tensor: Pre-allocated array structure to write tensors.
 * threads: The number of computer threads to use for computation.
 *
 * returns: void - computed strain tensors are written to pre-allocated memory.
 */
void compute_okada_strain(
    double *x_coords,
    double *y_coords,
    int n_coords,
    double *model_elements,
    int n_elements,
    double youngs_mod,
    double poisson,
    double calculation_depth,
    double *resultant_strain_tensor,
    int threads
)
{
    int i;

    double alpha = 1. / (2. * (1. - poisson));
    calculation_depth *= -1.;

    #pragma omp parallel for num_threads(threads)
    for (i = 0; i < n_coords; i++)
    {
        // === Variable declarations ===
        double *model_element, *u, *du, *ui, *e, *ei;//, *ei_transformed;
        double element_depth;

        e = &resultant_strain_tensor[i * 6];
        for (int j = 0; j < n_elements; j++)
        {
            double ui[3] = {0};
            double ei[6] = {0};
            model_element = &model_elements[j * 10];

            element_depth = (model_element[8] + model_element[9]) / 2.;

            // Convert from base coordinates to Okada fault coordinates
            double converted_coords[4] = {0};
            convert2okada_coordinates(
                x_coords[i],
                y_coords[i],
                model_element[0],
                model_element[1],
                model_element[2],
                model_element[3],
                model_element[8],
                model_element[9],
                model_element[7],
                converted_coords
            );

            double du[12] = {0};
            switch ((int) model_element[4])
            {
                case 100:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        -model_element[5],   // Strike-slip dislocation
                        model_element[6],    // Dip-slip dislocation
                        0.,                  // Tensile-fault dislocation
                        du
                    );
                    break;

                case 200:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        -model_element[5],   // Strike-slip dislocation
                        0.,                  // Dip-slip dislocation
                        model_element[6],    // Tensile-fault dislocation
                        du
                    );
                    break;

                case 300:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        0.,                  // Strike-slip dislocation
                        model_element[6],    // Dip-slip dislocation
                        model_element[5],    // Tensile-fault dislocation
                        du
                    );
                    break;

                case 400:
                    compute_ps_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        -model_element[5],   // Strike-slip potency
                        model_element[6],    // Dip-slip potency
                        0.,                  // Tensile-fault potency
                        0.,                  // Inflation potency
                        du
                    );
                    break;

                case 500:
                    compute_ps_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        0.,                  // Strike-slip potency
                        0.,                  // Dip-slip potency
                        model_element[5],    // Tensile-fault potency
                        model_element[6],    // Inflation potency
                        du
                    );
                    break;
            }

            // Convert displacement from Okada's field to given field
            double sw = sqrt(
                (model_element[3] - model_element[1])
                    * (model_element[3] - model_element[1])
                + (model_element[2] - model_element[0])
                    * (model_element[2] - model_element[0])
            );
            double sina = (model_element[3] - model_element[1]) / sw;
            double cosa = (model_element[2] - model_element[0]) / sw;

            // Convert displacement to strain
            ei[0] = du[3] * 0.001;
            ei[1] = du[7] * 0.001;
            ei[2] = du[11] * 0.001;
            ei[5] = (du[4] + du[6]) / 2. * 0.001;
            ei[4] = (du[5] + du[9]) / 2. * 0.001;
            ei[3] = (du[8] + du[10]) / 2. * 0.001;

            double ei_transformed[6] = {0};
            transform_tensor(sina, cosa, ei, ei_transformed);

            // Linearly stack with current strain tensor at this position
            for (int k = 0; k < 6; k++)
            {
                e[k] += ei_transformed[k];
            }
        }
    }
}

/*
 * Function: compute_okada_displacement
 * ------------------------------------
 * Compute the displacement using the equations set out in Okada, 1992.
 *
 * x_coords: x-component of coordinates at which to compute displacements.
 * y_coords: y-component of coordinates at which to compute displacements.
 * n_coords: The number of coordinates at which to compute displacements.
 * model_elements: Elements of the deformation model.
 * n_elements: Total number of elements in deformation model.
 * youngs_mod: The Young's modulus for the elastic half-space.
 * poisson: The Poisson's ratio for the elastic half-space.
 * calculation_depth: The depth at which to compute displacements.
 * resultant_displacement: Pre-allocated array structure to write tensors.
 * threads: The number of computer threads to use for computation.
 *
 * returns: void - computed displacements are written to pre-allocated memory.
 */
void compute_okada_displacement(
    double *x_coords,
    double *y_coords,
    int n_coords,
    double *model_elements,
    int n_elements,
    double youngs_mod,
    double poisson,
    double calculation_depth,
    double *resultant_displacement,
    int threads
)
{
    int i;

    double alpha = 1. / (2. * (1. - poisson));
    calculation_depth *= -1.;

    #pragma omp parallel for num_threads(threads)
    for (i = 0; i < n_coords; i++)
    {
        // === Variable declarations ===
        double *model_element, *u, *du, *ui, *s, *si, *sit;
        double element_depth;

        u = &resultant_displacement[i * 12];
        for (int j = 0; j < n_elements; j++)
        {
            double ui[12] = {0};
            model_element = &model_elements[j * 10];

            element_depth = (model_element[8] + model_element[9]) / 2.;

            // Convert from base coordinates to Okada fault coordinates
            double converted_coords[4] = {0};
            convert2okada_coordinates(
                x_coords[i],
                y_coords[i],
                model_element[0],
                model_element[1],
                model_element[2],
                model_element[3],
                model_element[8],
                model_element[9],
                model_element[7],
                converted_coords
            );

            double du[12] = {0};
            switch ((int) model_element[4])
            {
                case 100:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        -model_element[5],   // Strike-slip dislocation
                        model_element[6],    // Dip-slip dislocation
                        0.,                  // Tensile-fault dislocation
                        du
                    );
                    break;

                case 200:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        -model_element[5],   // Strike-slip dislocation
                        0.,                  // Dip-slip dislocation
                        model_element[6],    // Tensile-fault dislocation
                        du
                    );
                    break;

                case 300:
                    compute_frs_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        converted_coords[2], // Fault half length
                        converted_coords[3], // Fault half width
                        0.,                  // Strike-slip dislocation
                        model_element[6],    // Dip-slip dislocation
                        model_element[5],    // Tensile-fault dislocation
                        du
                    );
                    break;

                case 400:
                    compute_ps_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        -model_element[5],   // Strike-slip potency
                        model_element[6],    // Dip-slip potency
                        0.,                  // Tensile-fault potency
                        0.,                  // Inflation potency
                        du
                    );
                    break;

                case 500:
                    compute_ps_deformation(
                        alpha,
                        converted_coords[0], // x
                        converted_coords[1], // y
                        calculation_depth,
                        element_depth,
                        model_element[7],    // Dip
                        0.,                  // Strike-slip potency
                        0.,                  // Dip-slip potency
                        model_element[5],    // Tensile-fault potency
                        model_element[6],    // Inflation potency
                        du
                    );
                    break;
            }

            // Convert displacement from Okada's field to given field
            double sw = sqrt(
                (model_element[3] - model_element[1])
                    * (model_element[3] - model_element[1])
                + (model_element[2] - model_element[0])
                    * (model_element[2] - model_element[0])
            );
            double sina = (model_element[3] - model_element[1]) / sw;
            double cosa = (model_element[2] - model_element[0]) / sw;
            ui[0] = du[0] * cosa - du[1] * sina;
            ui[1] = du[0] * sina + du[1] * cosa;
            ui[2] = du[2];

            // Linearly stack with current displacement vector at this position
            for (int k = 0; k < 3; k++)
            {
                u[k] += ui[k];
            }
            for (int k = 3; k < 12; k++)
            {
                u[k] += du[k];
            }
        }
    }
}
