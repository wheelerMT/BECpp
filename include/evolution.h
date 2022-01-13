//
// Created by mattw on 11/01/2022.
//

#ifndef BECPP_EVOLUTION_H
#define BECPP_EVOLUTION_H

#include <complex>
#include "wavefunction.h"
#include "data.h"

constexpr std::complex<double> I{1, 0};

void apply_TF_density(Wavefunction &psi, const Parameters &params)
{
    double tf_density;

    double g = params.c0 + 4 * params.c2;
    double r_tf = std::pow(8 * g / PI, 0.25);

    for (int i = 0; i < psi.grid.nx; ++i)
    {
        for (int j = 0; j < psi.grid.ny; ++j)
        {
            double r2 = psi.grid.X[i][j] * psi.grid.X[i][j] + psi.grid.Y[i][j] * psi.grid.Y[i][j];

            if (r2 < r_tf)
            {
                tf_density = 15 / (8 * PI * r_tf) * (1 - r2 / (r_tf * r_tf));

            } else
            {
                tf_density = 0.;
            }

            psi.plus[j + i * psi.grid.nx] *= tf_density;
            psi.zero[j + i * psi.grid.nx] *= tf_density;
            psi.minus[j + i * psi.grid.nx] *= tf_density;
        }
    }

    psi.update_component_atom_num();
}

void fourier_step(Wavefunction &psi, const Parameters &params)
{
#pragma omp parallel for collapse(2) shared(psi, params) default(none)
    for (int i = 0; i < psi.grid.nx; ++i)
    {
        for (int j = 0; j < psi.grid.ny; ++j)
        {
            psi.plus_k[j + i * psi.grid.nx] *= exp(-0.25 * I * params.dt * (psi.grid.K[i][j] + 2 * params.q));
            psi.zero_k[j + i * psi.grid.nx] *= exp(-0.25 * I * params.dt * psi.grid.K[i][j]);
            psi.minus_k[j + i * psi.grid.nx] *= exp(-0.25 * I * params.dt * (psi.grid.K[i][j] + 2 * params.q));
        }
    }
}

void interaction_step(Wavefunction &psi, const Parameters &params)
{
#pragma omp parallel for collapse(2) shared(psi, params) default(none)
    for (int i = 0; i < psi.grid.nx; ++i)
    {
        for (int j = 0; j < psi.grid.ny; ++j)
        {
            // Calculate spin vector elements
            std::complex<double> f_perp =
                    sqrt(2.) * (std::conj(psi.plus[j + i * psi.grid.nx]) * psi.zero[j + i * psi.grid.nx]
                                + std::conj(psi.zero[j + i * psi.grid.nx]) * psi.minus[j + i * psi.grid.nx]);
            std::complex<double> f_z = std::pow(abs(psi.plus[j + i * psi.grid.nx]), 2) -
                                       std::pow(abs(psi.minus[j + i * psi.grid.nx]), 2);
            double F = sqrt(std::pow(abs(f_z), 2) + std::pow(abs(f_perp), 2));

            // Calculate trigonometric expressions
            double C = std::cos(params.c2 * F * params.dt);
            std::complex<double> S{};
            if (F > 1e-8)
            {
                S = I * std::sin(params.c2 * F * params.dt);
            }

            // Calculate density
            double n = std::pow(abs(psi.plus[j + i * psi.grid.nx]), 2) +
                       std::pow(abs(psi.zero[j + i * psi.grid.nx]), 2) +
                       std::pow(abs(psi.minus[j + i * psi.grid.nx]), 2);

            // Solve interaction part of flow
            std::complex<double> new_psi_plus = (C * psi.plus[j + i * psi.grid.nx] -
                                                 S * (f_z * psi.plus[j + i * psi.grid.nx]
                                                      + std::conj(f_perp) / sqrt(2.) * psi.zero[j + i * psi.grid.nx]))
                                                * exp(-I * params.dt * (params.V[i][j] - params.p + params.c0 * n));

            // Solve interaction part of flow
            std::complex<double> new_psi_zero = (C * psi.zero[j + i * psi.grid.nx] -
                                                 S / sqrt(2.) * (f_perp * psi.plus[j + i * psi.grid.nx]
                                                                 + std::conj(f_perp) * psi.minus[j + i * psi.grid.nx]))
                                                * exp(-I * params.dt * (params.V[i][j] + params.c0 * n));

            // Solve interaction part of flow
            std::complex<double> new_psi_minus = (C * psi.minus[j + i * psi.grid.nx] -
                                                  S * (f_perp / sqrt(2.) * psi.zero[j + i * psi.grid.nx]
                                                       - f_z * psi.minus[j + i * psi.grid.nx]))
                                                 * exp(-I * params.dt * (params.V[i][j] + params.p + params.c0 * n));

            // Update wavefunction
            psi.plus[j + i * psi.grid.nx] = new_psi_plus;
            psi.zero[j + i * psi.grid.nx] = new_psi_zero;
            psi.minus[j + i * psi.grid.nx] = new_psi_minus;

        }
    }
}

void renormalise_atom_num(Wavefunction &psi)
{
    double current_N_plus = psi.component_atom_number("plus");
    double current_N_zero = psi.component_atom_number("zero");
    double current_N_minus = psi.component_atom_number("minus");

    for (int i = 0; i < psi.grid.nx; ++i)
    {
        for (int j = 0; j < psi.grid.ny; ++j)
        {
            psi.plus[j + i * psi.grid.nx] *= sqrt(psi.N_plus) / sqrt(current_N_plus);
            psi.zero[j + i * psi.grid.nx] *= sqrt(psi.N_zero) / sqrt(current_N_zero);
            psi.minus[j + i * psi.grid.nx] *= sqrt(psi.N_minus) / sqrt(current_N_minus);
        }
    }
}

#endif //BECPP_EVOLUTION_H
