//
// Created by mattw on 11/01/2022.
//

#ifndef BECPP_EVOLUTION_H
#define BECPP_EVOLUTION_H

#include <complex>
#include "wavefunction.h"
#include "data.h"

constexpr std::complex<double> I{0, 1};

void fourier_step(Wavefunction &psi, const Parameters &params)
{
    for (int i = 0; i < psi.grid.nx; ++i)
    {
        for (int j = 0; j < psi.grid.ny; ++j)
        {
            psi.plus[j + i * psi.grid.nx] *= exp(-0.25 * I * params.dt * psi.grid.K[i][j] + 2 * params.q);
            psi.zero[j + i * psi.grid.nx] *= exp(-0.25 * I * params.dt * psi.grid.K[i][j]);
            psi.minus[j + i * psi.grid.nx] *= exp(-0.25 * I * params.dt * psi.grid.K[i][j] + 2 * params.q);
        }
    }
}


void interaction_step(Wavefunction &psi, const Parameters &params)
{

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
            psi.plus[j + i * psi.grid.nx] = new_psi_zero;
            psi.plus[j + i * psi.grid.nx] = new_psi_minus;

        }
    }
}

#endif //BECPP_EVOLUTION_H
