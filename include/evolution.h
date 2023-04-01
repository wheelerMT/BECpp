//
// Created by mattw on 11/01/2022.
//

#ifndef BECPP_EVOLUTION_H
#define BECPP_EVOLUTION_H

#include <complex>
#include "wavefunction.h"
#include "data.h"

constexpr std::complex<double> I{0, 1};

void apply_TF_density(Wavefunction2D &psi, const Parameters &params)
{
    double tf_density;

    double g = params.c0 + 4 * params.c2;
    double r_tf = std::pow(8 * g / PI, 0.25);

    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        for (int j = 0; j < psi.grid.m_yPoints; ++j)
        {
            double r2 = psi.grid.m_xMesh[i][j] * psi.grid.m_xMesh[i][j] + psi.grid.m_yMesh[i][j] * psi.grid.m_yMesh[i][j];

            if (r2 < r_tf)
            {
                tf_density = 15 / (8 * PI * r_tf) * (1 - r2 / (r_tf * r_tf));

            } else
            {
                tf_density = 0.;
            }

            psi.plus[j + i * psi.grid.m_yPoints] *= tf_density;
            psi.zero[j + i * psi.grid.m_yPoints] *= tf_density;
            psi.minus[j + i * psi.grid.m_yPoints] *= tf_density;
        }
    }

    psi.update_component_atom_num();
}

void fourier_step(Wavefunction2D &psi, const Parameters &params)
{
#pragma omp parallel for collapse(2) shared(psi, params, I) default(none)
    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        for (int j = 0; j < psi.grid.m_yPoints; ++j)
        {
            psi.plus_k[j + i * psi.grid.m_yPoints] *= exp(-0.25 * I * params.dt * (psi.grid.m_wavenumber[i][j] + 2 * params.q));
            psi.zero_k[j + i * psi.grid.m_yPoints] *= exp(-0.25 * I * params.dt * psi.grid.m_wavenumber[i][j]);
            psi.minus_k[j + i * psi.grid.m_yPoints] *= exp(-0.25 * I * params.dt * (psi.grid.m_wavenumber[i][j] + 2 * params.q));
        }
    }
}

void fourier_step(Wavefunction1D &psi, const Parameters &params)
{
#pragma omp parallel for shared(psi, params, I) default(none)
    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        psi.plus_k[i] *= exp(-0.25 * I * params.dt * (psi.grid.m_wavenumber[i] + 2 * params.q));
        psi.zero_k[i] *= exp(-0.25 * I * params.dt * psi.grid.m_wavenumber[i]);
        psi.minus_k[i] *= exp(-0.25 * I * params.dt * (psi.grid.m_wavenumber[i] + 2 * params.q));

    }
}

void fourier_step_KZ(Wavefunction1D &psi, const Parameters &params, int tau_q)
{
#pragma omp parallel for shared(psi, params, I, tau_q) default(none)
    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        psi.plus_k[i] *= exp(-0.25 * I * params.dt *
                             (psi.grid.m_wavenumber[i] + 2 * abs(params.c2) * (params.q - params.dt.real() / (2 * tau_q))));
        psi.zero_k[i] *= exp(-0.25 * I * params.dt * psi.grid.m_wavenumber[i]);
        psi.minus_k[i] *= exp(-0.25 * I * params.dt *
                              (psi.grid.m_wavenumber[i] + 2 * abs(params.c2) * (params.q - params.dt.real() / (2 * tau_q))));
    }
}

void interaction_step(Wavefunction2D &psi, const Parameters &params, const doubleArray_t &V)
{
#pragma omp parallel for collapse(2) shared(psi, params, I, V) default(none)
    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        for (int j = 0; j < psi.grid.m_yPoints; ++j)
        {
            // Calculate spin vector elements
            std::complex<double> f_perp =
                    sqrt(2.) * (std::conj(psi.plus[j + i * psi.grid.m_yPoints]) * psi.zero[j + i * psi.grid.m_yPoints]
                                + std::conj(psi.zero[j + i * psi.grid.m_yPoints]) * psi.minus[j + i * psi.grid.m_yPoints]);
            std::complex<double> f_z = std::pow(abs(psi.plus[j + i * psi.grid.m_yPoints]), 2) -
                                       std::pow(abs(psi.minus[j + i * psi.grid.m_yPoints]), 2);
            double F = sqrt(std::pow(abs(f_z), 2) + std::pow(abs(f_perp), 2));

            // Calculate trigonometric expressions
            std::complex<double> C = std::cos(params.c2 * F * params.dt);
            std::complex<double> S{};
            if (F > 1e-8)
            {
                S = I * std::sin(params.c2 * F * params.dt) / F;
            }

            // Calculate density
            double n = std::pow(abs(psi.plus[j + i * psi.grid.m_yPoints]), 2) +
                       std::pow(abs(psi.zero[j + i * psi.grid.m_yPoints]), 2) +
                       std::pow(abs(psi.minus[j + i * psi.grid.m_yPoints]), 2);

            // Solve interaction part of flow
            std::complex<double> new_psi_plus = (C * psi.plus[j + i * psi.grid.m_yPoints] -
                                                 S * (f_z * psi.plus[j + i * psi.grid.m_yPoints]
                                                      + std::conj(f_perp) / sqrt(2.) * psi.zero[j + i * psi.grid.m_yPoints]))
                                                * exp(-I * params.dt * (V[i][j] - params.p + params.c0 * n));

            // Solve interaction part of flow
            std::complex<double> new_psi_zero = (C * psi.zero[j + i * psi.grid.m_yPoints] -
                                                 S / sqrt(2.) * (f_perp * psi.plus[j + i * psi.grid.m_yPoints]
                                                                 + std::conj(f_perp) * psi.minus[j + i * psi.grid.m_yPoints]))
                                                * exp(-I * params.dt * (V[i][j] + params.c0 * n));

            // Solve interaction part of flow
            std::complex<double> new_psi_minus = (C * psi.minus[j + i * psi.grid.m_yPoints] -
                                                  S * (f_perp / sqrt(2.) * psi.zero[j + i * psi.grid.m_yPoints]
                                                       - f_z * psi.minus[j + i * psi.grid.m_yPoints]))
                                                 * exp(-I * params.dt * (V[i][j] + params.p + params.c0 * n));

            // Update wavefunction
            psi.plus[j + i * psi.grid.m_yPoints] = new_psi_plus;
            psi.zero[j + i * psi.grid.m_yPoints] = new_psi_zero;
            psi.minus[j + i * psi.grid.m_yPoints] = new_psi_minus;

        }
    }
}

void interaction_step(Wavefunction1D &psi, const Parameters &params, const std::vector<double> &V)
{
#pragma omp parallel for shared(psi, params, I, V) default(none)
    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        // Calculate spin vector elements
        std::complex<double> f_perp =
                sqrt(2.) * (std::conj(psi.plus[i]) * psi.zero[i] + std::conj(psi.zero[i]) * psi.minus[i]);
        std::complex<double> f_z = std::pow(abs(psi.plus[i]), 2) - std::pow(abs(psi.minus[i]), 2);
        double F = sqrt(std::pow(abs(f_z), 2) + std::pow(abs(f_perp), 2));

        // Calculate trigonometric expressions
        std::complex<double> C = std::cos(params.c2 * F * params.dt);
        std::complex<double> S{};
        if (F > 1e-8)
        {
            S = I * std::sin(params.c2 * F * params.dt) / F;
        }

        // Calculate density
        double n = std::pow(abs(psi.plus[i]), 2) +
                   std::pow(abs(psi.zero[i]), 2) +
                   std::pow(abs(psi.minus[i]), 2);

        // Solve interaction part of flow
        std::complex<double> new_psi_plus = (C * psi.plus[i] -
                                             S * (f_z * psi.plus[i] + std::conj(f_perp) / sqrt(2.) * psi.zero[i]))
                                            * exp(-I * params.dt * (V[i] - params.p + params.c0 * n));

        // Solve interaction part of flow
        std::complex<double> new_psi_zero = (C * psi.zero[i] -
                                             S / sqrt(2.) * (f_perp * psi.plus[i]
                                                             + std::conj(f_perp) * psi.minus[i]))
                                            * exp(-I * params.dt * (V[i] + params.c0 * n));

        // Solve interaction part of flow
        std::complex<double> new_psi_minus = (C * psi.minus[i] - S * (f_perp / sqrt(2.) * psi.zero[i]
                                                                      - f_z * psi.minus[i])) *
                                             exp(-I * params.dt * (V[i] + params.p + params.c0 * n));

        // Update wavefunction
        psi.plus[i] = new_psi_plus;
        psi.zero[i] = new_psi_zero;
        psi.minus[i] = new_psi_minus;
    }
}

void renormalise_atom_num(Wavefunction2D &psi)
{
    double current_N_plus = psi.component_atom_number("plus");
    double current_N_zero = psi.component_atom_number("zero");
    double current_N_minus = psi.component_atom_number("minus");

    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        for (int j = 0; j < psi.grid.m_yPoints; ++j)
        {
            psi.plus[j + i * psi.grid.m_yPoints] *= sqrt(psi.N_plus) / sqrt(current_N_plus);
            psi.zero[j + i * psi.grid.m_yPoints] *= sqrt(psi.N_zero) / sqrt(current_N_zero);
            psi.minus[j + i * psi.grid.m_yPoints] *= sqrt(psi.N_minus) / sqrt(current_N_minus);
        }
    }
}

void renormalise_atom_num(Wavefunction1D &psi)
{
    double current_N_plus = psi.component_atom_number("plus");
    double current_N_zero = psi.component_atom_number("zero");
    double current_N_minus = psi.component_atom_number("minus");

    for (int i = 0; i < psi.grid.m_xPoints; ++i)
    {
        psi.plus[i] *= sqrt(psi.N_plus) / sqrt(current_N_plus);
        psi.zero[i] *= sqrt(psi.N_zero) / sqrt(current_N_zero);
        psi.minus[i] *= sqrt(psi.N_minus) / sqrt(current_N_minus);
    }
}

#endif //BECPP_EVOLUTION_H
