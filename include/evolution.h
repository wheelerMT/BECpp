//
// Created by mattw on 11/01/2022.
//

#ifndef BECPP_EVOLUTION_H
#define BECPP_EVOLUTION_H

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

#endif //BECPP_EVOLUTION_H
