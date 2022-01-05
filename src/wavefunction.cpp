//
// Created by mattw on 05/01/2022.
//

#include "wavefunction.h"

void Wavefunction::generateInitialState(const std::string &gs_phase)
{
    if (gs_phase == "polar")
    {
        std::cout << "Generating polar initial state...\n";
        std::fill(plus.begin(), plus.end(), 0.);
        std::fill(zero.begin(), zero.end(), 1.);
        std::fill(minus.begin(), minus.end(), 0.);
    } else if (gs_phase == "FM")
    {
        std::cout << "Generating ferromagnetic initial state...\n";
        std::fill(plus.begin(), plus.end(), 1 / sqrt(2.));
        std::fill(zero.begin(), zero.end(), 0.);
        std::fill(minus.begin(), minus.end(), 1 / sqrt(2.));
    }
}

Wavefunction::Wavefunction(Grid &t_grid, const std::string &gs_phase) : grid{t_grid}
{
    // Size arrays
    plus.resize(grid.nx * grid.ny);
    zero.resize(grid.nx * grid.ny);
    minus.resize(grid.nx * grid.ny);
    plus_k.resize(grid.nx * grid.ny);
    zero_k.resize(grid.nx * grid.ny);
    minus_k.resize(grid.nx * grid.ny);

    // Populates wavefunction components
    generateInitialState(gs_phase);
}

void Wavefunction::add_noise(const std::string &components, double mean, double stddev)
{
    // Construct random generator
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator{seed1};
    std::normal_distribution<double> norm_dist{mean, stddev};

    if (components == "outer")
    {
        std::cout << "Adding noise...\n";

        // Add noise to outer components
        for (int i = 0; i < grid.nx; i++)
        {
            for (int j = 0; j < grid.ny; j++)
            {
                plus[j + i * grid.nx] +=
                        std::complex<double>{norm_dist(generator), norm_dist(generator)};
                minus[j + i * grid.nx] +=
                        std::complex<double>{norm_dist(generator), norm_dist(generator)};
            }
        }
    }
}

doubleArray_t Wavefunction::density()
{
    doubleArray_t density{};
    density.resize(grid.nx, std::vector<double>(grid.ny));

    for (int i = 0; i < grid.nx; i++)
    {
        for (int j = 0; j < grid.ny; j++)
        {
            density[i][j] += (std::pow(std::abs(plus[j + grid.nx * i]), 2)
                              + std::pow(std::abs(zero[j + grid.nx * i]), 2)
                              + std::pow(std::abs(minus[j + grid.nx * i]), 2));
        }
    }
    return density;
}
