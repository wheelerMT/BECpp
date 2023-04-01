//
// Created by mattw on 13/01/2022.
//

#ifndef BECPP_PHASE_H
#define BECPP_PHASE_H

#include <iostream>
#include <tuple>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include "grid.h"
#include "constants.h"

inline double heaviside(double x)
{
    return (x > 0) ? 1 : ((x < 0) ? 0 : 1);
}

std::vector<std::tuple<double, double>> generate_positions(const int n_vort, const double threshold, const Grid2D &grid,
                                                           const int max_iter)
{
    std::cout << "Finding " << n_vort << " vortex positions...\n";

    std::vector<std::tuple<double, double>> positions;
    int iterations = 0;

    // Construct random generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> uniform_distribution(-grid.m_xLength / 2, grid.m_xLength / 2);

    while (positions.size() < n_vort)
    {
        std::tuple<double, double> pos;
        pos = std::make_tuple(uniform_distribution(generator), uniform_distribution(generator));
        iterations += 1;

        bool triggered{false};

        for (const auto &accepted_pos: positions)
        {
            if (std::abs(std::get<0>(pos) - std::get<0>(accepted_pos)) < threshold)
            {
                if (std::abs(std::get<1>(pos) - std::get<1>(accepted_pos)) < threshold)
                {
                    triggered = true;
                    break;
                }
            }
        }

        if (!triggered)
        {
            positions.push_back(pos);
        }

        // If iterations exceed the maximum, return the current position list and continue
        if (iterations > max_iter)
        {
            std::cout << "WARNING: Max iterations exceeded, only found " << positions.size() << " suitable positions\n";
            return positions;
        }
    }

    std::cout << "Found " << n_vort << " positions in " << iterations << " iterations\n";

    return positions;
}

doubleArray_t construct_phase(const int n_vort, const double threshold, const Grid2D &grid, const int max_iter = 10000)
{
    std::cout << "Commencing construction of phase:\n";

    std::vector<std::tuple<double, double>> positions = generate_positions(n_vort, threshold, grid, max_iter);

    std::cout << "Constructing phase profile array...\n";

    // Phase array
    doubleArray_t theta;
    theta.resize(grid.m_xPoints, std::vector<double>(grid.m_yPoints));

    for (int num = 0; num < n_vort / 2; ++num)
    {
        doubleArray_t theta_k;
        theta_k.resize(grid.m_xPoints, std::vector<double>(grid.m_yPoints));

        // Extract positions
        auto[x_m, y_m] = positions[num];
        auto[x_p, y_p] = positions[n_vort / 2 + num];

        double x_m_tilde = 2 * PI * ((x_m + grid.m_xLength) / grid.m_xLength);
        double y_m_tilde = 2 * PI * ((y_m + grid.m_yLength) / grid.m_yLength);
        double x_p_tilde = 2 * PI * ((x_p + grid.m_xLength) / grid.m_xLength);
        double y_p_tilde = 2 * PI * ((y_p + grid.m_yLength) / grid.m_yLength);

#pragma omp parallel for collapse(2) shared(grid, theta, theta_k, y_m_tilde, x_m_tilde, y_p_tilde, x_p_tilde) default(none)
        for (int i = 0; i < grid.m_xPoints; ++i)
        {
            for (int j = 0; j < grid.m_yPoints; ++j)
            {
                double x_tilde = 2 * PI * ((grid.m_xMesh[i][j] + grid.m_xLength) / grid.m_xLength);
                double y_tilde = 2 * PI * ((grid.m_yMesh[i][j] + grid.m_xLength) / grid.m_xLength);

                // Aux variables
                double Y_minus = y_tilde - y_m_tilde;
                double X_minus = x_tilde - x_m_tilde;
                double Y_plus = y_tilde - y_p_tilde;
                double X_plus = x_tilde - x_p_tilde;

                double heav_xp = heaviside(X_plus);
                double heav_xm = heaviside(X_minus);

                for (int k = -5; k < 6; ++k)
                {
                    theta_k[i][j] += std::atan(std::tanh((Y_minus + 2 * PI * k) / 2)
                                               * std::tan((X_minus - PI) / 2))
                                     - std::atan(std::tanh((Y_plus + 2 * PI * k) / 2)
                                                 * std::tan((X_plus - PI) / 2))
                                     + PI * (heav_xp - heav_xm);
                }

                theta_k[i][j] -= y_tilde * (x_p_tilde - x_m_tilde) / (2 * PI);
                theta[i][j] += theta_k[i][j];
            }
        }

    }

    std::cout << "Phase constructed!\n";

    return theta;
}


#endif //BECPP_PHASE_H
