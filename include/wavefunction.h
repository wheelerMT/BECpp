//
// Created by mattw on 05/01/2022.
//

#ifndef BECPP_WAVEFUNCTION_H
#define BECPP_WAVEFUNCTION_H

#include <cmath>
#include <complex>
#include <chrono>
#include <string>
#include <iostream>
#include <random>
#include "grid.h"
#include "constants.h"

// Define arrays
using complexVector_t = std::vector<std::complex<double>>;
using complexArray_t = std::vector<complexVector_t>;
using doubleArray_t = std::vector<std::vector<double>>;

class Wavefunction
{
private:
    void generateInitialState(const std::string &gs_phase);

public:
    // Wavefunction components arrays
    complexVector_t plus{};
    complexVector_t zero{};
    complexVector_t minus{};

    // Wavefunction k-space arrays
    complexVector_t plus_k{};
    complexVector_t zero_k{};
    complexVector_t minus_k{};

    // Constructor
    Wavefunction(Grid &t_grid, const std::string &gs_phase);

    // Member functions
    void add_noise(std::string const &components, double mean, double stddev);

    doubleArray_t density();

    // Reference to grid object
    Grid &grid;
};


#endif //BECPP_WAVEFUNCTION_H
