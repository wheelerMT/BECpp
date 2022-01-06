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
#include <algorithm>
#include "grid.h"
#include "constants.h"
#include "fftw3.h"

// Define arrays
using complexVector_t = std::vector<std::complex<double>>;
using complexArray_t = std::vector<complexVector_t>;
using doubleArray_t = std::vector<std::vector<double>>;

class Wavefunction
{
private:
    void generateInitialState(const std::string &gs_phase);
    void generateFFTPlans();

    // FFT plans
    fftw_plan forward_plus{};
    fftw_plan forward_zero{};
    fftw_plan forward_minus{};
    fftw_plan backward_plus{};
    fftw_plan backward_zero{};
    fftw_plan backward_minus{};

public:
    // Wavefunction components arrays
    complexVector_t plus{};
    complexVector_t zero{};
    complexVector_t minus{};

    // Wavefunction k-space arrays
    complexVector_t plus_k{};
    complexVector_t zero_k{};
    complexVector_t minus_k{};

    // Reference to grid object
    Grid &grid;

    // Constructor
    Wavefunction(Grid &grid, const std::string &gs_phase);

    // FFT functions
    void fft();
    void ifft();

    // Member functions
    void add_noise(std::string const &components, double mean, double stddev);

    doubleArray_t density();


};


#endif //BECPP_WAVEFUNCTION_H
