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

class Wavefunction1D
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
    // Wavefunction2D components arrays
    complexVector_t plus{};
    complexVector_t zero{};
    complexVector_t minus{};

    // Wavefunction2D k-space arrays
    complexVector_t plus_k{};
    complexVector_t zero_k{};
    complexVector_t minus_k{};

    // Reference to grid object
    Grid1D &grid;

    // Atom numbers
    double N_plus{};
    double N_zero{};
    double N_minus{};
    double N{};

    // Constructor
    Wavefunction1D(Grid1D &grid, const std::string &gs_phase);

    // FFT functions
    void fft();
    void ifft();
    void destroy_fft_plans();

    // Member functions
    void add_noise(std::string const &components, double mean, double stddev);

    std::vector<double> density();

    double atom_number();

    double component_atom_number(const std::string &component);

    void update_component_atom_num();

};

class Wavefunction2D
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
    // Wavefunction2D components arrays
    complexVector_t plus{};
    complexVector_t zero{};
    complexVector_t minus{};

    // Wavefunction2D k-space arrays
    complexVector_t plus_k{};
    complexVector_t zero_k{};
    complexVector_t minus_k{};

    // Reference to grid object
    Grid2D &grid;

    // Atom numbers
    double N_plus{};
    double N_zero{};
    double N_minus{};
    double N{};

    // Constructor
    Wavefunction2D(Grid2D &grid, const std::string &gs_phase);

    // FFT functions
    void fft();
    void ifft();
    void destroy_fft_plans();

    // Member functions
    void add_noise(std::string const &components, double mean, double stddev);

    void apply_phase(const doubleArray_t &phase_profile);

    doubleArray_t density();

    double atom_number();

    double component_atom_number(const std::string &component);

    void update_component_atom_num();

};


#endif //BECPP_WAVEFUNCTION_H
