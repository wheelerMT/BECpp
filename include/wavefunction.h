#ifndef BECPP_WAVEFUNCTION_H
#define BECPP_WAVEFUNCTION_H

#include <cmath>
#include <complex>
#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include "grid.h"
#include "constants.h"
#include "fftw3.h"

struct FFTPlans
{
    fftw_plan plan_forward{};
    fftw_plan plan_backward{};
};

class Wavefunction1D
{
private:
    const Grid1D& m_grid;
    FFTPlans m_plans{};
    std::vector<std::complex<double>> m_component{};
    std::vector<std::complex<double>> m_fourierComponent{};
    double m_atomNumber{};

    void createFFTPlans(const Grid1D& grid);
    void destroyFFTPlans() const;

public:
    explicit Wavefunction1D(const Grid1D& grid);
    ~Wavefunction1D();

    [[nodiscard]] std::vector<std::complex<double>>& component();
    [[nodiscard]] std::vector<std::complex<double>>& fourierComponent();
    [[nodiscard]] std::vector<double> density() const;
    [[nodiscard]] double atomNumber() const;

    void fft() const;
    void ifft() const;
    void setComponent(std::vector<std::complex<double>> component);
};

#endif //BECPP_WAVEFUNCTION_H
