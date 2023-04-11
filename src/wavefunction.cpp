#include "wavefunction.h"

#include <utility>

Wavefunction1D::Wavefunction1D(const Grid1D& grid) : m_grid{grid}
{
    m_component.resize(grid.shape());
    m_fourierComponent.resize(grid.shape());

    createFFTPlans(grid);
}

void Wavefunction1D::createFFTPlans(const Grid1D& grid)
{
    m_plans.plan_forward = fftw_plan_dft_1d(
            static_cast<int>(grid.shape()),
            reinterpret_cast<fftw_complex*>(&m_component[0]),
            reinterpret_cast<fftw_complex*>(&m_fourierComponent[0]),
            FFTW_FORWARD, FFTW_MEASURE);
    m_plans.plan_backward = fftw_plan_dft_1d(
            static_cast<int>(grid.shape()),
            reinterpret_cast<fftw_complex*>(&m_fourierComponent[0]),
            reinterpret_cast<fftw_complex*>(&m_component[0]), FFTW_BACKWARD,
            FFTW_MEASURE);
}

Wavefunction1D::~Wavefunction1D() { destroyFFTPlans(); }

void Wavefunction1D::destroyFFTPlans() const
{
    fftw_destroy_plan(m_plans.plan_forward);
    fftw_destroy_plan(m_plans.plan_backward);
}

std::vector<std::complex<double>>& Wavefunction1D::component()
{
    return m_component;
}

std::vector<std::complex<double>>& Wavefunction1D::fourierComponent()
{
    return m_fourierComponent;
}

std::vector<double> Wavefunction1D::density() const
{
    std::vector<double> density{};
    density.resize(m_grid.shape());
    for (int i = 0; i < m_grid.shape(); ++i)
    {
        density[i] = std::pow(std::abs(m_component[i]), 2);
    }

    return density;
}

double Wavefunction1D::atomNumber() const { return m_atomNumber; }

void Wavefunction1D::fft() const {
    fftw_execute(m_plans.plan_forward);
}

void Wavefunction1D::ifft() const {
    fftw_execute(m_plans.plan_backward);
}

void Wavefunction1D::setComponent(std::vector<std::complex<double>> component)
{
    m_component = std::move(component);
}
