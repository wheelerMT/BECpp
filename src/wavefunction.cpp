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

void Wavefunction1D::destroyFFTPlans()
{
    fftw_destroy_plan(m_plans.plan_forward);
    fftw_destroy_plan(m_plans.plan_backward);
}

complexVector_t& Wavefunction1D::component()
{
    return m_component;
}

complexVector_t& Wavefunction1D::fourierComponent()
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

void Wavefunction1D::fft() { fftw_execute(m_plans.plan_forward); }

void Wavefunction1D::ifft()
{
    fftw_execute(m_plans.plan_backward);

    // Renormalise wavefunction
    for (size_t i = 0; i < m_grid.shape(); i++)
    {
        m_component[i] /= m_grid.shape();
    }
}

void Wavefunction1D::updateAtomNumber()
{
    for (size_t i = 0; i < m_grid.shape(); i++)
    {
        m_atomNumber +=
                std::pow(std::abs(density()[i]), 2) * m_grid.gridSpacing();
    }
}

void Wavefunction1D::setComponent(complexVector_t& component)
{
    m_component = component;

    // Update the Fourier-space wavefunction and atom number
    fft();
    updateAtomNumber();
}

Wavefunction2D::Wavefunction2D(const Grid2D& grid) : m_grid{grid}
{
    auto [xPoints, yPoints] = grid.shape();
    m_component.resize(xPoints * yPoints);
    m_fourierComponent.resize(xPoints * yPoints);

    createFFTPlans(grid);
}

void Wavefunction2D::createFFTPlans(const Grid2D& grid)
{
    auto [xPoints, yPoints] = grid.shape();
    m_plans.plan_forward = fftw_plan_dft_2d(
            static_cast<int>(xPoints), static_cast<int>(yPoints),
            reinterpret_cast<fftw_complex*>(&m_component[0]),
            reinterpret_cast<fftw_complex*>(&m_fourierComponent[0]),
            FFTW_FORWARD, FFTW_MEASURE);
    m_plans.plan_backward = fftw_plan_dft_2d(
            static_cast<int>(xPoints), static_cast<int>(yPoints),
            reinterpret_cast<fftw_complex*>(&m_fourierComponent[0]),
            reinterpret_cast<fftw_complex*>(&m_component[0]), FFTW_BACKWARD,
            FFTW_MEASURE);
}

Wavefunction2D::~Wavefunction2D() { destroyFFTPlans(); }

void Wavefunction2D::destroyFFTPlans()
{
    fftw_destroy_plan(m_plans.plan_forward);
    fftw_destroy_plan(m_plans.plan_backward);
}

complexVector_t& Wavefunction2D::component()
{
    return m_component;
}

complexVector_t& Wavefunction2D::fourierComponent()
{
    return m_fourierComponent;
}

std::vector<double> Wavefunction2D::density() const
{
    std::vector<double> density{};
    auto [xPoints, yPoints] = m_grid.shape();
    density.resize(xPoints * yPoints);
    for (size_t i = 0; i < xPoints; i++)
    {
        for (size_t j = 0; j < yPoints; j++)
        {
            density[j + i * yPoints] =
                    std::pow(std::abs(m_component[j + i * yPoints]), 2);
        }
    }

    return density;
}

double Wavefunction2D::atomNumber() const { return m_atomNumber; }

void Wavefunction2D::fft() { fftw_execute(m_plans.plan_forward); }

void Wavefunction2D::ifft()
{
    fftw_execute(m_plans.plan_backward);

    // Renormalise wavefunction
    auto [xPoints, yPoints] = m_grid.shape();
    for (size_t i = 0; i < xPoints; i++)
    {
        for (size_t j = 0; j < yPoints; j++)
        {
            m_component[j + i * yPoints] /= xPoints * yPoints;
        }
    }
}

void Wavefunction2D::updateAtomNumber()
{
    auto [xPoints, yPoints] = m_grid.shape();
    auto [xGridSpacing, yGridSpacing] = m_grid.gridSpacing();
    for (size_t i = 0; i < xPoints; i++)
    {
        for (size_t j = 0; j < yPoints; j++)
        {
            m_atomNumber += std::pow(std::abs(density()[j + i * yPoints]), 2) *
                            xGridSpacing * yGridSpacing;
        }
    }
}

void Wavefunction2D::setComponent(complexVector_t& component)
{
    m_component = component;

    // Update the Fourier-space wavefunction and atom number
    fft();
    updateAtomNumber();
}