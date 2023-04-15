#include "wavefunction.h"

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

complexVector_t& Wavefunction1D::component() { return m_component; }

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

void Wavefunction1D::fft() const { fftw_execute(m_plans.plan_forward); }

void Wavefunction1D::ifft()
{
    fftw_execute(m_plans.plan_backward);

    // Renormalise wavefunction
    for (int i = 0; i < m_grid.shape(); ++i)
    {
        m_component[i] /= static_cast<double>(m_grid.shape());
    }
}

void Wavefunction1D::updateAtomNumber()
{
    for (int i = 0; i < m_grid.shape(); ++i)
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

void Wavefunction2D::destroyFFTPlans() const
{
    fftw_destroy_plan(m_plans.plan_forward);
    fftw_destroy_plan(m_plans.plan_backward);
}

complexVector_t& Wavefunction2D::component() { return m_component; }

complexVector_t& Wavefunction2D::fourierComponent()
{
    return m_fourierComponent;
}

std::vector<double> Wavefunction2D::density() const
{
    std::vector<double> density{};
    auto [xPoints, yPoints] = m_grid.shape();
    density.resize(xPoints * yPoints);
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            density[j + i * yPoints] =
                    std::pow(std::abs(m_component[j + i * yPoints]), 2);
        }
    }

    return density;
}

double Wavefunction2D::atomNumber() const { return m_atomNumber; }

void Wavefunction2D::fft() const { fftw_execute(m_plans.plan_forward); }

void Wavefunction2D::ifft()
{
    fftw_execute(m_plans.plan_backward);

    // Renormalise wavefunction
    auto [xPoints, yPoints] = m_grid.shape();
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            m_component[j + i * yPoints] /=
                    static_cast<double>(xPoints * yPoints);
        }
    }
}

void Wavefunction2D::updateAtomNumber()
{
    auto [xPoints, yPoints] = m_grid.shape();
    auto [xGridSpacing, yGridSpacing] = m_grid.gridSpacing();
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
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

Wavefunction3D::Wavefunction3D(const Grid3D& grid) : m_grid{grid}
{
    auto [xPoints, yPoints, zPoints] = grid.shape();
    m_component.resize(xPoints * yPoints * zPoints);
    m_fourierComponent.resize(xPoints * yPoints * zPoints);

    createFFTPlans(grid);
}

void Wavefunction3D::createFFTPlans(const Grid3D& grid)
{
    auto [xPoints, yPoints, zPoints] = grid.shape();
    m_plans.plan_forward = fftw_plan_dft_3d(
            static_cast<int>(xPoints), static_cast<int>(yPoints),
            static_cast<int>(zPoints),
            reinterpret_cast<fftw_complex*>(&m_component[0]),
            reinterpret_cast<fftw_complex*>(&m_fourierComponent[0]),
            FFTW_FORWARD, FFTW_MEASURE);
    m_plans.plan_backward = fftw_plan_dft_3d(
            static_cast<int>(xPoints), static_cast<int>(yPoints),
            static_cast<int>(zPoints),
            reinterpret_cast<fftw_complex*>(&m_fourierComponent[0]),
            reinterpret_cast<fftw_complex*>(&m_component[0]), FFTW_BACKWARD,
            FFTW_MEASURE);
}

Wavefunction3D::~Wavefunction3D() { destroyFFTPlans(); }

void Wavefunction3D::destroyFFTPlans() const
{
    fftw_destroy_plan(m_plans.plan_forward);
    fftw_destroy_plan(m_plans.plan_backward);
}

complexVector_t& Wavefunction3D::component() { return m_component; }

complexVector_t& Wavefunction3D::fourierComponent()
{
    return m_fourierComponent;
}

std::vector<double> Wavefunction3D::density() const
{
    std::vector<double> density{};
    auto [xPoints, yPoints, zPoints] = m_grid.shape();
    density.resize(xPoints * yPoints * zPoints);
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                auto index = k + zPoints * (j + yPoints * i);
                density[index] = std::pow(std::abs(m_component[index]), 2);
            }
        }
    }

    return density;
}

double Wavefunction3D::atomNumber() const { return m_atomNumber; }

void Wavefunction3D::fft() const { fftw_execute(m_plans.plan_forward); }

void Wavefunction3D::ifft()
{
    fftw_execute(m_plans.plan_backward);

    // Renormalise wavefunction
    auto [xPoints, yPoints, zPoints] = m_grid.shape();
    auto n = static_cast<double>(xPoints * yPoints * zPoints);
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                auto index = k + zPoints * (j + yPoints * i);
                m_component[index] /= n;
            }
        }
    }
}

void Wavefunction3D::updateAtomNumber()
{
    auto [xPoints, yPoints, zPoints] = m_grid.shape();
    auto [xGridSpacing, yGridSpacing, zGridSpacing] = m_grid.gridSpacing();
    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                auto index = k + zPoints * (j + yPoints * i);
                m_atomNumber += std::pow(std::abs(density()[index]), 2) *
                                xGridSpacing * yGridSpacing * zGridSpacing;
            }
        }
    }
}

void Wavefunction3D::setComponent(complexVector_t& component)
{
    m_component = component;

    // Update the Fourier-space wavefunction and atom number
    fft();
    updateAtomNumber();
}