#include "grid.h"
#include "constants.h"
#include <algorithm>
#include <cmath>

void Grid2D::constructGridParams()
{
    // Calculate k-space grid spacing
    m_xFourierGridSpacing = PI / (m_xPoints / 2. * m_xGridSpacing);
    m_yFourierGridSpacing = PI / (m_yPoints / 2. * m_yGridSpacing);

    // Sets length of sides of box
    m_xLength = m_xPoints * m_xGridSpacing;
    m_yLength = m_yPoints * m_yGridSpacing;
}

void Grid2D::constructGrids()
{
    // Set grid sizes
    m_xMesh.resize(m_xPoints, std::vector<double>(m_yPoints));
    m_yMesh.resize(m_xPoints, std::vector<double>(m_yPoints));
    m_xFourierMesh.resize(m_xPoints, std::vector<double>(m_yPoints));
    m_yFourierMesh.resize(m_xPoints, std::vector<double>(m_yPoints));
    m_wavenumber.resize(m_xPoints, std::vector<double>(m_yPoints));

    // Construct grids
    for (int i = 0; i < m_xPoints; ++i)
    {
        for (int j = 0; j < m_yPoints; ++j)
        {
            m_xMesh[i][j] = (j - m_xPoints / 2.) * m_xGridSpacing;
            m_xFourierMesh[i][j] = (j - m_xPoints / 2.) * m_xFourierGridSpacing;
            m_yMesh[j][i] = (j - m_yPoints / 2.) * m_yGridSpacing;
            m_yFourierMesh[j][i] = (j - m_yPoints / 2.) * m_yFourierGridSpacing;
            m_wavenumber[i][j] = std::pow(m_xFourierMesh[i][j], 2) +
                                 std::pow(m_yFourierMesh[i][j], 2);
        }
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

void Grid2D::fftshift()
{
    std::vector<std::vector<double>> Kx_copy = m_xFourierMesh;
    std::vector<std::vector<double>> Ky_copy = m_yFourierMesh;

    // Reverse each row
    for (int i = 0; i < m_xPoints; i++)
    {
        for (int j = 0; j < m_yPoints; j++)
        {
            if (j < m_xPoints / 2)
            {
                m_xFourierMesh[i][j] = Kx_copy[i][m_xPoints / 2 + j];
                m_yFourierMesh[i][j] = Ky_copy[i][m_xPoints / 2 + j];
            } else if (j >= m_xPoints / 2)
            {
                m_xFourierMesh[i][j] = Kx_copy[i][j - m_xPoints / 2];
                m_yFourierMesh[i][j] = Ky_copy[i][j - m_xPoints / 2];
            }
        }
    }

    // Update copy arrays
    Kx_copy = m_xFourierMesh;
    Ky_copy = m_yFourierMesh;

    // Reverse each column
    for (int i = 0; i < m_xPoints; i++)
    {
        for (int j = 0; j < m_yPoints; j++)
        {
            if (j < m_xPoints / 2)
            {
                m_xFourierMesh[j][i] = Kx_copy[m_xPoints / 2 + j][i];
                m_yFourierMesh[j][i] = Ky_copy[m_xPoints / 2 + j][i];
            } else if (j >= m_xPoints / 2)
            {
                m_xFourierMesh[j][i] = Kx_copy[j - m_xPoints / 2][i];
                m_yFourierMesh[j][i] = Ky_copy[j - m_xPoints / 2][i];
            }
        }
    }

    // Update wave vector
    for (int i = 0; i < m_xPoints; ++i)
    {
        for (int j = 0; j < m_yPoints; ++j)
        {
            m_wavenumber[i][j] = std::pow(m_xFourierMesh[i][j], 2) +
                                 std::pow(m_yFourierMesh[i][j], 2);
        }
    }
}

Grid2D::Grid2D(unsigned int nx, unsigned int ny, double dx, double dy)
    : m_xPoints{nx}, m_yPoints{ny}, m_xGridSpacing{dx}, m_yGridSpacing{dy},
      m_gridSpacingProduct{dx * dy}
{
    Grid2D::constructGridParams();
    Grid2D::constructGrids();
}

unsigned int Grid2D::xPoints() const { return m_xPoints; }
unsigned int Grid2D::yPoints() const { return m_yPoints; }
double Grid2D::xGridSpacing() const { return m_xGridSpacing; }
double Grid2D::yGridSpacing() const { return m_yGridSpacing; }
double Grid2D::gridSpacingProduct() const { return m_gridSpacingProduct; }
double Grid2D::xFourierGridSpacing() const { return m_xFourierGridSpacing; }
double Grid2D::yFourierGridSpacing() const { return m_yFourierGridSpacing; }
double Grid2D::xLength() const { return m_xLength; }
double Grid2D::yLength() const { return m_yLength; }
doubleArray_t Grid2D::wavenumber() const { return m_wavenumber; }

void Grid1D::constructGridParams()
{
    m_xFourierGridSpacing = PI / (m_xPoints / 2. * m_xGridSpacing);
    m_xLength = m_xPoints * m_xGridSpacing;
}

void Grid1D::constructGrids()
{
    // Set grid sizes
    m_xMesh.resize(m_xPoints);
    m_xFourierMesh.resize(m_xPoints);
    m_wavenumber.resize(m_xPoints);

    // Construct grids
    for (int i = 0; i < m_xPoints; ++i)
    {
        m_xMesh[i] = (i - m_xPoints / 2.) * m_xGridSpacing;
        m_xFourierMesh[i] = (i - m_xPoints / 2.) * m_xFourierGridSpacing;
        m_wavenumber[i] = std::pow(m_xFourierMesh[i], 2);
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

Grid1D::Grid1D(unsigned int nx, double dx) : m_xPoints{nx}, m_xGridSpacing{dx}
{
    Grid1D::constructGridParams();
    Grid1D::constructGrids();
}

void Grid1D::fftshift()
{
    std::ranges::rotate(m_xFourierMesh.begin(),
                        m_xFourierMesh.begin() + m_xPoints / 2,
                        m_xFourierMesh.end());
    m_wavenumber = m_xFourierMesh;
}
unsigned int Grid1D::xPoints() const { return m_xPoints; }
double Grid1D::xGridSpacing() const { return m_xGridSpacing; }
double Grid1D::xFourierGridSpacing() const { return m_xFourierGridSpacing; }
double Grid1D::xLength() const { return m_xLength; }
std::vector<double> Grid1D::wavenumber() const { return m_wavenumber; }
