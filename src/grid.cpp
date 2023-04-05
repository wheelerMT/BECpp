#include "grid.h"
#include "constants.h"
#include <algorithm>
#include <cmath>
#include <utility>

Grid1D::Grid1D(unsigned int nx, double dx) : m_gridPoints{nx}, m_gridSpacing{dx}
{
    Grid1D::constructGridParams();
    Grid1D::constructMesh();
}

void Grid1D::constructGridParams()
{
    m_fourierGridSpacing = PI / (m_gridPoints / 2. * m_gridSpacing);
    m_length = m_gridPoints * m_gridSpacing;
}

void Grid1D::constructMesh()
{
    m_mesh.m_xMesh.resize(m_gridPoints);
    m_mesh.m_xFourierMesh.resize(m_gridPoints);
    m_mesh.m_wavenumber.resize(m_gridPoints);

    for (int i = 0; i < m_gridPoints; ++i)
    {
        m_mesh.m_xMesh[i] = (i - m_gridPoints / 2.) * m_gridSpacing;
        m_mesh.m_xFourierMesh[i] =
                (i - m_gridPoints / 2.) * m_fourierGridSpacing;
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

void Grid1D::fftshift()
{
    std::ranges::rotate(m_mesh.m_xFourierMesh.begin(),
                        m_mesh.m_xFourierMesh.begin() + m_gridPoints / 2,
                        m_mesh.m_xFourierMesh.end());

    for (int i = 0; i < m_gridPoints; ++i)
    {
        m_mesh.m_wavenumber[i] = std::pow(m_mesh.m_xFourierMesh[i], 2);
    }
}

unsigned int Grid1D::shape() const { return m_gridPoints; }

double Grid1D::gridSpacing() const { return m_gridSpacing; }

double Grid1D::fourierGridSpacing() const { return m_fourierGridSpacing; }

double Grid1D::gridLength() const { return m_length; }

std::vector<double> Grid1D::wavenumber() const { return m_mesh.m_wavenumber; }

Grid2D::Grid2D(std::tuple<unsigned int, unsigned int> points,
               std::tuple<double, double> gridSpacing)
    : m_gridPoints{std::move(points)}, m_gridSpacing{std::move(gridSpacing)}
{
    Grid2D::constructGridParams();
    Grid2D::constructMesh();
}

void Grid2D::constructGridParams()
{
    unsigned int xPoints{std::get<0>(m_gridPoints)};
    unsigned int yPoints{std::get<1>(m_gridPoints)};
    double xGridSpacing{std::get<0>(m_gridSpacing)};
    double yGridSpacing{std::get<1>(m_gridSpacing)};

    m_fourierGridSpacing = {PI / (xPoints / 2. * xGridSpacing),
                            PI / (yPoints / 2. * yGridSpacing)};
    m_gridLength = {xPoints * xGridSpacing, yPoints * yGridSpacing};
}

void Grid2D::constructMesh()
{
    unsigned int xPoints{std::get<0>(m_gridPoints)};
    unsigned int yPoints{std::get<1>(m_gridPoints)};
    double xGridSpacing{std::get<0>(m_gridSpacing)};
    double yGridSpacing{std::get<1>(m_gridSpacing)};
    double xFourierGridSpacing{std::get<0>(m_fourierGridSpacing)};
    double yFourierGridSpacing{std::get<1>(m_fourierGridSpacing)};

    m_mesh.xMesh.resize(xPoints, std::vector<double>(yPoints));
    m_mesh.yMesh.resize(xPoints, std::vector<double>(yPoints));
    m_mesh.xFourierMesh.resize(xPoints, std::vector<double>(yPoints));
    m_mesh.yFourierMesh.resize(xPoints, std::vector<double>(yPoints));
    m_mesh.wavenumber.resize(xPoints, std::vector<double>(yPoints));

    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            m_mesh.xMesh[i][j] = (j - xPoints / 2.) * xGridSpacing;
            m_mesh.xFourierMesh[i][j] = (j - xPoints / 2.) * xFourierGridSpacing;
            m_mesh.yMesh[j][i] = (j - yPoints / 2.) * yGridSpacing;
            m_mesh.yFourierMesh[j][i] = (j - yPoints / 2.) * yFourierGridSpacing;
        }
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

void Grid2D::fftshift()
{
    // To implement
}

std::tuple<unsigned int, unsigned int> Grid2D::ndim() const
{
    return m_gridPoints;
}

std::tuple<double, double> Grid2D::gridSpacing() const { return m_gridSpacing; }

std::tuple<double, double> Grid2D::fourierGridSpacing() const
{
    return m_fourierGridSpacing;
}

std::tuple<double, double> Grid2D::gridLength() const { return m_gridLength; }

doubleArray_t Grid2D::wavenumber() const { return m_mesh.wavenumber; }

void Grid3D::constructGridParams()
{
    m_gridSpacingProduct = m_xGridSpacing * m_yGridSpacing * m_zGridSpacing;

    m_xFourierGridSpacing = PI / (m_xPoints / 2. * m_xGridSpacing);
    m_yFourierGridSpacing = PI / (m_yPoints / 2. * m_yGridSpacing);
    m_zFourierGridSpacing = PI / (m_zPoints / 2. * m_zGridSpacing);

    m_xLength = m_xPoints * m_xGridSpacing;
    m_yLength = m_yPoints * m_yGridSpacing;
    m_zLength = m_zPoints * m_zGridSpacing;
}

void Grid3D::constructMesh()
{
    // Set grid sizes
    m_xMesh.resize(m_xPoints,
                   std::vector<std::vector<double>>(
                           m_yPoints, std::vector<double>(m_zPoints)));
    m_yMesh.resize(m_xPoints,
                   std::vector<std::vector<double>>(
                           m_yPoints, std::vector<double>(m_zPoints)));
    m_xFourierMesh.resize(m_xPoints,
                          std::vector<std::vector<double>>(
                                  m_yPoints, std::vector<double>(m_zPoints)));
    m_yFourierMesh.resize(m_xPoints,
                          std::vector<std::vector<double>>(
                                  m_yPoints, std::vector<double>(m_zPoints)));
    m_wavenumber.resize(m_xPoints,
                        std::vector<std::vector<double>>(
                                m_yPoints, std::vector<double>(m_zPoints)));

    // Construct grids
    for (int i = 0; i < m_xPoints; ++i)
    {
        for (int j = 0; j < m_yPoints; ++j)
        {
            for (int k = 0; k < m_zPoints; ++k)
            {
                m_xMesh[i][j][k] = (j - m_xPoints / 2.) * m_xGridSpacing;
                m_xFourierMesh[i][j][k] =
                        (j - m_xPoints / 2.) * m_xFourierGridSpacing;
                m_yMesh[j][i][k] = (j - m_yPoints / 2.) * m_yGridSpacing;
                m_yFourierMesh[j][i][k] =
                        (j - m_yPoints / 2.) * m_yFourierGridSpacing;
            }
        }
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

void Grid3D::fftshift()
{
    // TODO: Implement using std::rotate
}

Grid3D::Grid3D(std::tuple<unsigned int, unsigned int, unsigned int> points,
               std::tuple<double, double, double> gridSpacing)
    : m_xPoints{std::get<0>(points)}, m_yPoints{std::get<1>(points)},
      m_zPoints{std::get<2>(points)}, m_xGridSpacing{std::get<0>(gridSpacing)},
      m_yGridSpacing{std::get<1>(gridSpacing)},
      m_zGridSpacing{std::get<2>(gridSpacing)}
{
    Grid3D::constructGridParams();
    Grid3D::constructMesh();
}

unsigned int Grid3D::xPoints() const { return m_xPoints; }

unsigned int Grid3D::yPoints() const { return m_yPoints; }

unsigned int Grid3D::zPoints() const { return m_zPoints; }

double Grid3D::xGridSpacing() const { return m_xGridSpacing; }

double Grid3D::yGridSpacing() const { return m_yGridSpacing; }

double Grid3D::zGridSpacing() const { return m_zGridSpacing; }

double Grid3D::gridSpacingProduct() const { return m_gridSpacingProduct; }

double Grid3D::xFourierGridSpacing() const { return m_xFourierGridSpacing; }

double Grid3D::yFourierGridSpacing() const { return m_yFourierGridSpacing; }

double Grid3D::zFourierGridSpacing() const { return m_zFourierGridSpacing; }

double Grid3D::xLength() const { return m_xLength; }

double Grid3D::yLength() const { return m_yLength; }

double Grid3D::zLength() const { return m_zLength; }

std::vector<doubleArray_t> Grid3D::wavenumber() const { return m_wavenumber; }
