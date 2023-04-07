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
            m_mesh.xFourierMesh[i][j] =
                    (j - xPoints / 2.) * xFourierGridSpacing;
            m_mesh.yMesh[j][i] = (j - yPoints / 2.) * yGridSpacing;
            m_mesh.yFourierMesh[j][i] =
                    (j - yPoints / 2.) * yFourierGridSpacing;
        }
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

void Grid2D::fftshift()
{
    // To implement
}

std::tuple<unsigned int, unsigned int> Grid2D::shape() const
{
    return m_gridPoints;
}

std::tuple<double, double> Grid2D::gridSpacing() const { return m_gridSpacing; }

std::tuple<double, double> Grid2D::fourierGridSpacing() const
{
    return m_fourierGridSpacing;
}

std::tuple<double, double> Grid2D::gridLength() const { return m_gridLength; }

vector2D_t Grid2D::wavenumber() const { return m_mesh.wavenumber; }

Grid3D::Grid3D(std::tuple<unsigned int, unsigned int, unsigned int> points,
               std::tuple<double, double, double> gridSpacing)
    : m_gridPoints{std::move(points)}, m_gridSpacing{std::move(gridSpacing)}
{
    Grid3D::constructGridParams();
    Grid3D::constructMesh();
}

void Grid3D::constructGridParams()
{
    auto [xPoints, yPoints, zPoints] = shape();
    auto [xGridSpacing, yGridSpacing, zGridSpacing] = gridSpacing();

    m_fourierGridSpacing = {PI / (xPoints / 2. * xGridSpacing),
                            PI / (yPoints / 2. * yGridSpacing),
                            PI / (zPoints / 2. * zGridSpacing)};
    m_gridLength = {xPoints * xGridSpacing, yPoints * zGridSpacing,
                    yPoints * zGridSpacing};
}

void Grid3D::constructMesh()
{
    auto [xPoints, yPoints, zPoints] = shape();
    auto [xGridSpacing, yGridSpacing, zGridSpacing] = gridSpacing();
    auto [xFourierGridSpacing, yFourierGridSpacing, zFourierGridSpacing] =
            fourierGridSpacing();

    m_mesh.xMesh.resize(xPoints,
                        vector2D_t(yPoints, std::vector<double>(zPoints)));
    m_mesh.yMesh.resize(xPoints,
                        vector2D_t(yPoints, std::vector<double>(zPoints)));
    m_mesh.zMesh.resize(xPoints,
                        vector2D_t(yPoints, std::vector<double>(zPoints)));
    m_mesh.xFourierMesh.resize(
            xPoints, vector2D_t(yPoints, std::vector<double>(zPoints)));
    m_mesh.yFourierMesh.resize(
            xPoints, vector2D_t(yPoints, std::vector<double>(zPoints)));
    m_mesh.zFourierMesh.resize(
            xPoints, vector2D_t(yPoints, std::vector<double>(zPoints)));
    m_mesh.wavenumber.resize(xPoints,
                             vector2D_t(yPoints, std::vector<double>(zPoints)));

    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int k = 0; k < zPoints; ++k)
            {
                m_mesh.xMesh[i][j][k] = (j - xPoints / 2.) * xGridSpacing;
                m_mesh.xFourierMesh[i][j][k] =
                        (j - xPoints / 2.) * xFourierGridSpacing;
                m_mesh.yMesh[j][i][k] = (j - yPoints / 2.) * yGridSpacing;
                m_mesh.yFourierMesh[j][i][k] =
                        (j - yPoints / 2.) * yFourierGridSpacing;
                m_mesh.zMesh[k][j][i] = (k - zPoints / 2.) * zGridSpacing;
                m_mesh.zFourierMesh[k][j][i] =
                        (k - zPoints / 2.) * zFourierGridSpacing;
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

std::tuple<unsigned int, unsigned int, unsigned int> Grid3D::shape() const
{
    return m_gridPoints;
}

std::tuple<double, double, double> Grid3D::gridSpacing() const
{
    return m_gridSpacing;
}

std::tuple<double, double, double> Grid3D::fourierGridSpacing() const
{
    return m_fourierGridSpacing;
}

std::tuple<double, double, double> Grid3D::gridLength() const
{
    return m_gridLength;
}

vector3D_t Grid3D::wavenumber() const { return m_mesh.wavenumber; }
