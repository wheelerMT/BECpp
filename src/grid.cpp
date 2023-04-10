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
    auto xPoints = shape();
    auto xGridSpacing = gridSpacing();
    auto xFourierGridSpacing = fourierGridSpacing();
    m_mesh.xMesh.resize(xPoints);
    m_mesh.xFourierMesh.resize(xPoints);
    m_mesh.wavenumber.resize(xPoints);

    for (int i = 0; i < xPoints; ++i)
    {
        m_mesh.xMesh[i] = (i - xPoints / 2.) * xGridSpacing;
        if (i < xPoints / 2)
        {
            m_mesh.xFourierMesh[i] = i * xFourierGridSpacing;
        } else
        {
            m_mesh.xFourierMesh[i] = (i - xPoints) * xFourierGridSpacing;
        }

        m_mesh.wavenumber[i] = std::pow(m_mesh.xFourierMesh[i], 2);
    }
}

unsigned int Grid1D::shape() const { return m_gridPoints; }

double Grid1D::gridSpacing() const { return m_gridSpacing; }

double Grid1D::fourierGridSpacing() const { return m_fourierGridSpacing; }

double Grid1D::gridLength() const { return m_length; }

std::vector<double> Grid1D::wavenumber() const { return m_mesh.wavenumber; }

Grid2D::Grid2D(std::tuple<unsigned int, unsigned int> points,
               std::tuple<double, double> gridSpacing)
    : m_gridPoints{std::move(points)}, m_gridSpacing{std::move(gridSpacing)}
{
    Grid2D::constructGridParams();
    Grid2D::constructMesh();
}

void Grid2D::constructGridParams()
{
    auto [xPoints, yPoints] = shape();
    auto [xGridSpacing, yGridSpacing] = gridSpacing();

    m_fourierGridSpacing = {PI / (xPoints / 2. * xGridSpacing),
                            PI / (yPoints / 2. * yGridSpacing)};
    m_gridLength = {xPoints * xGridSpacing, yPoints * yGridSpacing};
}

void Grid2D::constructMesh()
{
    auto [xPoints, yPoints] = shape();
    auto [xGridSpacing, yGridSpacing] = gridSpacing();
    auto [xFourierGridSpacing, yFourierGridSpacing] = fourierGridSpacing();

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
            m_mesh.yMesh[i][j] = (j - yPoints / 2.) * yGridSpacing;
            if (i < xPoints / 2)
            {
                m_mesh.xFourierMesh[i][j] = i * xFourierGridSpacing;
                m_mesh.yFourierMesh[i][j] = j * yFourierGridSpacing;
            } else
            {
                m_mesh.xFourierMesh[i][j] = (i - xPoints) * xFourierGridSpacing;
                m_mesh.yFourierMesh[i][j] = (j - yPoints) * yFourierGridSpacing;
            }

            m_mesh.wavenumber[i][j] =
                    std::pow(m_mesh.xFourierMesh[i][j], 2) +
                    std::pow(m_mesh.yFourierMesh[i][j], 2);
        }
    }
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

    for (int k = 0; k < zPoints; ++k)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int i = 0; i < xPoints; ++i)
            {
                m_mesh.xMesh[i][j][k] = (i - xPoints / 2.) * xGridSpacing;
                m_mesh.yMesh[i][j][k] = (j - yPoints / 2.) * yGridSpacing;
                m_mesh.zMesh[i][j][k] = (k - zPoints / 2.) * zGridSpacing;

                if (i < xPoints / 2)
                {
                    m_mesh.xFourierMesh[i][j][k] = i * xFourierGridSpacing;
                    m_mesh.yFourierMesh[i][j][k] = j * yFourierGridSpacing;
                    m_mesh.zFourierMesh[i][j][k] = k * yFourierGridSpacing;
                } else
                {
                    m_mesh.xFourierMesh[i][j][k] =
                            (i - xPoints) * xFourierGridSpacing;
                    m_mesh.yFourierMesh[i][j][k] =
                            (j - yPoints) * yFourierGridSpacing;
                    m_mesh.zFourierMesh[i][j][k] =
                            (k - zPoints) * yFourierGridSpacing;
                }

                m_mesh.wavenumber[i][j][k] =
                        std::pow(m_mesh.xFourierMesh[i][j][k], 2) +
                        std::pow(m_mesh.yFourierMesh[i][j][k], 2) +
                        std::pow(m_mesh.zFourierMesh[i][j][k], 2);
            }
        }
    }
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
