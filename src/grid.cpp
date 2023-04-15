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

void resizeMesh1D(Mesh1D& mesh, unsigned int xPoints)
{
    mesh.xMesh.resize(xPoints);
    mesh.xFourierMesh.resize(xPoints);
    mesh.wavenumber.resize(xPoints);
}

void Grid1D::constructMesh()
{
    auto xPoints = shape();
    auto xGridSpacing = gridSpacing();
    auto xFourierGridSpacing = fourierGridSpacing();
    resizeMesh1D(m_mesh, xPoints);

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

void resizeMesh2D(Mesh2D& mesh, unsigned int xPoints, unsigned int yPoints)
{
    mesh.xMesh.resize(xPoints * yPoints);
    mesh.yMesh.resize(xPoints * yPoints);
    mesh.xFourierMesh.resize(xPoints * yPoints);
    mesh.yFourierMesh.resize(xPoints * yPoints);
    mesh.wavenumber.resize(xPoints * yPoints);
}

void Grid2D::constructMesh()
{
    auto [xPoints, yPoints] = shape();
    auto [xGridSpacing, yGridSpacing] = gridSpacing();
    auto [xFourierGridSpacing, yFourierGridSpacing] = fourierGridSpacing();
    resizeMesh2D(m_mesh, xPoints, yPoints);

    for (int i = 0; i < xPoints; ++i)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            auto index = j + i * yPoints;
            m_mesh.xMesh[index] = (j - xPoints / 2.) * xGridSpacing;
            m_mesh.yMesh[index] = (j - yPoints / 2.) * yGridSpacing;
            if (i < xPoints / 2)
            {
                m_mesh.xFourierMesh[index] = i * xFourierGridSpacing;
                m_mesh.yFourierMesh[index] = j * yFourierGridSpacing;
            } else
            {
                m_mesh.xFourierMesh[index] =
                        (i - xPoints) * xFourierGridSpacing;
                m_mesh.yFourierMesh[index] =
                        (j - yPoints) * yFourierGridSpacing;
            }

            m_mesh.wavenumber[index] = std::pow(m_mesh.xFourierMesh[index], 2) +
                                       std::pow(m_mesh.yFourierMesh[index], 2);
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

std::vector<double> Grid2D::wavenumber() const { return m_mesh.wavenumber; }

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

void resizeMesh3D(Mesh3D& mesh, unsigned int xPoints, unsigned int yPoints,
                  unsigned int zPoints)
{
    mesh.xMesh.resize(xPoints * yPoints * zPoints);
    mesh.yMesh.resize(xPoints * yPoints * zPoints);
    mesh.zMesh.resize(xPoints * yPoints * zPoints);
    mesh.xFourierMesh.resize(xPoints * yPoints * zPoints);
    mesh.yFourierMesh.resize(xPoints * yPoints * zPoints);
    mesh.zFourierMesh.resize(xPoints * yPoints * zPoints);
    mesh.wavenumber.resize(xPoints * yPoints * zPoints);
}

void Grid3D::constructMesh()
{
    auto [xPoints, yPoints, zPoints] = shape();
    auto [xGridSpacing, yGridSpacing, zGridSpacing] = gridSpacing();
    auto [xFourierGridSpacing, yFourierGridSpacing, zFourierGridSpacing] =
            fourierGridSpacing();
    resizeMesh3D(m_mesh, xPoints, yPoints, zPoints);

    for (int k = 0; k < zPoints; ++k)
    {
        for (int j = 0; j < yPoints; ++j)
        {
            for (int i = 0; i < xPoints; ++i)
            {
                auto index = k + zPoints * (j + yPoints * i);
                m_mesh.xMesh[index] = (i - xPoints / 2.) * xGridSpacing;
                m_mesh.yMesh[index] = (j - yPoints / 2.) * yGridSpacing;
                m_mesh.zMesh[index] = (k - zPoints / 2.) * zGridSpacing;

                if (i < xPoints / 2)
                {
                    m_mesh.xFourierMesh[index] = i * xFourierGridSpacing;
                    m_mesh.yFourierMesh[index] = j * yFourierGridSpacing;
                    m_mesh.zFourierMesh[index] = k * yFourierGridSpacing;
                } else
                {
                    m_mesh.xFourierMesh[index] =
                            (i - xPoints) * xFourierGridSpacing;
                    m_mesh.yFourierMesh[index] =
                            (j - yPoints) * yFourierGridSpacing;
                    m_mesh.zFourierMesh[index] =
                            (k - zPoints) * yFourierGridSpacing;
                }

                m_mesh.wavenumber[index] =
                        std::pow(m_mesh.xFourierMesh[index], 2) +
                        std::pow(m_mesh.yFourierMesh[index], 2) +
                        std::pow(m_mesh.zFourierMesh[index], 2);
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

std::vector<double> Grid3D::wavenumber() const { return m_mesh.wavenumber; }
