#include "grid.h"
#include "constants.h"
#include <algorithm>
#include <cmath>
#include <utility>

Grid1D::Grid1D(unsigned long nx, double dx)
    : m_gridPoints{nx}, m_gridSpacing{dx}
{
    Grid1D::constructGridParams();
    Grid1D::constructMesh();
}

void Grid1D::constructGridParams()
{
    m_fourierGridSpacing = PI / (m_gridPoints / 2. * m_gridSpacing);
    m_length = m_gridPoints * m_gridSpacing;
}

void resizeMesh1D(Mesh1D& mesh, unsigned long xPoints)
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

unsigned long Grid1D::shape() const { return m_gridPoints; }

double Grid1D::gridSpacing() const { return m_gridSpacing; }

double Grid1D::fourierGridSpacing() const { return m_fourierGridSpacing; }

double Grid1D::gridLength() const { return m_length; }

std::vector<double> Grid1D::wavenumber() const { return m_mesh.wavenumber; }

Grid2D::Grid2D(std::tuple<unsigned long, unsigned long> points,
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

void resizeMesh2D(Mesh2D& mesh, unsigned long xPoints, unsigned long yPoints)
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
            m_mesh.xMesh[j + i * yPoints] = (j - xPoints / 2.) * xGridSpacing;
            m_mesh.yMesh[j + i * yPoints] = (j - yPoints / 2.) * yGridSpacing;
            if (i < xPoints / 2)
            {
                m_mesh.xFourierMesh[j + i * yPoints] = i * xFourierGridSpacing;
                m_mesh.yFourierMesh[j + i * yPoints] = j * yFourierGridSpacing;
            } else
            {
                m_mesh.xFourierMesh[j + i * yPoints] =
                        (i - xPoints) * xFourierGridSpacing;
                m_mesh.yFourierMesh[j + i * yPoints] =
                        (j - yPoints) * yFourierGridSpacing;
            }

            m_mesh.wavenumber[j + i * yPoints] =
                    std::pow(m_mesh.xFourierMesh[j + i * yPoints], 2) +
                    std::pow(m_mesh.yFourierMesh[j + i * yPoints], 2);
        }
    }
}

std::tuple<unsigned long, unsigned long> Grid2D::shape() const
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

Grid3D::Grid3D(std::tuple<unsigned long, unsigned long, unsigned long> points,
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

void resizeMesh3D(Mesh3D& mesh, unsigned long xPoints, unsigned long yPoints,
                  unsigned long zPoints)
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
                m_mesh.xMesh[k + zPoints * (j + yPoints * i)] =
                        (i - xPoints / 2.) * xGridSpacing;
                m_mesh.yMesh[k + zPoints * (j + yPoints * i)] =
                        (j - yPoints / 2.) * yGridSpacing;
                m_mesh.zMesh[k + zPoints * (j + yPoints * i)] =
                        (k - zPoints / 2.) * zGridSpacing;

                if (i < xPoints / 2)
                {
                    m_mesh.xFourierMesh[k + zPoints * (j + yPoints * i)] =
                            i * xFourierGridSpacing;
                    m_mesh.yFourierMesh[k + zPoints * (j + yPoints * i)] =
                            j * yFourierGridSpacing;
                    m_mesh.zFourierMesh[k + zPoints * (j + yPoints * i)] =
                            k * yFourierGridSpacing;
                } else
                {
                    m_mesh.xFourierMesh[k + zPoints * (j + yPoints * i)] =
                            (i - xPoints) * xFourierGridSpacing;
                    m_mesh.yFourierMesh[k + zPoints * (j + yPoints * i)] =
                            (j - yPoints) * yFourierGridSpacing;
                    m_mesh.zFourierMesh[k + zPoints * (j + yPoints * i)] =
                            (k - zPoints) * yFourierGridSpacing;
                }

                m_mesh.wavenumber[k + zPoints * (j + yPoints * i)] =
                        std::pow(m_mesh.xFourierMesh[k + zPoints * (j + yPoints * i)],
                                 2) +
                        std::pow(m_mesh.yFourierMesh[k + zPoints * (j + yPoints * i)],
                                 2) +
                        std::pow(m_mesh.zFourierMesh[k + zPoints * (j + yPoints * i)],
                                 2);
            }
        }
    }
}

std::tuple<unsigned long, unsigned long, unsigned long> Grid3D::shape() const
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
