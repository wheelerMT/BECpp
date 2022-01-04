//
// Created by mattw on 04/01/2022.
//


#include "Grid.h"
#include "constants.h"

void Grid::constructGridParams()
{
    m_nx = std::get<0>(m_points);
    m_ny = std::get<1>(m_points);

    // Gets and calculates grid spacings
    m_dx = std::get<0>(m_grid_spacing);
    m_dy = std::get<1>(m_grid_spacing);
    m_dkx = PI / (m_nx / 2. * m_dx);
    m_dky = PI / (m_ny / 2. * m_dy);

    // Sets length of sides of box
    m_len_x = m_nx * m_dx;
    m_len_y = m_ny * m_dy;

}

void Grid::constructGrids()
{
    // Set grid sizes
    X.resize(m_nx, std::vector<double>(m_ny));
    Y.resize(m_nx, std::vector<double>(m_ny));
    Kx.resize(m_nx, std::vector<double>(m_ny));
    Ky.resize(m_nx, std::vector<double>(m_ny));

    // Construct grids
    for (int i = 0; i < m_nx; ++i)
    {
        for (int j = 0; j < m_ny; ++j)
        {
            X[i][j] = (j - m_nx / 2.) * m_dx;
            Kx[i][j] = (j - m_nx / 2.) * m_dkx;
            Y[j][i] = (j - m_ny / 2.) * m_dy;
            Ky[j][i] = (j - m_ny / 2.) * m_dky;
        }
    }

    // Shift the k-space grids so they are in the right order
    fftshift();
}

void Grid::fftshift()
{
    /*
        Shifts the zero-frequency component to the center
        of the spectrum.
    */

    /*
        This is a poor implementation.
        It would be nice to make use of std::rotate in a future implementation.
    */

    std::vector<std::vector<double>> Kx_copy = Kx;
    std::vector<std::vector<double>> Ky_copy = Ky;

    // Reverse each row
    for (int i = 0; i < m_nx; i++)
    {
        for (int j = 0; j < m_ny; j++)
        {
            if (j < m_nx / 2)
            {
                Kx[i][j] = Kx_copy[i][m_nx / 2 + j];
                Ky[i][j] = Ky_copy[i][m_nx / 2 + j];
            }
            else if (j >= m_nx / 2)
            {
                Kx[i][j] = Kx_copy[i][j - m_nx / 2];
                Ky[i][j] = Ky_copy[i][j - m_nx / 2];
            }
        }
    }

    // Update copy arrays
    Kx_copy = Kx;
    Ky_copy = Ky;

    // Reverse each column
    for (int i = 0; i < m_nx; i++)
    {
        for (int j = 0; j < m_ny; j++)
        {
            if (j < m_nx / 2)
            {
                Kx[j][i] = Kx_copy[m_nx / 2 + j][i];
                Ky[j][i] = Ky_copy[m_nx / 2 + j][i];
            }
            else if (j >= m_nx / 2)
            {
                Kx[j][i] = Kx_copy[j - m_nx / 2][i];
                Ky[j][i] = Ky_copy[j - m_nx / 2][i];
            }
        }
    }
}

Grid::Grid(const std::pair<int, int> points, const std::pair<double, double> grid_spacing)
: m_points{points}, m_grid_spacing{grid_spacing}
{
    constructGridParams();
    constructGrids();
}