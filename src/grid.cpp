//
// Created by mattw on 04/01/2022.
//


#include "grid.h"
#include "constants.h"

void Grid::constructGridParams()
{
    // Calculate k-space grid spacing
    dkx = PI / (nx / 2. * dx);
    dky = PI / (ny / 2. * dy);

    // Sets length of sides of box
    len_x = nx * dx;
    len_y = ny * dy;
}

void Grid::constructGrids()
{
    // Set grid sizes
    X.resize(nx, std::vector<double>(ny));
    Y.resize(nx, std::vector<double>(ny));
    Kx.resize(nx, std::vector<double>(ny));
    Ky.resize(nx, std::vector<double>(ny));

    // Construct grids
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            X[i][j] = (j - nx / 2.) * dx;
            Kx[i][j] = (j - nx / 2.) * dkx;
            Y[j][i] = (j - ny / 2.) * dy;
            Ky[j][i] = (j - ny / 2.) * dky;
        }
    }

    // Shift the k-space grids, so they are in the right order
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
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (j < nx / 2)
            {
                Kx[i][j] = Kx_copy[i][nx / 2 + j];
                Ky[i][j] = Ky_copy[i][nx / 2 + j];
            } else if (j >= nx / 2)
            {
                Kx[i][j] = Kx_copy[i][j - nx / 2];
                Ky[i][j] = Ky_copy[i][j - nx / 2];
            }
        }
    }

    // Update copy arrays
    Kx_copy = Kx;
    Ky_copy = Ky;

    // Reverse each column
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            if (j < nx / 2)
            {
                Kx[j][i] = Kx_copy[nx / 2 + j][i];
                Ky[j][i] = Ky_copy[nx / 2 + j][i];
            } else if (j >= nx / 2)
            {
                Kx[j][i] = Kx_copy[j - nx / 2][i];
                Ky[j][i] = Ky_copy[j - nx / 2][i];
            }
        }
    }
}

Grid::Grid(int nx, int ny, double dx, double dy)
        : nx{nx}, ny{ny}, dx{dx}, dy{dy}
{
    constructGridParams();
    constructGrids();
}

Grid::Grid(const Grid &grid) : nx{grid.nx}, ny{grid.ny}, dx{grid.dx}, dy{grid.dy}
{
    constructGridParams();
    constructGrids();
}
