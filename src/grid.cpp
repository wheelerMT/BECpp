//
// Created by mattw on 04/01/2022.
//

#include <cmath>
#include <algorithm>
#include "grid.h"
#include "constants.h"

void Grid2D::construct_grid_params()
{
    // Calculate k-space grid spacing
    dkx = PI / (nx / 2. * dx);
    dky = PI / (ny / 2. * dy);

    // Sets length of sides of box
    len_x = nx * dx;
    len_y = ny * dy;
}

void Grid2D::construct_grids()
{
    // Set grid sizes
    X.resize(nx, std::vector<double>(ny));
    Y.resize(nx, std::vector<double>(ny));
    Kx.resize(nx, std::vector<double>(ny));
    Ky.resize(nx, std::vector<double>(ny));
    K.resize(nx, std::vector<double>(ny));

    // Construct grids
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            X[i][j] = (j - nx / 2.) * dx;
            Kx[i][j] = (j - nx / 2.) * dkx;
            Y[j][i] = (j - ny / 2.) * dy;
            Ky[j][i] = (j - ny / 2.) * dky;
            K[i][j] = std::pow(Kx[i][j], 2) + std::pow(Ky[i][j], 2);
        }
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

void Grid2D::fftshift()
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

    // Update wavevector
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            K[i][j] = std::pow(Kx[i][j], 2) + std::pow(Ky[i][j], 2);
        }
    }
}

Grid2D::Grid2D(unsigned int nx, unsigned int ny, double dx, double dy)
        : nx{nx}, ny{ny}, dx{dx}, dy{dy}
{
    construct_grid_params();
    construct_grids();
}

Grid2D::Grid2D(const Grid2D &grid) : nx{grid.nx}, ny{grid.ny}, dx{grid.dx}, dy{grid.dy}
{
    construct_grid_params();
    construct_grids();
}

void Grid1D::construct_grid_params()
{
    // Calculate k-space grid spacing
    dkx = PI / (nx / 2. * dx);

    // Sets length of side of box
    len_x = nx * dx;
}

void Grid1D::construct_grids()
{
    // Set grid sizes
    X.resize(nx);
    Kx.resize(nx);
    K.resize(nx);

    // Construct grids
    for (int i = 0; i < nx; ++i)
    {
        X[i] = (i - nx / 2.) * dx;
        Kx[i] = (i - nx / 2.) * dkx;
        K[i] = std::pow(Kx[i], 2);
    }

    // Shift the k-space grids, so they are in the right order
    fftshift();
}

Grid1D::Grid1D(unsigned int nx, double dx) : nx{nx}, dx{dx}
{
    construct_grid_params();
    construct_grids();
}

Grid1D::Grid1D(const Grid1D &grid) : nx{grid.nx}, dx{grid.dx}
{
    construct_grid_params();
    construct_grids();
}

void Grid1D::fftshift()
{
    std::rotate(Kx.begin(), Kx.begin() + nx / 2, Kx.end());
    K = Kx;
}
