//
// Created by mattw on 04/01/2022.
//

#ifndef BECPP_GRID_H
#define BECPP_GRID_H

#include <utility>
#include <vector>

using doubleArray_t = std::vector<std::vector<double>>;

class Grid2D
{
private:
    // Grid2D & parameter construction functions
    void construct_grid_params();

    void construct_grids();

    void fftshift();

public:
    // Grid2D points
    const unsigned int nx{};
    const unsigned int ny{};

    // Grid2D spacing
    double dx{};
    double dy{};
    double dkx{};
    double dky{};

    // Grid2D lengths
    double len_x{};
    double len_y{};

    // Grids
    doubleArray_t X{};
    doubleArray_t Y{};
    doubleArray_t Kx{};
    doubleArray_t Ky{};
    doubleArray_t K{};

    // Constructors
    Grid2D(unsigned int nx, unsigned int ny, double dx, double dy);

    Grid2D(const Grid2D &grid);  // Copy constructor

    // Declare wavefunction class a friend
    friend class Wavefunction;
};

#endif //BECPP_GRID_H
