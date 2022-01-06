//
// Created by mattw on 04/01/2022.
//

#ifndef BECPP_GRID_H
#define BECPP_GRID_H

#include <utility>
#include <vector>

class Grid
{
private:
    // Grid & parameter construction functions
    void constructGridParams();

    void constructGrids();

    void fftshift();

public:
    // Grid points
    const unsigned int nx{};
    const unsigned int ny{};

    // Grid spacing
    double dx{};
    double dy{};
    double dkx{};
    double dky{};

    // Grid lengths
    double len_x{};
    double len_y{};

    // Grids
    std::vector<std::vector<double>> X{};
    std::vector<std::vector<double>> Y{};
    std::vector<std::vector<double>> Kx{};
    std::vector<std::vector<double>> Ky{};

    // Constructors
    Grid(unsigned int nx, unsigned int ny, double dx, double dy);

    Grid(const Grid &grid);  // Copy constructor

    // Declare wavefunction class a friend
    friend class Wavefunction;
};

#endif //BECPP_GRID_H
