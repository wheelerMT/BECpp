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
    int nx{};
    int ny{};

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

    // Constructor
    Grid(int t_nx, int t_ny, double t_dx, double t_dy);

    // Declare wavefunction class a friend
    friend class Wavefunction;
};

#endif //BECPP_GRID_H
