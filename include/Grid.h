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
    // Grid points
    std::pair<int, int> m_points{};
    int m_nx{};
    int m_ny{};

    // Grid spacing
    std::pair<double, double> m_grid_spacing{};
    double m_dx{};
    double m_dy{};
    double m_dkx{};
    double m_dky{};

    // Grid lengths
    double m_len_x{};
    double m_len_y{};

    // Grid & parameter construction functions
    void constructGridParams();
    void constructGrids();
    void fftshift();

public:
    // Grids
    std::vector<std::vector<double>> X{};
    std::vector<std::vector<double>> Y{};
    std::vector<std::vector<double>> Kx{};
    std::vector<std::vector<double>> Ky{};

    // Constructor
    Grid(std::pair<int, int> points, std::pair<double, double> grid_spacing);

    // Declare wavefunction class a friend
    friend class Wavefunction;
};
#endif //BECPP_GRID_H
