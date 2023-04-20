#include "BECpp.h"
#include "data.h"
#include "evolution.h"
#include "grid.h"
#include "wavefunction.h"

constexpr auto GRID_POINTS_X = 128;
constexpr auto GRID_POINTS_Y = 128;
constexpr auto GRID_SPACING_X = 0.5;
constexpr auto GRID_SPACING_Y = 0.5;

Parameters createParams()
{
    Parameters params{};
    params.intStrength = 1.0;
    params.trap = std::vector<double>(GRID_POINTS_X * GRID_POINTS_Y, 0.0);
    params.numTimeSteps = 250;
    params.timeStep = std::complex<double>{0.0, -1e-2};

    return params;
}

int main(int argc, char* argv[])
{

    // Create grid object
    std::tuple<unsigned int, unsigned int> points{GRID_POINTS_X, GRID_POINTS_Y};
    std::tuple<double, double> gridSpacing{GRID_SPACING_X, GRID_SPACING_Y};
    Grid2D grid{points, gridSpacing};

    // Create wavefunction and gaussian initial state
    complexVector_t initialState(GRID_POINTS_X * GRID_POINTS_Y);
    for (int i = 0; i < GRID_POINTS_X; ++i)
    {
        for (int j = 0; j < GRID_POINTS_Y; ++j)
        {
            auto index = j + i * GRID_POINTS_Y;
            initialState[index] = exp((-pow(grid.xMesh()[index], 2) -
                                       pow(grid.yMesh()[index], 2)) /
                                      500);
        }
    }
    Wavefunction2D wavefunction{grid};
    wavefunction.setComponent(initialState);

    // Create parameters
    Parameters params = createParams();

    // Create data manager
    DataManager2D dm{"groundState.h5", params, grid};

    // Evolution loop
    for (int i = 0; i < params.numTimeSteps; ++i)
    {
        std::cout << "On iteration " << i << "\n";
        fourierStep(wavefunction, params);
        wavefunction.ifft();
        interactionStep(wavefunction, params);
        wavefunction.fft();
        fourierStep(wavefunction, params);
        if (params.timeStep.imag() != 0.0) { renormaliseAtomNum(wavefunction); }

        // Save wavefunction data every 10 time steps
        if (i % 50 == 0)
        {
            wavefunction.ifft();
            dm.saveWavefunctionData(wavefunction);
        }
    }

    return EXIT_SUCCESS;
}
