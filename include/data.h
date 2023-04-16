#ifndef BECPP_DATA_H
#define BECPP_DATA_H

#include <string>
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"
#include "highfive/H5File.hpp"
#include "grid.h"
#include "wavefunction.h"

struct Parameters
{
    double intStrength{}; // Interaction strength
    std::vector<double> trap{};
    int numTimeSteps{};
    std::complex<double> timeStep{};  // Time step
    double currentTime{};
};

class DataManager1D
{
private:
    unsigned int m_saveIndex{0};

    void saveParameters(const Parameters &params, const Grid1D &grid);
    void generateWavefunctionDatasets(const Grid1D &grid);

public:
    DataManager1D(const std::string &filename, const Parameters &params, const Grid1D &grid);

    void saveWavefunctionData(Wavefunction1D &psi);

    std::string filename;
    HighFive::File file;
};

class DataManager2D
{
private:
    unsigned int m_saveIndex{0};

    void saveParameters(const Parameters &params, const Grid2D &grid);
    void generateWavefunctionDatasets(const Grid2D &grid);

public:
    DataManager2D(const std::string &filename, const Parameters &params, const Grid2D &grid);
    void saveWavefunctionData(Wavefunction2D &psi);

    std::string filename;
    HighFive::File file;
};

#endif //BECPP_DATA_H
