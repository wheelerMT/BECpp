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
    unsigned int saveIndex{0};

    void saveParameters(const Parameters &params, const Grid1D &grid);
    void generateWavefunctionDatasets(const Grid1D &grid);

public:
    DataManager1D(const std::string &filename, const Parameters &params, const Grid1D &grid);

    void saveWavefunctionData(Wavefunction1D &psi);

    std::string filename;
    HighFive::File file;
};


#endif //BECPP_DATA_H
