//
// Created by mattw on 06/01/2022.
//

#ifndef BECPP_DATA_H
#define BECPP_DATA_H

#include <string>
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"
#include "highfive/H5File.hpp"
#include "grid.h"
#include "wavefunction.h"

// Define arrays
using doubleArray_t = std::vector<std::vector<double>>;

struct Parameters
{
    Parameters(double c0, double c2, int nt, int nframe, double dt, doubleArray_t &V);

    double c0;
    double c2;
    int nt;
    int nframe;
    double dt;
    doubleArray_t V;

};

class DataManager
{
private:

    void save_parameters(const Parameters &params, const Grid &grid);

    void generate_wfn_datasets(const Grid &grid);

public:
    // Constructor
    DataManager(const std::string &filename, const Parameters &params, const Grid &grid);

    void save_wavefunction_data(Wavefunction &psi);

    std::string filename;
    HighFive::File file;
};


#endif //BECPP_DATA_H
