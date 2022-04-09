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
    Parameters(double c0, double c2, double p, double q, int nt, int nframe, double dt);

    double c0;
    double c2;
    double p;
    double q;
    int nt;
    int nframe;
    std::complex<double> dt;

    void imaginary_time(const std::string &toggle);

};

class DataManager
{
private:
    unsigned int save_index{0};

    void save_parameters(const Parameters &params, const Grid2D &grid);

    void generate_wfn_datasets(const Grid2D &grid);

public:
    // Constructor
    DataManager(const std::string &filename, const Parameters &params, const Grid2D &grid);

    void save_wavefunction_data(Wavefunction2D &psi);

    std::string filename;
    HighFive::File file;
};


#endif //BECPP_DATA_H
