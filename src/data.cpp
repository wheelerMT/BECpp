//
// Created by mattw on 06/01/2022.
//

#include "data.h"

Parameters::Parameters(double c0, double c2, int nt, int nframe, double dt, doubleArray_t &V)
        : c0{c0}, c2{c2}, nt{nt}, nframe{nframe}, dt{dt}, V{V}
{
}

DataManager::DataManager(const std::string &filename, const Parameters &params, const Grid &grid) :
        filename{filename},
        file{filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate}
{
    save_parameters(params);
    generate_wfn_datasets(grid);

}

void DataManager::save_parameters(const Parameters &params)
{
    // Save condensate and time parameters to file
    file.createDataSet("/parameters/c0", params.c0);
    file.createDataSet("/parameters/c2", params.c2);
    file.createDataSet("/time/nt", params.nt);
    file.createDataSet("/time/nframe", params.nframe);
    file.createDataSet("/time/dt", params.dt);
}

void DataManager::generate_wfn_datasets(const Grid &grid)
{
    // Define dataspaces with arbitrary length of last dimension
    HighFive::DataSpace ds_plus = HighFive::DataSpace({grid.nx * grid.ny, 1},
                                                      {grid.nx * grid.ny, HighFive::DataSpace::UNLIMITED});
    HighFive::DataSpace ds_zero = HighFive::DataSpace({grid.nx * grid.ny, 1},
                                                      {grid.nx * grid.ny, HighFive::DataSpace::UNLIMITED});
    HighFive::DataSpace ds_minus = HighFive::DataSpace({grid.nx * grid.ny, 1},
                                                       {grid.nx * grid.ny, HighFive::DataSpace::UNLIMITED});

    // Use chunking
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{2, 1}));

    // Create wavefunction datasets
    file.createDataSet("/wavefunction/psi_plus", ds_plus,
                       HighFive::AtomicType<std::complex<double>>(),props);
    file.createDataSet("/wavefunction/psi_zero", ds_zero,
                       HighFive::AtomicType<std::complex<double>>(),props);
    file.createDataSet("/wavefunction/psi_minus", ds_minus,
                       HighFive::AtomicType<std::complex<double>>(),props);
}
