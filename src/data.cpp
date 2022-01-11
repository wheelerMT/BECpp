//
// Created by mattw on 06/01/2022.
//

#include "data.h"

Parameters::Parameters(double c0, double c2, double q, int nt, int nframe, double dt, doubleArray_t &V)
        : c0{c0}, c2{c2}, q{q}, nt{nt}, nframe{nframe}, dt{dt}, V{V}
{
}

DataManager::DataManager(const std::string &filename, const Parameters &params, const Grid &grid) :
        filename{filename},
        file{filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate}
{
    save_parameters(params, grid);
    generate_wfn_datasets(grid);

}

void DataManager::save_parameters(const Parameters &params, const Grid &grid)
{
    // Save condensate and time parameters to file
    file.createDataSet("/parameters/c0", params.c0);
    file.createDataSet("/parameters/c2", params.c2);
    file.createDataSet("/time/nt", params.nt);
    file.createDataSet("/time/nframe", params.nframe);
    file.createDataSet("/time/dt", params.dt);

    // Save grid parameters to file
    file.createDataSet("/grid/nx", grid.nx);
    file.createDataSet("/grid/ny", grid.ny);
    file.createDataSet("/grid/dx", grid.dx);
    file.createDataSet("/grid/dy", grid.dy);

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
    props.add(HighFive::Chunking(std::vector<hsize_t>{(grid.nx * grid.ny) / 4, 1}));

    // Create wavefunction datasets
    file.createDataSet("/wavefunction/psi_plus", ds_plus,
                       HighFive::AtomicType<std::complex<double>>(), props);
    file.createDataSet("/wavefunction/psi_zero", ds_zero,
                       HighFive::AtomicType<std::complex<double>>(), props);
    file.createDataSet("/wavefunction/psi_minus", ds_minus,
                       HighFive::AtomicType<std::complex<double>>(), props);

}

void DataManager::save_wavefunction_data(Wavefunction &psi)
{
    // Load in datasets
    HighFive::DataSet ds_plus = file.getDataSet("/wavefunction/psi_plus");
    HighFive::DataSet ds_zero = file.getDataSet("/wavefunction/psi_zero");
    HighFive::DataSet ds_minus = file.getDataSet("/wavefunction/psi_minus");

    // Resize datasets
    ds_plus.resize({psi.grid.nx * psi.grid.ny, save_index + 1});
    ds_zero.resize({psi.grid.nx * psi.grid.ny, save_index + 1});
    ds_minus.resize({psi.grid.nx * psi.grid.ny, save_index + 1});

    // Save new wavefunction data
    ds_plus.select({0, save_index}, {psi.grid.nx * psi.grid.ny, 1}).write(psi.plus);
    ds_zero.select({0, save_index}, {psi.grid.nx * psi.grid.ny, 1}).write(psi.zero);
    ds_minus.select({0, save_index}, {psi.grid.nx * psi.grid.ny, 1}).write(psi.minus);

    // Increase save index
    save_index += 1;
}
