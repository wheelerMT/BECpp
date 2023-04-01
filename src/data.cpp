//
// Created by mattw on 06/01/2022.
//

#include "data.h"

Parameters::Parameters(double c0, double c2, double p, double q, int nt, int nframe, double dt)
        : c0{c0}, c2{c2}, p{p}, q{q}, nt{nt}, nframe{nframe}, dt{dt, 0}
{
}

void Parameters::imaginary_time(const std::string &toggle)
{
    assert(toggle == "on" or toggle == "off"); // Check correct toggle string is passed
    if (toggle == "on") dt *= std::complex<double>{0, -1};
    else if (toggle == "off") dt /= std::complex<double>{0, -1};
}

DataManager::DataManager(const std::string &filename, const Parameters &params, const Grid2D &grid) :
        filename{filename},
        file{filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate}
{
    save_parameters(params, grid);
    generate_wfn_datasets(grid);

}

void DataManager::save_parameters(const Parameters &params, const Grid2D &grid)
{
    // Save condensate and time parameters to file
    file.createDataSet("/parameters/c0", params.c0);
    file.createDataSet("/parameters/c2", params.c2);
    file.createDataSet("/time/nt", params.nt);
    file.createDataSet("/time/nframe", params.nframe);
    file.createDataSet("/time/dt", params.dt.real());

    // Save grid parameters to file
    file.createDataSet("/grid/m_xPoints", grid.m_xPoints);
    file.createDataSet("/grid/m_yPoints", grid.m_yPoints);
    file.createDataSet("/grid/m_xGridSpacing", grid.m_xGridSpacing);
    file.createDataSet("/grid/m_yGridSpacing", grid.m_yGridSpacing);

}

void DataManager::generate_wfn_datasets(const Grid2D &grid)
{
    // Define dataspaces with arbitrary length of last dimension
    HighFive::DataSpace ds_plus = HighFive::DataSpace({grid.m_xPoints * grid.m_yPoints, 1},
                                                      {grid.m_xPoints * grid.m_yPoints, HighFive::DataSpace::UNLIMITED});
    HighFive::DataSpace ds_zero = HighFive::DataSpace({grid.m_xPoints * grid.m_yPoints, 1},
                                                      {grid.m_xPoints * grid.m_yPoints, HighFive::DataSpace::UNLIMITED});
    HighFive::DataSpace ds_minus = HighFive::DataSpace({grid.m_xPoints * grid.m_yPoints, 1},
                                                       {grid.m_xPoints * grid.m_yPoints, HighFive::DataSpace::UNLIMITED});

    // Use chunking
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{(grid.m_xPoints * grid.m_yPoints) / 4, 1}));

    // Create wavefunction datasets
    file.createDataSet("/wavefunction/psi_plus", ds_plus,
                       HighFive::AtomicType<std::complex<double>>(), props);
    file.createDataSet("/wavefunction/psi_zero", ds_zero,
                       HighFive::AtomicType<std::complex<double>>(), props);
    file.createDataSet("/wavefunction/psi_minus", ds_minus,
                       HighFive::AtomicType<std::complex<double>>(), props);

}

void DataManager::save_wavefunction_data(Wavefunction2D &psi)
{
    // Load in datasets
    HighFive::DataSet ds_plus = file.getDataSet("/wavefunction/psi_plus");
    HighFive::DataSet ds_zero = file.getDataSet("/wavefunction/psi_zero");
    HighFive::DataSet ds_minus = file.getDataSet("/wavefunction/psi_minus");

    // Resize datasets
    ds_plus.resize({psi.grid.m_xPoints * psi.grid.m_yPoints, save_index + 1});
    ds_zero.resize({psi.grid.m_xPoints * psi.grid.m_yPoints, save_index + 1});
    ds_minus.resize({psi.grid.m_xPoints * psi.grid.m_yPoints, save_index + 1});

    // FFT so we update real-space arrays
    psi.ifft();

    // Save new wavefunction data
    ds_plus.select({0, save_index}, {psi.grid.m_xPoints * psi.grid.m_yPoints, 1}).write(psi.plus);
    ds_zero.select({0, save_index}, {psi.grid.m_xPoints * psi.grid.m_yPoints, 1}).write(psi.zero);
    ds_minus.select({0, save_index}, {psi.grid.m_xPoints * psi.grid.m_yPoints, 1}).write(psi.minus);

    // Increase save index
    save_index += 1;
}
