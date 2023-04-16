#include "data.h"

DataManager1D::DataManager1D(const std::string& filename,
                             const Parameters& params, const Grid1D& grid)
    : filename{filename},
      file{filename, HighFive::File::ReadWrite | HighFive::File::Create |
                             HighFive::File::Truncate}
{
    saveParameters(params, grid);
    generateWavefunctionDatasets(grid);
}

void DataManager1D::saveParameters(const Parameters& params, const Grid1D& grid)
{
    // Save condensate and time parameters to file
    file.createDataSet("/parameters/intStrength", params.intStrength);
    file.createDataSet("/parameters/numTimeSteps", params.numTimeSteps);
    file.createDataSet("/parameters/dt", params.timeStep);

    // Save grid parameters to file
    file.createDataSet("/grid/xPoints", grid.shape());
    file.createDataSet("/grid/xGridSpacing", grid.gridSpacing());
}

void DataManager1D::generateWavefunctionDatasets(const Grid1D& grid)
{
    // Define data space with arbitrary length of last dimension
    HighFive::DataSpace dsWavefunction = HighFive::DataSpace(
            {grid.shape(), 1}, {grid.shape(), HighFive::DataSpace::UNLIMITED});

    // Use chunking
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{(grid.shape()) / 4, 1}));

    // Create wavefunction dataset
    file.createDataSet("wavefunction", dsWavefunction,
                       HighFive::AtomicType<std::complex<double>>(), props);
}

void DataManager1D::saveWavefunctionData(Wavefunction1D& wfn)
{
    // Load in datasets
    HighFive::DataSet dsWavefunction = file.getDataSet("wavefunction");

    // Resize datasets
    dsWavefunction.resize({wfn.grid().shape(), saveIndex + 1});

    // FFT so we update real-space arrays
    wfn.ifft();

    // Save new wavefunction data
    dsWavefunction.select({0, saveIndex}, {wfn.grid().shape(), 1})
            .write(wfn.component());

    saveIndex += 1;
}
