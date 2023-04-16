#include <gtest/gtest.h>
#include "data.h"
#include "grid.h"

constexpr auto GRID_LENGTH = 32;
constexpr auto GRID_SPACING = 0.5;

Parameters parameters()
{
    Parameters params{};
    params.intStrength = 1.0;
    params.numTimeSteps = 100;
    params.timeStep = std::complex<double> {1e-2, 0};
    params.currentTime = 0.0;

    return params;
}

class DataManager1DTest : public ::testing::Test
{
public:
    Grid1D grid{GRID_LENGTH, GRID_SPACING};
    Parameters params = parameters();
    DataManager1D dm{"test_file.h5", params, grid};
};

TEST_F(DataManager1DTest, TestDataSetsCorrect)
{
    ASSERT_TRUE(dm.file.exist("parameters/intStrength"));
    ASSERT_TRUE(dm.file.exist("parameters/numTimeSteps"));
    ASSERT_TRUE(dm.file.exist("parameters/dt"));
    ASSERT_TRUE(dm.file.exist("grid/xPoints"));
    ASSERT_TRUE(dm.file.exist("grid/xGridSpacing"));
    ASSERT_TRUE(dm.file.exist("wavefunction"));
}

TEST_F(DataManager1DTest, TestWavefunctionSaved)
{
    Wavefunction1D wfn{grid};
    complexVector_t initialState{};
    initialState.resize(GRID_LENGTH, std::complex<double>{1.0, 0.0});
    wfn.setComponent(initialState);
    dm.saveWavefunctionData(wfn);

    auto wfnDataSet = dm.file.getDataSet("wavefunction");
    auto loadedWfn = wfnDataSet.read<std::vector<std::complex<double>>>();
    wfnDataSet.read(loadedWfn);

    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        ASSERT_EQ(loadedWfn[i], wfn.component()[i]);
    }
}