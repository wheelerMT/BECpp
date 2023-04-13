#include "wavefunction.h"
#include <gtest/gtest.h>

constexpr auto GRID_LENGTH = 128;
constexpr auto GRID_SPACING = 0.5;

class Wavefunction1DTest : public ::testing::Test
{
public:
    Grid1D grid{GRID_LENGTH, GRID_SPACING};
    Wavefunction1D wavefunction{grid};
};

TEST_F(Wavefunction1DTest, SetComponentCorrect)
{
    std::vector<std::complex<double>> initialState{};
    std::complex<double> real{1.0, 0.0};
    initialState.resize(GRID_LENGTH, real);

    wavefunction.setComponent(initialState);
    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        ASSERT_EQ(wavefunction.component()[i], real);
    }
}

TEST_F(Wavefunction1DTest, AtomNumberCorrect)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    ASSERT_EQ(GRID_LENGTH * GRID_SPACING, wavefunction.atomNumber());
}

TEST_F(Wavefunction1DTest, DensityCorrect)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    for (size_t i = 0; i < GRID_LENGTH; i++)
    {
        ASSERT_EQ(wavefunction.density()[i], 1.0);
    }
}

TEST_F(Wavefunction1DTest, SameWavefunctionAfterFFT)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    wavefunction.fft();
    wavefunction.ifft();

    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        ASSERT_EQ(initialState[i], wavefunction.component()[i]);
    }
}

TEST_F(Wavefunction1DTest, FourierSpaceWavefunctionUpdated)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH, {1.0, 2.0});
    wavefunction.setComponent(initialState);

    std::complex<double> zero{0.0, 0.0};

    // Sufficient to check one element non-zero
    wavefunction.fft();
    ASSERT_NE(zero, wavefunction.fourierComponent()[0]);
}
