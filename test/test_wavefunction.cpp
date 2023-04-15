#include "wavefunction.h"
#include <gtest/gtest.h>

constexpr auto GRID_LENGTH = 16;
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

    for (int i = 0; i < GRID_LENGTH; i++)
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

class Wavefunction2DTest : public ::testing::Test
{
public:
    std::tuple<unsigned long, unsigned long> points{GRID_LENGTH, GRID_LENGTH};
    std::tuple<double, double> spacing{GRID_SPACING, GRID_SPACING};
    Grid2D grid{points, spacing};
    Wavefunction2D wavefunction{grid};
};

TEST_F(Wavefunction2DTest, SetComponentCorrect)
{
    std::vector<std::complex<double>> initialState{};
    std::complex<double> real{1.0, 0.0};
    initialState.resize(GRID_LENGTH * GRID_LENGTH, real);

    wavefunction.setComponent(initialState);
    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        for (int j = 0; j < GRID_LENGTH; ++j)
        {
            auto index = j + i * GRID_LENGTH;
            ASSERT_EQ(wavefunction.component()[index], real);
        }
    }
}

TEST_F(Wavefunction2DTest, AtomNumberCorrect)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    double atomNumber = GRID_LENGTH * GRID_LENGTH * GRID_SPACING * GRID_SPACING;
    ASSERT_EQ(atomNumber, wavefunction.atomNumber());
}

TEST_F(Wavefunction2DTest, DensityCorrect)
{
    complexVector_t initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        for (int j = 0; j < GRID_LENGTH; ++j)
        {
            auto index = j + i * GRID_LENGTH;
            ASSERT_EQ(wavefunction.density()[index], 1.0);
        }
    }
}

TEST_F(Wavefunction2DTest, SameWavefunctionAfterFFT)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    wavefunction.fft();
    wavefunction.ifft();

    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        for (int j = 0; j < GRID_LENGTH; ++j)
        {
            auto index = j + i * GRID_LENGTH;
            ASSERT_EQ(initialState[index], wavefunction.component()[index]);
        }
    }
}

TEST_F(Wavefunction2DTest, FourierSpaceWavefunctionUpdated)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH, {1.0, 2.0});
    wavefunction.setComponent(initialState);

    std::complex<double> zero{0.0, 0.0};

    // Sufficient to check one element non-zero
    wavefunction.fft();
    ASSERT_NE(zero, wavefunction.fourierComponent()[0]);
}

class Wavefunction3DTest : public ::testing::Test
{
public:
    std::tuple<unsigned long, unsigned long, unsigned long> points{
            GRID_LENGTH, GRID_LENGTH, GRID_LENGTH};
    std::tuple<double, double, double> spacing{GRID_SPACING, GRID_SPACING,
                                               GRID_SPACING};
    Grid3D grid{points, spacing};
    Wavefunction3D wavefunction{grid};
};

TEST_F(Wavefunction3DTest, SetComponentCorrect)
{
    std::vector<std::complex<double>> initialState{};
    std::complex<double> real{1.0, 0.0};
    initialState.resize(GRID_LENGTH * GRID_LENGTH * GRID_LENGTH, real);

    wavefunction.setComponent(initialState);
    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        for (int j = 0; j < GRID_LENGTH; ++j)
        {
            for (int k = 0; k < GRID_LENGTH; ++k)
            {
                auto index = k + GRID_LENGTH * (j + i * GRID_LENGTH);
                ASSERT_EQ(wavefunction.component()[index], real);
            }
        }
    }
}

TEST_F(Wavefunction3DTest, AtomNumberCorrect)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH * GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    double atomNumber = std::pow(GRID_LENGTH, 3) * std::pow(GRID_SPACING, 3);
    ASSERT_EQ(atomNumber, wavefunction.atomNumber());
}

TEST_F(Wavefunction3DTest, DensityCorrect)
{
    complexVector_t initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH * GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        for (int j = 0; j < GRID_LENGTH; ++j)
        {
            for (int k = 0; k < GRID_LENGTH; ++k)
            {
                auto index = k + GRID_LENGTH * (j + i * GRID_LENGTH);
                ASSERT_EQ(wavefunction.density()[index], 1.0);
            }
        }
    }
}

TEST_F(Wavefunction3DTest, SameWavefunctionAfterFFT)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH * GRID_LENGTH, {1.0, 0.0});
    wavefunction.setComponent(initialState);

    wavefunction.fft();
    wavefunction.ifft();

    for (int i = 0; i < GRID_LENGTH; ++i)
    {
        for (int j = 0; j < GRID_LENGTH; ++j)
        {
            for (int k = 0; k < GRID_LENGTH; ++k)
            {
                auto index = k + GRID_LENGTH * (j + i * GRID_LENGTH);
                ASSERT_EQ(initialState[index], wavefunction.component()[index]);
            }
        }
    }
}

TEST_F(Wavefunction3DTest, FourierSpaceWavefunctionUpdated)
{
    std::vector<std::complex<double>> initialState{};
    initialState.resize(GRID_LENGTH * GRID_LENGTH * GRID_LENGTH, {1.0, 2.0});
    wavefunction.setComponent(initialState);

    std::complex<double> zero{0.0, 0.0};

    // Sufficient to check one element non-zero
    wavefunction.fft();
    ASSERT_NE(zero, wavefunction.fourierComponent()[0]);
}