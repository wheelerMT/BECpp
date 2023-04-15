#include "constants.h"
#include "grid.h"
#include <cmath>
#include <gtest/gtest.h>

class Grid1DTest : public ::testing::Test
{
public:
    Grid1D grid{128, 0.5};
};

TEST_F(Grid1DTest, PointsSetCorrectly) { ASSERT_EQ(grid.shape(), 128); }

TEST_F(Grid1DTest, SpacingSetCorrectly) { ASSERT_EQ(grid.gridSpacing(), 0.5); }

TEST_F(Grid1DTest, FourierSpacingSetCorrectly)
{
    ASSERT_EQ(grid.fourierGridSpacing(), PI / 32.0);
}

TEST_F(Grid1DTest, LengthSetCorrectly) { ASSERT_EQ(grid.gridLength(), 64); }

TEST_F(Grid1DTest, WavenumberSetCorrectly)
{
    ASSERT_EQ(grid.wavenumber()[0], 0.0);
}

class Grid2DTest : public ::testing::Test
{
public:
    std::tuple<unsigned int, unsigned int> points{128, 128};
    std::tuple<double, double> gridSpacing{0.5, 0.5};
    Grid2D grid{points, gridSpacing};
};

TEST_F(Grid2DTest, PointsSetCorrectly)
{
    auto [xPoints, yPoints] = grid.shape();
    ASSERT_EQ(xPoints, 128);
    ASSERT_EQ(yPoints, 128);
}

TEST_F(Grid2DTest, SpacingSetCorrectly)
{
    auto [xGridSpacing, yGridSpacing] = grid.gridSpacing();
    ASSERT_EQ(xGridSpacing, 0.5);
    ASSERT_EQ(yGridSpacing, 0.5);
}

TEST_F(Grid2DTest, FourierSpacingSetCorrectly)
{
    auto [xFourierGridSpacing, yFourierGridSpacing] = grid.fourierGridSpacing();
    ASSERT_EQ(xFourierGridSpacing, PI / 32.0);
    ASSERT_EQ(yFourierGridSpacing, PI / 32.0);
}

TEST_F(Grid2DTest, LengthSetCorrectly)
{
    auto [xLength, yLength] = grid.gridLength();
    ASSERT_EQ(xLength, 64);
    ASSERT_EQ(yLength, 64);
}

TEST_F(Grid2DTest, WavenumberSetCorrectly)
{
    ASSERT_EQ(grid.wavenumber()[0], 0.0);
}

class Grid3DTest : public ::testing::Test
{
public:
    std::tuple<unsigned int, unsigned int, unsigned int> points{64, 64, 64};
    std::tuple<double, double, double> gridSpacings{0.5, 0.5, 0.5};
    Grid3D grid{points, gridSpacings};
};

TEST_F(Grid3DTest, PointsSetCorrectly)
{
    auto [xPoints, yPoints, zPoints] = grid.shape();
    ASSERT_EQ(xPoints, 64);
    ASSERT_EQ(yPoints, 64);
    ASSERT_EQ(zPoints, 64);
}

TEST_F(Grid3DTest, SpacingSetCorrectly)
{
    auto [xGridSpacing, yGridSpacing, zGridSpacing] = grid.gridSpacing();
    ASSERT_EQ(xGridSpacing, 0.5);
    ASSERT_EQ(yGridSpacing, 0.5);
    ASSERT_EQ(zGridSpacing, 0.5);
}

TEST_F(Grid3DTest, FourierSpacingSetCorrectly)
{
    auto [xFourierGridSpacing, yFourierGridSpacing, zFourierGridSpacing] =
            grid.fourierGridSpacing();
    ASSERT_EQ(xFourierGridSpacing, PI / 16.0);
    ASSERT_EQ(yFourierGridSpacing, PI / 16.0);
    ASSERT_EQ(zFourierGridSpacing, PI / 16.0);
}

TEST_F(Grid3DTest, LengthSetCorrectly)
{
    auto [xLength, yLength, zLength] = grid.gridLength();
    ASSERT_EQ(xLength, 32);
    ASSERT_EQ(yLength, 32);
    ASSERT_EQ(zLength, 32);
}

TEST_F(Grid3DTest, WavenumberSetCorrectly)
{
    ASSERT_EQ(grid.wavenumber()[0], 0.0);
}