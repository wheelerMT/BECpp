#include <gtest/gtest.h>
#include <cmath>
#include "grid.h"
#include "constants.h"

class Grid1DTest : public ::testing::Test
{
public:
    Grid1D grid{128, 0.5};
};

TEST_F(Grid1DTest, PointsSetCorrectly)
{
    ASSERT_EQ(grid.xPoints(), 128);
}

TEST_F(Grid1DTest, SpacingSetCorrectly)
{
    ASSERT_EQ(grid.xGridSpacing(), 0.5);
}

TEST_F(Grid1DTest, FourierSpacingSetCorrectly)
{
    ASSERT_EQ(grid.xFourierGridSpacing(), PI / 32.0);
}

TEST_F(Grid1DTest, LengthSetCorrectly)
{
    ASSERT_EQ(grid.xLength(), 64);
}

TEST_F(Grid1DTest, WavenumberSetCorrectly)
{
    ASSERT_EQ(grid.wavenumber()[0], 0.0);
}

class Grid2DTest : public ::testing::Test
{
public:
    Grid2D grid{128, 128, 0.5, 0.5};
};

TEST_F(Grid2DTest, PointsSetCorrectly)
{
    ASSERT_EQ(grid.xPoints(), 128);
    ASSERT_EQ(grid.yPoints(), 128);
}

TEST_F(Grid2DTest, SpacingSetCorrectly)
{
    ASSERT_EQ(grid.xGridSpacing(), 0.5);
    ASSERT_EQ(grid.yGridSpacing(), 0.5);
}

TEST_F(Grid2DTest, FourierSpacingSetCorrectly)
{
    ASSERT_EQ(grid.xFourierGridSpacing(), PI / 32.0);
    ASSERT_EQ(grid.yFourierGridSpacing(), PI / 32.0);
}

TEST_F(Grid2DTest, LengthSetCorrectly)
{
    ASSERT_EQ(grid.xLength(), 64);
    ASSERT_EQ(grid.yLength(), 64);
}

TEST_F(Grid2DTest, WavenumberSetCorrectly)
{
    ASSERT_EQ(grid.wavenumber()[0][0], 0.0);
}
