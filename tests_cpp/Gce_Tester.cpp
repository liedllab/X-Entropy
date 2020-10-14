#include "pch.h"

#include "../src/kde.h"

TEST(GceTest, ConstructorTest)
{
  Gce gce{ {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 10 };
  ASSERT_EQ(gce.getHistogram().size(), 10);
  ASSERT_EQ(gce.getGridLength(), 11);
  ASSERT_EQ(gce.getGrid().size(), 11);
  EXPECT_DOUBLE_EQ( gce.getHistogram()[0], 1.0 / 11.0);
  EXPECT_DOUBLE_EQ( gce.getHistogram()[4], 2.0 / 11.0); // Middle is double
  EXPECT_DOUBLE_EQ(gce.getGrid()[0], -0.1);
  EXPECT_DOUBLE_EQ(gce.getGrid()[10], 1.1);
  EXPECT_THROW(Gce( std::vector<double>(), 10), EmptyListError);
}


TEST(GceTest, CalculateTest)
{
  Gce gce{ {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}, 10 };
  ASSERT_EQ(gce.getHistogram().size(), 10);
  ASSERT_EQ(gce.getGridLength(), 11);
  ASSERT_EQ(gce.getGrid().size(), 11);
  ASSERT_EQ(gce.getCenters().size(), 10);
  ASSERT_DOUBLE_EQ(gce.getCenters().at(0), -0.04);
  ASSERT_DOUBLE_EQ(gce.getCenters().at(9), 1.04);
  gce.calculate();
  ASSERT_EQ(gce.getDensityEstimation().size(), 10);
  EXPECT_DOUBLE_EQ(gce.getDensityEstimation()[0] / gce.getDensityEstimation()[4], 0.5);
  EXPECT_DOUBLE_EQ(gce.integrate("Simpson", 0.0, 1.0), 0.74242424242424232);
}
