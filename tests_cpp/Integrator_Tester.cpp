#include "pch.h"
#include "../src/Integrators.h"


TEST(Integrator_Test, Riemann_Integral)
{
  std::vector<double> f1{ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> f2{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> f3{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> f4{ 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  EXPECT_DOUBLE_EQ(Riemann()(f1, 1.0), 0.1);
  EXPECT_DOUBLE_EQ(Riemann()(f2, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(Riemann()(f3, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(Riemann()(f4, 1.0), 0.2);
}

TEST(Integrator_Test, Simpson_Integral)
{
  std::vector<double> f1{ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> f2{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::vector<double> f3{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  std::vector<double> f4{ 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // First test case is 4.0 / 30.0, this happens due to the term 4 * x_{2j - 1}
  // You can check the calculation, the index is not even, so you get 4 * 1.0
  // This is multiplied by the range (1.0) and divided by three times the number
  // of sub-intervals (i.e., 30)
  EXPECT_DOUBLE_EQ(Simpson()(f1, 1.0), 0.13333333333333333); 
  EXPECT_DOUBLE_EQ(Simpson()(f2, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(Simpson()(f3, 1.0), 1.0);
  EXPECT_DOUBLE_EQ(Simpson()(f4, 1.0), 0.2);
}

TEST(Integrator_Test, Find_Integrator)
{
  EXPECT_THROW(getFunction("Random"), UnknownIntegrator);
  // Weird people might do that!
  auto test = getFunction("SiMpSoN");
  EXPECT_NO_THROW(dynamic_cast<Simpson&>(*test));
  test = getFunction("Riemann");
  EXPECT_NO_THROW(dynamic_cast<Riemann&>(*test));
}