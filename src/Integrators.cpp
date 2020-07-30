#include "Integrators.h"

std::unique_ptr<IIntegration> getFunction(const std::string& type) {
  if (type == "Simpson") {
    return std::make_unique<Simpson>();
  } else if (type ==  "Riemann") {
      return std::make_unique<Riemann>();
  }

  throw UnknownIntegrator( "The function you searched for is not yet implemented." );
}


double Simpson::operator() (const std::vector<double> &function, int steps, double range) const noexcept{
  if (steps % 2) {
    steps -= 1;
  }
  double ret = function.at(0);
  double h = range / (3 * steps);
  ret += function.at(function.size() - 1);
  double pcalc[omp_get_max_threads()] = { 0 };
  double calc = 0;


  // Check this!!
  #pragma omp parallel for
  for (int j = 1; j < (steps / 2); ++j) {
    pcalc[omp_get_thread_num()] += function.at(2 * j - 1);
  }


  for (int i = 0; i < omp_get_max_threads(); ++i) {
    calc += pcalc[i];
    pcalc[i] = 0;
  }


  ret += 4 * calc;
  calc = 0;

  // Check this!
  #pragma omp parallel for
  for (int j = 1; j < ((steps / 2) - 1); ++j) {
    pcalc[omp_get_thread_num()] += function.at(2 * j);
  }


  for (int i = 0; i < omp_get_max_threads(); ++i) {
    calc += pcalc[i];
  }
  ret += 2 * calc;
  return ret * h;
}

double Riemann::operator() (const std::vector<double> &function, int steps, double range) const noexcept{
		double dx{ range / static_cast<double>( function.size() ) };
		double ret{ 0 };
		for (auto val : function) {
			ret += val * dx;
		}
		return ret;
	}