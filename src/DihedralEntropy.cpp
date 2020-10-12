#include "DihedralEntropy.h"

/********************************************************************************
 * 										                                                          *
 * Python ports.								                                                *
 *										                                                          *
 ********************************************************************************/

namespace py = boost::python;

DihedralEntropy::DihedralEntropy(py::list &l, int n) : m_entropy(0), m_res(n) {
  for (int i = 0; i < static_cast<int>(py::len(l)); ++i) {
    m_angles.push_back(py::extract<double>(l[i]));
  }
  m_length = static_cast<int>( m_angles.size() );
  integrate();
}

DihedralEntropy::DihedralEntropy(py::list &l, int n, py::str &numericalIntegral)
: m_entropy{ 0.0 }, m_res{ n } {
  for (int i = 0; i < static_cast<int>( py::len(l) ); ++i) {
    m_angles.push_back(py::extract<double>(l[i]));
  }
  m_length = m_angles.size();
  std::string numInt{ py::extract<char *>(numericalIntegral) };
  integrate(getFunction(numInt));
}

DihedralEntropy::DihedralEntropy(py::list &l, py::str &numericalIntegral) : 
  DihedralEntropy{ l, (2 << 12), numericalIntegral } {}

/****
 * Used to calculate the integral of an dihedral angle distribution. If another
 * Entropy is desired, the GCE (or kde in python) has to be used directly and
 * everything encoded here has to be done by python. Or you are also welcome to
 * just change the Code here and send me the copy ;)
 */
void DihedralEntropy::integrate(const std::unique_ptr<IIntegration>& inte) {
  std::vector<double> mirrored(m_angles.size() * 3.0);
  double dx{ 0.0 };
  double norm{ 0.0 };
  // Mirrors the dihedrals to get rid of boundary problems
  #pragma omp parallel for
  for (int i = 0; i < (int) m_angles.size(); ++i) {
    mirrored.at(i) = m_angles.at(i) - 360.0;
    mirrored.at(i + m_angles.size()) = m_angles.at(i);
    mirrored.at(i + 2 * m_angles.size()) = m_angles.at(i) + 360.0;
  }
  // Uses the GCE to calculate the Density.
  Gce kde{ mirrored, m_res };
  kde.calculate();
  dx = kde.getGrid().at(1) - kde.getGrid().at(0);
  // This is defined here, because the handling of those is easier.
  auto density = kde.getDensityEstimation();
  auto xgrid{ kde.getGrid() };
  std::vector<double> finalDens;
  // Integrate only over the angles between -180 and 180 degrees.
  for (int i = 0; i < (int) density.size(); ++i) {
    if (xgrid.at(i) >= -180.0 && xgrid.at(i) <= 180.0) {
      finalDens.push_back(density.at(i));
    }
    if (xgrid.at(i) > 180.0) {
      break;
    }
  }
  // Integration part.
  norm = (*inte)(finalDens, 360);
  //for (int i = 0; i < (int) finalDens.size(); ++i) {
  //  norm += finalDens.at(i);
  //}

  OMPExceptionHandler except;

  #pragma omp parallel for
  for (int i = 0; i < (int) finalDens.size(); ++i) {

    except.Run([&] {
      finalDens.at(i) /= norm;
      if (finalDens.at(i) > 0) {
        finalDens.at(i) = finalDens.at(i) * log(finalDens.at(i));
      } else if (finalDens.at(i) != finalDens.at(i) ) {
        throw IntegrationError("Created a NAN during integration!");
      }
    });
  }
  except.Rethrow();
  m_entropy = (*inte)(finalDens, 360);
  // Finally done, just multiply with the universal gasconstant 
  // (which you will find to be defined in the header)
  m_entropy *= GASCONSTANT;

}

void DihedralEntropy::integrate() {
  integrate(std::make_unique<Riemann>());
}

/****
 * Get the entropy.
 * @return: The value for the entropy.
 */
double DihedralEntropy::getEntropy() {
  return m_entropy;
}

