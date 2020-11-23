#include "Integrators.h"

std::unique_ptr<IIntegration> getFunction(const std::string &type)
{
    std::string newstr{type};
    transform(newstr.begin(), newstr.end(), newstr.begin(), ::tolower);
    if (newstr == "simpson" || newstr == "kepler")
    {
        return std::make_unique<Simpson>();
    }
    else if (newstr == "riemann")
    {
        return std::make_unique<Riemann>();
    }

    throw UnknownIntegrator("The function you searched for is not yet implemented.");
}

double Simpson::operator()(const std::vector<double> &function, double range) const noexcept
{
    // TODO: Should actually throw an error
    int steps{static_cast<int>(function.size())};
    if (steps % 2)
    {
        steps -= 1;
    }
    double ret = function.at(0);
    double h = range / (3 * steps);
    ret += function.at(function.size() - 1);
    double calc = 0.0;

// C
#pragma omp parallel for reduction(+ \
                                   : calc)
    for (int j = 1; j <= (steps / 2); ++j)
    {
        calc += function.at(2 * j - 1);
    }

    ret += 4.0 * calc;
    calc = 0.0;

// Check this!
#pragma omp parallel for reduction(+ \
                                   : calc)
    for (int j = 1; j <= ((steps / 2) - 1); ++j)
    {
        calc += function.at(2 * j);
    }

    ret += 2.0 * calc;
    return ret * h;
}

double Riemann::operator()(const std::vector<double> &function, double range) const noexcept
{
    double dx{range / static_cast<double>(function.size())};
    double ret{0};
    for (auto val : function)
    {
        ret += val * dx;
    }
    return ret;
}