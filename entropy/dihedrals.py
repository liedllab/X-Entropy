from entropy.kde import DihedralEntropy as ent
from entropy.kde import GCE
import numpy as np
import warnings

# We want to to change that default, since ignoring warnings is ultimately the users decision:
warnings.simplefilter("always")

# Boltzmann constant in kCal/(mol*K):
KB_kcal = 0.001987204118


def calculateEntropy(dihedralArr, resolution=4096, method="Simpson"):
    """Calculate the dihedral entropy of a trajectory.

    The dihedral entropy of a number of different dihedral angles can be calculated using this
    function. The output will be a number for the entropy in each direction. This is calculated
    using a generalized cross entropy method, published by Y.Botev et al.

    >>> calculateEntropy([[1] * 600000, [1] * 600000])
    [0.00627780151892203, 0.00627780151892203]

    Parameters
    ----------
    dihedralArr: list(list(float)) or list(float)
        A 2D array holding all the dihedrals of a simulation (number Atoms,
        number Dihedrals)
    resolution: int, optional
        The resolution used for the initial construction of the histogram
        bin size (default is 16,000)
    method: str, optional
        The method for the numerical integration scheme. Can be one of
        "Riemann" or "Simpson" (default is "Simpson")

    Returns
    -------
    list:
    A list of floats that are the entropies for the different dihedrals
    """
    if not isinstance(resolution, int):
        print("Resolution is not of type int. Trying to cast it to int...")
        resolution = int(resolution)
    if resolution < 100:
        warn_msg = "You are using a rather small resolution. " \
                   "This may potentially lead to inaccurate results..."
        warnings.warn(warn_msg, RuntimeWarning)
    elif resolution > 10000:
        warn_msg = "You are using a rather large resolution. " \
                   "Amongst other things, this may potentially lead to very long runtimes " \
                   "without necessarily improving the accuracy of the result..."
        warnings.warn(warn_msg, RuntimeWarning)

    values = []

    if isinstance(dihedralArr[0], (float)):
        dihedralArr = [dihedralArr]

    for dihedrals in dihedralArr:
        # The mirroring is not necessary, all of this is done via the ent module
        entropyCalculator = ent(list(dihedrals), resolution, method)
        values.append(entropyCalculator.getResult() * -1)

    return values


def calculateReweightedEntropy(dihedralArr, weightArr, resolution=16000, method="Simpson", id_gas=8.3145):
    """Calculate the dihedral entropy of an accelerated MD (aMD) trajectory using Maclaurin series
       expansion for reweighing (see https://doi.org/10.1021/ct500090q).

    Parameters
    ----------
    dihedralArr: list(list(float)) or list(float)
        A 2D-array, holding all the dihedrals of the trajectory.
    weightArr: list(list(double))
        The weights of the aMD trajectory.
    resolution: int, optional
        The resolution of the initial Histogram approximation (default is 16000)
    method: str, optional
        The method for the numerical integral, at the moment only
        Riemann and Simpson is available (default is Simpson)
    id_gas: float, optional
        Default is 8.314 J/mol/K

    Returns
    -------
    values:
    A list of floats that are the entropies for the different dihedrals in J/mol/K
    """

    values = []

    if isinstance(dihedralArr[0], float):
        dihedralArr = [dihedralArr]

    for dihedral in dihedralArr:
        inHist = reweighting(dihedral, weightArr, resolution=resolution)
        kernel = GCE(list(inHist), list(np.linspace(-180 - 360, 180 + 360, num=resolution)), len(dihedral))
        kernel.calculate()
        values.append(kernel.integrate(method, -180, 180) * - id_gas)

    return values


def reweighting(diheds, weights, bins=None, resolution=(2 << 12)):
    """Reweight the histogram using the given weights.

    Parameters
    ----------
    diheds: list(list(float))
        Array, which holds the data values to be histogrammed
    weights: list(list(float))
        The weights used for the histogram
    bins: list or None, optional
        The bins used for histogramming (default is None)
    resolution: resolution used to construct the bins if bins is None (default 2<<12)

    Returns
    -------
    The reweighed histogram of the input data.
    """

    def mirror(arr):
        return list(arr) + list(arr) + list(arr)

    diheds = np.array(mirror(diheds))
    weights = np.array(mirror(weights))

    if bins is None:
        bins = np.linspace(-180 - 360, 180 + 360, num=resolution)

    hist = np.histogram(diheds, bins=bins, weights=weights)[0]
    norm = np.linalg.norm(hist)
    hist = np.divide(hist, norm)
    return hist


def maclaurin_series(xs, mac_order=10):
    """Maclaurin series to estimate e^{x}.

    Parameters
    ----------
    xs: array of floats
    mc_order: int
        max order for the Maclaurin series

    Returns
    -------

    """
    return np.sum([xs ** i / np.math.factorial(i) for i in np.arange(mac_order)], axis=0)


def calculate_amd_weight(boostEne, mac_order=10, T=300.0, kb=KB_kcal):
    """Calculate weights from aMD simulation boost energies for histogram
    reweighing using Maclaurin series expansion (see https://doi.org/10.1021/ct500090q).

    Paramteres:
    -----------------------------
    boostEne: list(list(float))
        Boost energies from the aMD simulation
    mc_order: int, optional
        The order of the Maclaurin series (default is 10)
    T: float, optional
        The simulation temperature, used to calculate the Boltzmann factor from the energies (default is 300.0)
    kb: float, optional
        Boltzmann constant. Default: 0.001987204118 kCal/(mol*K)

    Returns
    -----------------------------
    Weights for the histogram calculation
    """

    scaled_biasE = np.array(boostEne) / (T * kb)
    # MCweight = np.zeroslike(boostEne)
    # for x in range(mc_order+1):
    #     MCweight = np.add(MCweight, (np.divide(np.power(boostEne, x), math.factorial(x))))
    # return MCweight
    return maclaurin_series(scaled_biasE, mc_order=mac_order)
