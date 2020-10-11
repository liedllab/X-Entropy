from entropy.kde import DihedralEntropy as ent
from entropy.kde import GCE
import numpy as np
import warnings

# We want to to change that default, since ignoring warnings is ultimately the users decision:
warnings.simplefilter("always")


def calculateEntropy(dihedralArr, resolution=2160, method="Simpson"):
    """Calculate the dihedral entropy of a trajectory.

    The dihedral entropy of a number of different dihedral angles can be calculated using this
    function. The output will be a number for the entropy in each direction. This is calcualted
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

# TODO temperature is of no use here
def calculateReweightedEntropy(dihedralArr, weightArr, resolution = 16000, method = "Simpson", mc_order=10, temp=300):
    """Calculate the dihedral entropy of an accelerated MD (aMD) trajectory using Maclaurin series
       expansion for reweighting (see https://doi.org/10.1021/ct500090q).

    Parameters
    ----------
    dihedralArr: list(list(float)) or list(float)
        A 2D-array, holding all the dihedrals of the trajectory.
    weightArr: list(list(double))
        The weights of the aMD trajctory.
    resolution: int, optional
        The resolution of the initial Histogram approximation (default is 16000)
    method: str, optional
        The method for the numerical integral, at the moment only
        Riemann and Simpson is available (default is Simpson)
    mc_order: int, optional
        The order of the MCLaurin series expansion (default is 10)
    temp: float, optional
        The temperature for the reweighting (default is 300)

    Returns
    -------
    values:
    A list holding the values for each angle.
    """

    values = []

    if isinstance(dihedralArr[0], (float)):
        dihedralArr = [dihedralArr]

    for dihedral in dihedralArr:
        inHist = reweighting(dihedral, weightArr, mc_order=mc_order, temp=temp, resolution=resolution)
        kernel = GCE(list(inHist), list(np.linspace(-180 - 360, 180 + 360, num = resolution)), len(dihedral))
        kernel.calculate()
        values.append(kernel.integrate(method, -180, 180) * -8.3145)

    return values


# TODO temp is not used in the function below.
def reweighting (diheds, weights, mc_order = 10, temp = 300, binsX = None, resolution=(2 << 12)):
    """Reweight the histogram of a dihedral with a Maclaurin series expansion (https://doi.org/10.1021/ct500090q). 

    Parameters
    ----------
    diheds: list(list(float))
        A 2D-array, holding a dihedral
    weights: list(list(float))
        The weights of the aMD trajectory
    mc_order: int, optional
        The order of the MCLaurin series (default is 10)
    temp: float, optional
        The simulation temperature (default is 300)
    binsX: list or None, optional
        The bins used for histogramming (default is None)

    Returns
    -------
    The reweighted histogram of the dihedral.
    """

    def mirror (arr):
      return list(arr) + list(arr) + list(arr)
    def fact(x):
      if x == 0:
        return 1
      return fact(x - 1) * x

    diheds = np.array(mirror(diheds))
    weights = np.array(mirror(weights))
    if binsX is None:
        binsX = np.linspace(-180 - 360, 180 + 360, num = resolution)
    MCweight = np.ones(len(weights))

    for x in range(mc_order):
        MCweight = np.add(MCweight, np.divide(np.power(weights, x + 1), fact(x + 1)))
    hist = np.histogram(diheds, bins = binsX, weights = MCweight)[0]
    norm = np.linalg.norm(hist)
    hist = np.divide(hist, norm)
    return hist
