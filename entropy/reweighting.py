import numpy as np


# Boltzmann constant in kCal/(mol*K):
KB_kcal = 0.001987204118


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
    # going to the nth order includes [0,n]!!
    return np.sum([xs ** i / np.math.factorial(i) for i in np.arange(mac_order+1)], axis=0)


def calculate_amd_weight_maclau(boostEne, mac_order=10, T=300.0, kb=KB_kcal):
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
    return maclaurin_series(scaled_biasE, mac_order=mac_order)


def calculate_amd_weight_exp(boostEne, T=300.0, kb=KB_kcal):
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
    # TODO
    print("WIP")

    # scaled_biasE = np.array(boostEne) / (T * kb)
    #
    # return maclaurin_series(scaled_biasE, mac_order=mac_order)
    pass
