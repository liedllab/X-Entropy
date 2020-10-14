
"""
part of the entroPy module

@author: paq
"""
import numpy as np
from entropy.kde_kernel import __kde_kernel as kernel
from .resolution import process_resolution_argument
from .reweighting import reweighting, calculate_amd_weight
import warnings

# We want to to change that default, since ignoring warnings is ultimately the users decision:
warnings.simplefilter("always")


def process_method_argument(method, available_methods=("simpson", "riemann")):
    """Catches unknown integration methods on a python
    level. Deals with upper and lower case writing.

    """
    if not (method.lower() in available_methods):
        err_msg = "{} is not a valid integration method\nChoose from {}.".format(method, available_methods)
        raise ValueError(err_msg)
    return method.lower().capitalize()  # first letter caps, rest lower case


def process_weights_argument(weights, verbose=False):
    """This function will mainly return a switch, whether
    or whether not weights will be passed on to the kde"""
    weight_switch = True
    if weights is None:
        if verbose:
            print("No weights been given.")
        weight_switch = False
    else:
        if verbose:
            print("Weights have been given.")
        weights = np.array(weights)

    return weights, weight_switch


def preprocess_dihedral_data(data):
    # TODO
    return data


class dihedralEntropy(object):
    __bins = None
    __data = None
    __weights= None
    __has_weights = None
    __resolution = None
    __kde_done = False
    __successful = None
    __verbose = None
    __pdf = None
    __pdf_x = None
    __bandwidth = None
    __is_finished = None

    def __init__(self, data, weights=None, resolution=4096, verbose=False):
        self.__data = preprocess_dihedral_data(data)
        weights, weight_switch = process_weights_argument(weights, verbose=verbose)
        self.__has_weights = weight_switch
        self.__weights = weights
        self.__resolution = process_resolution_argument(resolution)
        self.__verbose = verbose
        self.__is_finished = False

    def calculate_legacy(self, verbose=None, ):
        """
        copy from kde module, to sneak a peak for me

        """

        if verbose:
            print("Initializing C++ kernel for kde...")
        k = kernel(self.data, self.resolution)
        k.calculate()
        self.set_is_finished(True)
        if verbose:
            print("KDE finished.")
        self.set_pdf_x(center_grid(k.get_grid()))
        self.set_bandwidth(k.get_bandwidth())
        self.set_pdf(k.get_pdf())

    def calculate(self, dihedralArr, resolution=4096, method="Simpson", verbose=False):
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
            The resolution for the estimation of the probability density function.
            This value is applied on the mirrored data. Therefore the resolution in
            real space will be $resolution / 3.
            If no power of two is given, the next higher power of two
            is picked.(default is 4,096)
        method: str, optional
            The method for the numerical integration scheme. Can be one of
            "Riemann" or "Simpson" (default is "Simpson")

        Returns
        -------
        list:
        A list of floats that are the entropies for the different dihedrals
        """
        verbose = verbose or self.verbose
        if not (resolution is None):
            new_res = process_resolution_argument(resolution)
            if verbose:
                print("Using resolution of {}".format(new_res))
            self.set_resolution(new_res)
        method = process_method_argument(method)
        if verbose:
            print("Using the following method for integration: {}.".format(method))

        if self.has_weights:
            k = kernel_weights(self.data, self.weights, self.resolution)
        else:
            k = kernel(self.data, self.resolution)
        k.calculate()
        # integrate and calculat entropy from that...

        # old
        values = []

        if isinstance(dihedralArr[0], float):
            dihedralArr = [dihedralArr]

        for dihedrals in dihedralArr:
            # The mirroring is not necessary, all of this is done via the ent module
            entropyCalculator = ent(list(dihedrals), resolution, method)
            values.append(entropyCalculator.getResult() * -1)

        return values

    # Getter #
    def get_resolution(self):
        return self.__resolution

    def get_weights(self):
        return self.__weights

    def get_has_weights(self):
        return self.__has_weights

    def get_verbose(self):
        return self.__verbose

    def get_data(self):
        return self.__data

    def get_is_finished(self):
        return self.__is_finished

    def get_pdf_x(self):
        if not self.is_finished:
            self.calculate()
        return self.__pdf_x

    def get_bandwidth(self):
        if not self.is_finished:
            self.calculate()
        return self.__bandwidth

    def get_pdf(self):
        if not self.is_finished:
            self.calculate()
        return self.__pdf

    # setter #
    def set_resolution(self, value):
        self.__resolution = value

    def set_weights(self, value):
        print("Weights cannot be changed after initialization...")
        pass

    def set_has_weights(self, value):
        print("This flag cannot be changed after initialization...")
        pass

    def set_verbose(self, value):
        print("Verbosity cannot be changed after initialization...")
        pass

    def set_data(self, value):
        print("Data cannot be changed after initialization...")
        pass

    def set_is_finished(self, value):
        self.__is_finished = value

    def set_pdf_x(self, value):
        self.__pdf_x = value

    def set_bandwidth(self, value):
        self.__bandwidth = value

    def set_pdf(self, value):
        self.__pdf = value

    resolution = property(get_resolution, set_resolution, None, "resolution for kde")
    verbose = property(get_verbose, set_verbose, None, "Extend of print messages")
    data = property(get_data, set_data, None, "Data to do kde on.")
    weights = property(get_weights, set_weights, None, "Weights for the reweighting.")
    has_weights = property(get_has_weights, set_has_weights, None, "Flag, whether weights for the reweighting "
                                                                   "have been given for initialization.")
    is_finished = property(get_is_finished, set_is_finished, None, "Data to do kde on.")
    bandwidth = property(get_bandwidth, set_bandwidth, None, "bandwidth")
    pdf_x = property(get_pdf_x, set_pdf_x, None, "Pprobability densite sxfvsj_x ")
    pdf = property(get_pdf, set_pdf, None, "Pprobability densite sxfvsj ")

