
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


def start_end_from_grid(grid):
    return np.nanmin(grid), np.nanmax(grid)


def process_shapes(data, weights=None):
    if not (weights is None):
        data, weights = np.array(data), np.array(weights)
        if data.shape != weights.shape:
            err_msg = "Shapes of data and weights is inconsistent!\n" \
                      "data: {}, weights: {}".format(data.shape,weights.shape)
            raise ValueError(err_msg)

    if len(data.shape) == 0:
        err_msg = "Shape of data is suspicious\n" \
                  "{}".format(data.shape)
        raise ValueError(err_msg)
    elif len(data.shape) == 1:
        pass  # CONTINUE HERE
    elif len(data.shape) == 2:
        pass  # CONTINUE HERE
    else:
        err_msg = "Shape of data is suspicious\n" \
                  "{}".format(data.shape)
        raise ValueError(err_msg)



class dihedralEntropy(object):

    def __init__(self, data, weights=None, resolution=4096, verbose=False, method="Simpson"):
        # flags and output
        self.__pdf = None
        self.__pdf_x = None
        self.__bandwidth = None
        self.__kde_is_calculated = False
        self.__entropy_is_calculated = False
        self.__entropies = None
        # input
        self.__data = preprocess_dihedral_data(data)
        weights, weight_switch = process_weights_argument(weights, verbose=verbose)
        self.__has_weights = weight_switch
        self.__weights = weights
        self.__method = process_method_argument(method)
        self.__resolution = process_resolution_argument(resolution)
        self.__verbose = verbose

    def calculate(self, resolution=4096, method=None, verbose=None, id_gas=8.314):
        """Calculate the dihedral entropy of a set of dihedrals.
        # TODO Docstring
        The dihedral entropy of a number of different dihedral angles can be calculated using this
        function. The output will be a number for the entropy in each direction. This is calculated
        using a generalized cross entropy method, published by Y.Botev et al.

        >>> calculateEntropy([[1] * 600000, [1] * 600000])
        [0.00627780151892203, 0.00627780151892203]

        Parameters
        ----------
        resolution: int, optional
            The resolution for the estimation of the probability density function.
            This value is applied on the mirrored data. Therefore the resolution in
            real space will be $resolution / 3.
            If no power of two is given, the next higher power of two
            is picked.(default is 4,096)
        method: str, optional
            The method for the numerical integration scheme. Can be one of
            "Riemann" or "Simpson" (default is "Simpson")
        verbose: bool
        id_gas: float
            ideal gas constants in J/(mol*K)
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
        method = method or process_method_argument(method)
        if verbose:
            print("Using the following method for integration: {}.".format(method))

        if verbose:
            print("Initializing C++ kernel for kde...")
        if self.has_weights:
            k = kernel_weights(self.data, self.weights, self.resolution)
            if verbose:
                print("Weights have been given for the calculation of the histograms.")
        else:
            k = kernel(self.data, self.resolution)
        k.calculate()
        self.set_kde_is_calculated(True)
        if verbose:
            print("KDE finished.")
        self.set_pdf_x(k.get_grid())
        self.set_bandwidth(k.get_bandwidth())
        self.set_pdf(k.get_pdf())
        start, end = start_end_from_grid(self.pdf_x)
        integral = k.integrate(start, end, method=method)
        entropies = integral * id_gas
        self.set_entropies(entropies)
        self.set_entropy_is_calculated(True)
        ### old
        # values = []
        #
        # if isinstance(dihedralArr[0], float):
        #     dihedralArr = [dihedralArr]
        #
        # for dihedrals in dihedralArr:
        #     # The mirroring is not necessary, all of this is done via the ent module
        #     entropyCalculator = ent(list(dihedrals), resolution, method)
        #     values.append(entropyCalculator.getResult() * -1)
        #
        # return values

        return entropies

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

    def get_kde_is_calculated(self):
        return self.__kde_is_calculated

    def get_method(self):
        return self.__is_finished

    def get_pdf_x(self):
        if not self.kde_is_calculated:
            self.calculate()
        return self.__pdf_x

    def get_bandwidth(self):
        if not self.kde_is_calculated:
            self.calculate()
        return self.__bandwidth

    def get_pdf(self):
        if not self.kde_is_calculated:
            self.calculate()
        return self.__pdf

    def get_entropies(self):
        if not self.entropy_is_calculated:
            self.calculate()
        return self.__entropies

    def get_entropy_is_calculated(self):
        return self.__entropy_is_calculated

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

    def set_method(self, value):
        self.__method = value

    def set_kde_is_calculated(self, value):
        self.__kde_is_calculated = value

    def set_pdf_x(self, value):
        self.__pdf_x = value

    def set_bandwidth(self, value):
        self.__bandwidth = value

    def set_pdf(self, value):
        self.__pdf = value

    def set_entropies(self, value):
        self.__entropies = value

    def set_entropy_is_calculated(self, value):
        self.__entropy_is_calculated = value

    resolution = property(get_resolution, set_resolution, None, "resolution for kde")
    verbose = property(get_verbose, set_verbose, None, "Extend of print messages")
    data = property(get_data, set_data, None, "Data to do kde on.")
    weights = property(get_weights, set_weights, None, "Weights for the reweighting.")
    has_weights = property(get_has_weights, set_has_weights, None, "Flag, whether weights for the reweighting "
                                                                   "have been given for initialization.")
    kde_is_calculated = property(get_kde_is_calculated, set_kde_is_calculated, None, "Has the kde calculation "
                                                                                     "already been done?")
    entropy_is_calculated = property(get_entropy_is_calculated, set_entropy_is_calculated, None, "Has the entropy "
                                                                                                 "calculation already "
                                                                                                 "been done?")
    bandwidth = property(get_bandwidth, set_bandwidth, None, "bandwidth")
    pdf_x = property(get_pdf_x, set_pdf_x, None, "Grid for probability density function.")
    pdf = property(get_pdf, set_pdf, None, "Probability density function.")
    method = property(get_method, set_method, None, "Method for integrating the probability density function.")
    entropies = property(get_entropies, set_entropies, None, "Calculated entropies for the data sets.")
