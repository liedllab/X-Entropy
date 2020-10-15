
"""
part of the entroPy module

@author: paq
"""
import numpy as np
from entropy.kde_kernel import _kde_kernel
from .resolution import process_resolution_argument
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
    # TODO mirroring and so on...
    return data


def start_end_from_grid(grid):
    return np.nanmin(grid), np.nanmax(grid)


def process_data_shapes(data, weights=None, weight_switch=True):
    data = np.array(data)
    if False in np.isfinite(data):
        err_msg = "Non-finite values in data!"
        raise ValueError(err_msg)

    if weight_switch:
        weights = np.array(weights)
        if False in np.isfinite(weights):
            err_msg = "Non-finite values in weights!"
            raise ValueError(err_msg)
        if data.shape != weights.shape:
            err_msg = "Shapes of data and weights is inconsistent!\n" \
                      "data: {}, weights: {}".format(data.shape,weights.shape)
            raise ValueError(err_msg)

    if len(data.shape) == 0:
        err_msg = "Shape of data is suspicious\n" \
                  "{}".format(data.shape)
        raise ValueError(err_msg)
    elif len(data.shape) == 1:
        data = np.array([data])
        if weight_switch:
            weights = np.array([weights])
    elif len(data.shape) == 2:
        pass  # This ia all good
    else:
        err_msg = "Shape of data is suspicious\n" \
                  "{}".format(data.shape)
        raise ValueError(err_msg)
    return data, weights


class dihedralEntropy(object):

    def __init__(self, data, weights=None, resolution=4096, verbose=False, method="Simpson"):
        # flags and output
        self.__pdf = None
        self.__pdf_x = None
        self.__bandwidth = None
        self.__kde_is_calculated = False
        self.__entropy_is_calculated = False
        self.__entropies = None
        # input data
        weights, weight_switch = process_weights_argument(weights, verbose=verbose)
        self.__has_weights = weight_switch
        data, weights = process_data_shapes(data, weights, weight_switch)
        self.__data = preprocess_dihedral_data(data)
        self.__weights = weights
        # other input
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
            if verbose:
                print("Weights have been given for the calculation of the histograms.")
            iterable = zip(self.data, self.weights)
        else:
            iterable = zip(self.data, np.full(None,len(self.data)))

        bws, pdfxs, pdfs, ents, ints = [], [], [], [], []
        for dat, ws in iterable:
            # depending on whether weights are None or not, you will get a weighted pdf or a simple pdf
            kernel = _kde_kernel(data, weights, resolution)
            kernel.calculate()

            bws.append(kernel.get_bandwidth())
            pdfxs.append(kernel.get_grid())
            pdfs.append(kernel.get_pdf())

            start, end = start_end_from_grid(self.pdf_x)
            integral = kernel.integrate(start, end, method=method)
            ints.append(integral)
            entropy = integral * id_gas
            ents.append(entropy)
        if verbose:
            print("KDE finished.")

        self.__pdf_x = np.array(pdfxs)
        self.__bandwidth = np.array(bws)
        self.__pdf = np.array(pdfs)
        self.set_entropies(np.array(ent))

        self.set_kde_is_calculated(True)
        self.set_entropy_is_calculated(True)
        return np.array(ents)

    @property
    def resolution(self):
        """Resolution for kde. Needs to be a power of 2. If no power of two is given,
        the next higher power of two is set automatically."""
        return self.__resolution

    @resolution.setter
    def resolution(self, value):
        if self.is_finished:
            print("After changing the resolution you should use .calculate() again, before accessing any results...\n"
                  "It is probably better to explicitly call calculate with a specific resolution.")
        self.__resolution = process_resolution_argument(value)

    @property
    def verbose(self):
        """Extend of print messages. """
        return self.__verbose

    @verbose.setter
    def verbose(self, value):
        print("Verbosity cannot be changed after initialization...")
        pass

    @property
    def is_finished(self):
        """Is the kde finished? """
        return self.__is_finished

    @is_finished.setter
    def is_finished(self, value):
        print("You really shouldn't change this flag yourself. Use .calculate()")
        pass

    @property
    def has_weights(self):
        """Have weights been given?"""
        return self.__has_weights

    @has_weights.setter
    def has_weights(self, value):
        print("This flag cannot be changed after initialization.")
        pass

    @property
    def bandwidth(self):
        """Bandwidth from the kde. """
        if not self.is_finished:
            self.calculate()
        return self.__bandwidth

    @bandwidth.setter
    def bandwidth(self, value):
        print("You really shouldn't change this property yourself. Use .calculate()")
        pass

    @property
    def pdf(self):
        """Probability density function. """
        if not self.is_finished:
            self.calculate()
        return self.__pdf

    @pdf.setter
    def pdf(self, value):
        print("You really shouldn't change .pdf yourself. Use .calculate()")
        pass

    @property
    def data(self):
        """"Data to do kde on. """
        return self.__data

    @data.setter
    def data(self, value):
        print("Data cannot be changed after initialization...")
        pass

    @property
    def weights(self):
        """Weights for the pdf calculation. """
        return self.__weights

    @weights.setter
    def weights(self, value):
        print("Weights cannot be changed after initialization...")
        pass

    @property
    def pdf_x(self):
        """Grid for probability density function. """
        if not self.is_finished:
            self.calculate()
        return self.__pdf_x

    @pdf_x.setter
    def pdf_x(self, value):
        print("You really shouldn't change .pdf_x yourself. Use .calculate()")
        pass

    # Getter #
    def get_kde_is_calculated(self):
        return self.__kde_is_calculated

    def get_method(self):
        return self.__is_finished

    def get_entropies(self):
        if not self.entropy_is_calculated:
            self.calculate()
        return self.__entropies

    def get_entropy_is_calculated(self):
        return self.__entropy_is_calculated

    # setter #

    def set_method(self, value):
        self.__method = value

    def set_kde_is_calculated(self, value):
        self.__kde_is_calculated = value

    def set_entropies(self, value):
        self.__entropies = value

    def set_entropy_is_calculated(self, value):
        self.__entropy_is_calculated = value

    kde_is_calculated = property(get_kde_is_calculated, set_kde_is_calculated, None, "Has the kde calculation "
                                                                                     "already been done?")
    entropy_is_calculated = property(get_entropy_is_calculated, set_entropy_is_calculated, None, "Has the entropy "
                                                                                                 "calculation already "
                                                                                                 "been done?")
    method = property(get_method, set_method, None, "Method for integrating the probability density function.")
    entropies = property(get_entropies, set_entropies, None, "Calculated entropies for the data sets.")
