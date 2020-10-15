
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


def preprocess_dihedral_data(diheds, weights, weight_switch):
    # mirror data
    diheds_out = []
    for dihed in diheds:
        diheds_out.append(np.concatenate([dihed-360, dihed, dihed+360]))
    diheds = np.array(diheds_out)
    if weight_switch:  #
        weights_out = []
        for weight in weights:
            weights_out.append(np.conatenate([weight, weight, weight]))
        weights = np.array(weights_out)

    return diheds, weights


def postprocess_dihedral_pdf(pdf, pdf_x):
    # mirror data
    lower_idx = np.argmin(np.abs(pdf_x + 180))
    if pdf_x[lower_idx] < -180:
        lower_idx = lower_idx + 1
    upper_idx = np.argmin(np.abs(pdf_x - 180))
    if pdf_x[upper_idx] > 180:
        upper_idx = upper_idx - 1

    pdf_out = pdf[lower_idx:upper_idx]
    pdf_out *= 3  # we normalized in the "mirrored data", which is 3 times the actual data
    pdf_x_out = pdf_x[lower_idx:upper_idx]
    assert len(pdf_out)==len(pdf_x_out)
    return pdf_out, pdf_x_out


def start_end_from_grid(grid):
    # You will need this for non dihedrals...
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
        self.__pdfs = None
        self.__pdf_xs = None
        self.__bandwidth = None
        self.__is_finished = False
        self.__entropies = None
        # input data
        weights, weight_switch = process_weights_argument(weights, verbose=verbose)
        self.__has_weights = weight_switch
        data, weights = process_data_shapes(data, weights, weight_switch)
        data, weights = preprocess_dihedral_data(data, weights, weight_switch)
        self.__data = data
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
        if not (verbose is None):  # it is possible to overwrite kde.verbose for this subroutine
            verbose = verbose
        else:
            verbose = self.verbose

        if not (resolution is None):  # reset resolution eventually
            new_res = process_resolution_argument(resolution)
            self.__resolution = new_res
            if verbose:
                print("Using resolution of {}".format(self.resolution))

        if not (method is None):    # reset method eventually
            new_method = process_method_argument(method)
            self.__method = new_method
            if verbose:
                print("Using the following method for integration: {}.".format(self.method))

        if verbose:
            print("Initializing C++ kernel for kde...")
        # if weights are given, zip with weights, if not, give None
        if self.has_weights:
            if verbose:
                print("Weights have been given for the calculation of the histograms.")
            iterable = zip(self.data, self.weights)
        else:
            iterable = zip(self.data, np.full(len(self.data), None))

        bws, pdf_xs, pdfs, ents = [], [], [], []
        for dat, ws in iterable:
            # depending on whether weights are None or not, you will get a weighted pdf or a simple pdf
            kernel = _kde_kernel(dat, self.resolution, ws)
            kernel.calculate()

            bws.append(kernel.get_bandwidth())

            pdf_temp = kernel.get_pdf()
            pdf_x_temp = kernel.get_grid()

            pdf_temp, pdf_x_temp = postprocess_dihedral_pdf(pdf_temp, pdf_x_temp)
            # start, end = start_end_from_grid(pdf_x_temp)

            pdf_xs.append(pdf_x_temp)
            pdfs.append(pdf_temp)

            # integral = kernel.integrate(start, end, method=self.method)
            p_logp = kernel.calculate_entropy(-180, 180, method=self.method)

            entropy = - p_logp * id_gas
            ents.append(entropy)
        if verbose:
            print("KDE finished.")

        self.__pdf_xs = np.array(pdf_xs)
        self.__bandwidth = np.array(bws)
        self.__pdfs = np.array(pdfs)
        self.__entropies = np.array(ents)

        self.__is_finished = True
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
    def pdfs(self):
        """Probability density function. """
        if not self.is_finished:
            self.calculate()
        return self.__pdfs

    @pdfs.setter
    def pdfs(self, value):
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
    def pdf_xs(self):
        """Grid for probability density function. """
        if not self.is_finished:
            self.calculate()
        return self.__pdf_xs

    @pdf_xs.setter
    def pdf_xs(self, value):
        print("You really shouldn't change .pdf_x yourself. Use .calculate()")
        pass

    #getter

    def get_method(self):
        return self.__method

    def get_entropies(self):
        if not self.is_finished:
            self.calculate()
        return self.__entropies

    # setter #

    def set_method(self, value):
        self.__method = value

    def set_entropies(self, value):
        self.__entropies = value


    method = property(get_method, set_method, None, "Method for integrating the probability density function.")
    entropies = property(get_entropies, set_entropies, None, "Calculated entropies for the data sets.")
