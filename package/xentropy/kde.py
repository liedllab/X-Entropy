"""
part of the entroPy module

@author: paq
"""
import numpy as np
from xentropy.kde_kernel import _kde_kernel
from .internal.resolution import process_resolution_argument
from .internal.pre_post_processing import sanity_check_input_data, process_weights_argument, reshape_arrays_eventually

def check_dims(data):
    data = np.squeeze(data)
    if len(np.shape(data))>1:
        err_msg = "Kde is designed for single data sets only!\n" \
                  "You provided data of the following shape\n" \
                  "data: {}".format(data.shape)
        raise ValueError(err_msg)

class kde(object):

    def __init__(self, data, weights=None, resolution="auto", verbose=False):
        # to be set after calculation
        self.__is_finished = False
        self.__pdf = None
        self.__pdf_x = None
        self.__bandwidth = None
        # input
        # TODO preprocess and sanitycheck data here, too (code in dihedrals.py should be applicable)
        weights, weight_switch = process_weights_argument(weights, verbose=verbose)
        sanity_check_input_data(data, weights, weight_switch)
        self.__has_weights = weight_switch
        check_dims(data)  # kde can only take single data sets currently
        self.__data = data
        self.__weights = weights
        self.__verbose = verbose
        self.__resolution = process_resolution_argument(resolution, self.data, verbose=self.verbose)

    def calculate(self, resolution=None, verbose=None):
        """Perform the kde

        Parameters
        ----------
        verbose: bool
            You may overwrite the kde.verbose property for this function specifically
        resolution: int
            You may specifically reset the resolution for the kde calculation.

        Returns
        -------

        """
        if not (verbose is None):  # it is possible to overwrite kde.verbose for this subroutine
            verbose = verbose
        else:
            verbose = self.verbose
        if not (resolution is None):
            new_res = process_resolution_argument(resolution, self.data, verbose=verbose)
            if verbose:
                print("Using resolution of {}".format(new_res))
            self.__resolution = new_res

        if verbose:
            print("Initializing C++ kernel for kde...")
        if self.has_weights:
            # print(self.weights)  # DEBUG
            kernel = _kde_kernel(self.data, self.resolution, self.weights)
        else:
            kernel = _kde_kernel(self.data, self.resolution)
        kernel.calculate()
        # print(kernel.get_pdf()) # DEBUG
        self.__is_finished = True
        if verbose:
            print("KDE finished.")
        self.__pdf_x = reshape_arrays_eventually(kernel.get_grid())
        self.__bandwidth = kernel.get_bandwidth()
        self.__pdf = reshape_arrays_eventually(kernel.get_pdf())

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
        self.__resolution = process_resolution_argument(value, self.data)

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
