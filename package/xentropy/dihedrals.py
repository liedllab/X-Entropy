
"""
part of the X-EntroPy module

@author: paq
"""
import numpy as np
from xentropy.kde_kernel import _kde_kernel

from xentropy.internal.resolution import minim_or_sqrt
from .internal.resolution import process_resolution_argument, rules_of_thumb_dihedral

from .internal.pre_post_processing import preprocess_dihedral_data, process_data_shapes, \
    process_weights_argument, process_method_argument, reshape_arrays_eventually, postprocess_dihedral_pdf
from .internal.constants import id_gas_SI, PI
import warnings

# We want to to change that default, since ignoring warnings is ultimately the users decision:
warnings.simplefilter("always")


class dihedralEntropy(object):

    def __init__(self, data, weights=None, resolution="auto", verbose=False, method="Simpson"):
        # flags and output
        self.__verbose = verbose
        self.__pdf = None
        self.__pdf_x = None
        self.__bandwidth = None
        self.__is_finished = False
        self.__entropy = None
        # input data
        # determine resolution from non mirrored data!!
        self.__resolution = process_resolution_argument(resolution, data, verbose=self.verbose,
                                                        rules_of_thumb=rules_of_thumb_dihedral())
        # process input arrays
        weights, weight_switch = process_weights_argument(weights, verbose=verbose)
        self.__has_weights = weight_switch
        data, weights = process_data_shapes(data, weights, weight_switch)
        self.__data = data
        self.__weights = weights
        # other input
        self.__method = process_method_argument(method)

    def calculate(self, resolution=None, method=None, verbose=None, id_gas=id_gas_SI):
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
            new_res = process_resolution_argument(resolution, self.data, verbose=verbose,
                                                  rules_of_thumb=rules_of_thumb_dihedral())
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

        bws, pdf_xs, pdfs, ents = [], [], [], []
        for dat in self.data:
            # apply periodic copies here
            dat, ws = preprocess_dihedral_data(dat, self.weights, self.has_weights)
            # depending on whether weights are None or not, you will get a weighted pdf or a simple pdf
            kernel = _kde_kernel(dat, self.resolution, ws)
            kernel.calculate()

            bws.append(kernel.get_bandwidth())

            pdf_temp = kernel.get_pdf()
            pdf_x_temp = kernel.get_grid()

            norm_for_mirrored_data = kernel.integrate(-PI, PI, method=self.method)
            pdf_temp, pdf_x_temp = postprocess_dihedral_pdf(pdf_temp, pdf_x_temp, norm_for_mirrored_data)
            # start, end = start_end_from_grid(pdf_x_temp)

            pdf_xs.append(pdf_x_temp)
            pdfs.append(pdf_temp)

            # integral = kernel.integrate(start, end, method=self.method)
            p_logp = kernel.calculate_entropy(-PI, PI, method=self.method)

            entropy = - p_logp * id_gas
            ents.append(entropy)
        if verbose:
            print("KDE finished.")

        ents = reshape_arrays_eventually(np.array(ents))
        self.__pdf_x = reshape_arrays_eventually(np.array(pdf_xs))
        self.__bandwidth = np.array(bws)
        self.__pdf = reshape_arrays_eventually(np.array(pdfs))
        self.__entropy = ents

        self.__is_finished = True
        return ents

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
        self.__resolution = process_resolution_argument(value, self.data, rules_of_thumb=rules_of_thumb_dihedral())

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

    @property
    def method(self):
        return self.__method

    @method.setter
    def method(self, value):
        """Method for integrating the probability density function."""
        if self.is_finished:
            print("After changing the method you should use .calculate() again, before accessing any results...\n"
                  "It is probably better to explicitly call calculate with a specific method.")
        self.__method = process_method_argument(value)

    @property
    def entropy(self):
        """Calculated entropies for the data sets."""
        if not self.is_finished:
            self.calculate()
        return self.__entropy

    @entropy.setter
    def entropy(self, value):
        print("You really shouldn't change .entropies yourself. Use .calculate()")
        pass
