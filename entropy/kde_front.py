"""
part of the entroPy module

@author: paq
"""
import numpy as np
from entropy import kde as kde_kernel
from .resolution import process_resolution_argument



class kde(object):
    __bins = None
    __data = None
    __resolution = None
    __kde_done = False
    __successful = None
    __verbose = None
    __pdf = None
    __pdf_x = None
    __is_finished = None

    def __init__(self, data, resolution=4096, verbose=False):
        self.__data = data
        self.__resolution = process_resolution_argument(resolution)
        self.__verbose = verbose
        self.__is_finished = False


    @classmethod
    def generate_bins(self, data):
        pass
        return bins

    def calculate(self, verbose=None):
        """

        Parameters
        ----------
        data
        nbins
        verbose

        Returns
        -------

        """
        verbose = verbose or self.verbose
        if verbose:
            print("Initializing C++ kernel for kde...")
        kernel = kde_kernel.GCE(list(self.data), self.resolution)
        kernel.calculate()
        self.set_is_finished(True)

    def pdf(self, xdata):
        pass

    # Getter #
    def get_resolution(self):
        return self.__resolution

    def get_verbose(self):
        return self.__verbose

    def get_data(self):
        return self.__data

    def get_is_finished(self):
        return self.__is_finished

    # setter #
    def set_resolution(self, value):
        print("Resolution cannot be changed after initialization...")
        pass

    def set_verbose(self, value):
        print("Verbosity cannot be changed after initialization...")
        pass

    def set_data(self, value):
        print("Data cannot be changed after initialization...")
        pass

    def set_is_finished(self, value):
        self.__is_finished = value

    def contains_atom(self, atom):
        if atom.get_nr() == self.get_a_1().get_nr():
            return True
        elif atom.get_nr() == self.get_a_2().get_nr():
            return True
        elif atom.get_nr() == self.get_a_3().get_nr():
            return True
        else:
            return False

    resolution = property(get_resolution, set_resolution, None, "resolution for kde")
    verbose = property(get_verbose, set_verbose, None, "Extend of print messages")
    data = property(get_data, set_data, None, "Data to do kde on.")
    data = property(get_data, set_data, None, "Data to do kde on.")
    is_finished = property(get_is_finished, set_is_finished, None, "Data to do kde on.")
