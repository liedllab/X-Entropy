"""
part of the entroPy module

@author: paq
"""
import numpy as np
from entropy.kde_kernel import __kde_kernel as kernel
from .resolution import process_resolution_argument


def center_grid(grid):
    # TODO center grid
    # return np.mean(list(zip(grid, grid[1:])), axis=1)
    return grid


class kde(object):
    __bins = None
    __data = None
    __resolution = None
    __kde_done = False
    __successful = None
    __verbose = None
    __pdf = None
    __pdf_x = None
    __bandwidth = None
    __is_finished = None

    def __init__(self, data, resolution=4096, verbose=False):
        self.__data = data
        self.__resolution = process_resolution_argument(resolution)
        self.__verbose = verbose
        self.__is_finished = False

    def calculate(self, verbose=None, resolution=None):
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
        if not (resolution is None):
            new_res = process_resolution_argument(resolution)
            if verbose:
                print("Using resolution of {}".format(new_res))
            self.set_resolution(new_res)

        if verbose:
            print("Initializing C++ kernel for kde...")
        k = kernel(self.data, self.resolution)
        k.calculate()
        self.set_is_finished(True)
        if verbose:
            print("KDE finished.")
        self.set_pdf_x(center_grid(k.getGrid()))
        self.set_bandwidth(k.getBandwidth())
        self.set_pdf(k.getPDF())

    # Getter #
    def get_resolution(self):
        return self.__resolution

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
    is_finished = property(get_is_finished, set_is_finished, None, "Data to do kde on.")
    bandwidth = property(get_bandwidth, set_bandwidth, None, "bandwidth")
    pdf_x = property(get_pdf_x, set_pdf_x, None, "Pprobability densite sxfvsj_x ")
    pdf = property(get_pdf, set_pdf, None, "Pprobability densite sxfvsj ")

