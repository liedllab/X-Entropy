import warnings
import numpy as np
from .pre_post_processing import start_end_from_grid


def is_power_of_two(val):
    """This function will evaluate whether $val is
    a power of two between 1 and 524288 or not.
    Higher powers are not tested here.

    Parameters
    ----------
    val : numeric

    Returns
    -------
    bool

    """

    pows_of_two = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
                   4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288]
    val = int(val)
    return val in pows_of_two


def next_power_of_two(val):
    """Returns the next higher power of two.

    Parameters
    ----------
    val : numeric

    Returns
    -------
    pow_of_two : int

    """
    return int(2**(np.log(val) // np.log(2) + 1))


def process_resolution_argument(resolution, data, verbose=False):
    """Warns about potentially too high or too low
    values and picks the next higher power of two,
    if it was no power of two initially.

    Parameters
    ----------
    resolution
    data

    Returns
    -------

    """
    if isinstance(resolution, str):
        resolution = resolution_from_rule_of_thumb(resolution, data, verbose=verbose)
        # will be checked for whether it is a power of two or not below
    elif isinstance(resolution, int):
        pass
    elif isinstance(resolution, float):
        print("Resolution is not of type int. Trying to cast it to int...")
        resolution = int(resolution)
    else:
        err_msg = "Cannot interpret given argument for resolution:\n{}\n" \
                  "Please give either a single integer or a string.".format(resolution)
        raise ValueError(err_msg)
    if resolution < 100:
        warn_msg = "You are using a rather small resolution. " \
                   "This may potentially lead to inaccurate results..."
        warnings.warn(warn_msg, RuntimeWarning)
    elif resolution > 10000:
        warn_msg = "You are using a rather large resolution. " \
                   "Amongst other things, this may potentially lead to very long runtimes " \
                   "without necessarily improving the accuracy of the result..."
        warnings.warn(warn_msg, RuntimeWarning)

    if not is_power_of_two(resolution):
        resolution = next_power_of_two(resolution)
    return resolution


def resolution_from_rule_of_thumb(resolution, data, verbose=False):
    rules_of_thumb = {"auto": minim_or_sqrt,
                      "freedman_diaconis": freedman_diaconis, "fd": freedman_diaconis,
                      "sturges": sturges,
                      "doane": doane,
                      "sqrt": square_root_choice,
                      "scott": scott}
    resolution = resolution.lower()
    resolution = resolution.replace(" ", "_")

    if not (resolution in list(rules_of_thumb.keys())):
        err_msg = "Cannot interpret given argument for resolution. " \
                  "Give either an integer, or choose of the following:\n{}".format(rules_of_thumb.keys())
        raise ValueError(err_msg)
    data = np.squeeze(data)  # you need to do this, because otherwise, you will have issues with single data sets...
    squeezed_shape = data.shape
    if len(squeezed_shape) == 1:
        return rules_of_thumb[resolution](data)
    elif len(squeezed_shape) == 2:
        if verbose:
            print("Found multiple data sets. Applying rule of thumb on all, and take the maximum resolution estimated.")
        return np.max([rules_of_thumb[resolution](dat) for dat in data])
    else:  # bad paq. you really should handle this properly...
        print("Suspicious data shape...")
        return 4096


def minim_or_sqrt(data, minim=16):
    return np.max([minim, square_root_choice(data)])


def interquartiles(data):
    data = np.sort(data)
    data_quarts = np.array_split(data, 4)
    return data_quarts[1][0], data_quarts[-2][-1]


def scott(data):
    bin_edges = np.histogram_bin_edges(data, bins="scott")
    return len(bin_edges) - 1


def freedman_diaconis(data):
    bin_edges = np.histogram_bin_edges(data, bins="fd")
    return len(bin_edges)-1


def square_root_choice(data):
    return int(np.ceil(np.sqrt(len(data))))


def sturges(data):
    bin_edges = np.histogram_bin_edges(data, bins="sturges")
    return len(bin_edges)-1


def doane(data):
    bin_edges = np.histogram_bin_edges(data, bins="doane")
    return len(bin_edges)-1


# below we have legacy functions, which we do not use
def silverman(data):
    """This is a legacy function"""
    n_dat = len(data)
    iqr = np.diff(interquartiles(data))[0]
    either_or = np.min([np.std(data), iqr / 1.34])
    return 0.9 * either_or * n_dat ** (-1 / 5)


def res_from_silverman(data):
    """This is a legacy function"""
    start, end = start_end_from_grid(data)
    data_range = end-start
    predicted_bandw = silverman(data)
    return data_range / predicted_bandw
