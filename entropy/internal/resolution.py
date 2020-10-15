import warnings
import numpy as np


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


def process_resolution_argument(resolution):
    """Warns about potentially too high or too low
    values and picks the next higher power of two,
    if it was no power of two initially.

    Parameters
    ----------
    resolution

    Returns
    -------

    """
    if not isinstance(resolution, int):
        print("Resolution is not of type int. Trying to cast it to int...")
        resolution = int(resolution)
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