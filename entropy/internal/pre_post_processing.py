import numpy as np


def reshape_arrays_eventually(some_array):
    out = np.squeeze(some_array)
    if out.shape == ():  # special case: single number
        out = out[()]  # weird way of accessing a squeezed single number for np arrays
    return out


def start_end_from_grid(grid):
    # You will need this for non dihedrals...
    return np.nanmin(grid), np.nanmax(grid)


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


def sanity_check_input_data(data, weights=None, weight_switch=True):
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
                      "data: {}, weights: {}".format(data.shape, weights.shape)
            raise ValueError(err_msg)


def process_data_shapes(data, weights=None, weight_switch=True):
    data = np.array(data)
    sanity_check_input_data(data, weights, weight_switch)

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


def process_method_argument(method, available_methods=("simpson", "riemann")):
    """Catches unknown integration methods on a python
    level. Deals with upper and lower case writing.

    """
    if not (method.lower() in available_methods):
        err_msg = "{} is not a valid integration method\nChoose from {}.".format(method, available_methods)
        raise ValueError(err_msg)
    return method.lower().capitalize()  # first letter caps, rest lower case
