import numpy as np
import warnings

from .constants import PI, TWOPI


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


def deg_to_rad(dat_deg):
    return dat_deg/180*PI


def rad_to_deg(dat_rad):
    return dat_rad*180/PI


def sanity_check_dihedral_units(diheds):
    if np.max(diheds) > TWOPI or np.min(diheds) < -TWOPI:
        warn_msg = "X-Entropy detected dihedral values that are much larger than pi. " \
                   "Your values are assumed to haven been given in degree and not radian. " \
                   "X-Entropy expects dihedrals in radian and will transform the unit of your input. " \
                   "Did you give dihedrals in degree?"
        warnings.warn(warn_msg, RuntimeWarning)
        diheds = deg_to_rad(diheds)
    return diheds


def preprocess_dihedral_data(diheds, weights, weight_switch):
    # mirror data
    # # <- for multiple data sets <- # #
    # diheds_out = []
    # for dihed in diheds:
    #     diheds_out.append(np.concatenate([dihed-360, dihed, dihed+360]))
    # diheds = np.array(diheds_out)
    # if weight_switch:  #
    #     weights_out = []
    #     for weight in weights:
    #         weights_out.append(np.concatenate([weight, weight, weight]))
    #     weights = np.array(weights_out)
    diheds = sanity_check_dihedral_units(diheds)
    diheds = np.concatenate([diheds - TWOPI, diheds, diheds + TWOPI])
    if weight_switch:  #
        weights = np.concatenate([weights, weights, weights])

    return diheds, weights


def postprocess_dihedral_pdf(pdf, pdf_x, norm_for_mirrored_data=1/3):
    # mirror data
    lower_idx = np.argmin(np.abs(pdf_x + PI))
    if pdf_x[lower_idx] < -PI:
        lower_idx = lower_idx + 1
    upper_idx = np.argmin(np.abs(pdf_x - PI))
    if pdf_x[upper_idx] > PI:
        upper_idx = upper_idx - 1

    pdf_out = pdf[lower_idx:upper_idx]
    pdf_out /= norm_for_mirrored_data  # we normalized in the "mirrored data", which is 3 times the actual data
    pdf_x_out = pdf_x[lower_idx:upper_idx]
    assert len(pdf_out) == len(pdf_x_out)
    return pdf_out, pdf_x_out


def sanity_check_input_data(data, weights=None, weight_switch=True):
    if False in np.isfinite(data):
        err_msg = "Non-finite values in data!"
        raise ValueError(err_msg)

    if weight_switch:
        weights = np.array(weights)
        if False in np.isfinite(weights):
            err_msg = "Non-finite values in weights!"
            raise ValueError(err_msg)
        if data.shape[-1] != weights.shape[-1]:
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
            weights = np.array(weights)
            # this is as good as pass...
    elif len(data.shape) == 2:
        if data.shape[0] > data.shape[1]:
            warn_msg = "Your input data has suspicious shape. First dim is expected to be smaller than second dim. " \
                       "Data should have shape: (N_features, N_frames)."
            warnings.warn(warn_msg, RuntimeWarning)  # potentially transposed data
        pass  # This is all good
    else:  # This 3 dim or higher...
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
