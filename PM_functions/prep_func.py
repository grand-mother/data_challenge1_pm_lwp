import numpy as np
import math
import electronic_chain.ec_config as ec_config
import PM_functions.config as PM_config
from PM_functions.LNA import LNA_get
from PM_functions.filters import filter_get


def prep_func(**kwargs):
    PM_config.PM_files_path = "PM_functions/PM_files"

    # Use (i)rfft by default
    ret_dict = {"rfft": True, "irfft": True}

    # ToDo: !!!!  Need to be read before! Probably require another type of module - in the file loop. Or some other logic change
    sampling_time = 0.5
    fs = 1 / sampling_time * 1000  # sampling frequency, MHZ
    N = math.ceil(fs)
    freqs_tot = np.fft.rfftfreq(N, d=sampling_time * 1e-9) / 1e6  # MHz

    ret_dict["freqs"] = freqs_tot

    # Read LNA coefficient
    ret_dict.update(LNA_get(freqs_tot, 1))

    # Read cable and filter coefficient
    ret_dict.update(filter_get(freqs_tot, 1))

    return ret_dict