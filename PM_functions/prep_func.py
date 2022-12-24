import numpy as np
import math
import electronic_chain.ec_config as ec_config
import PM_functions.config as PM_config


def prep_func(**kwargs):
    PM_config.PM_files_path = "PM_functions/PM_files"

    # Use (i)rfft by default
    ret_dict = {"rfft": True, "irfft": True}

    return ret_dict