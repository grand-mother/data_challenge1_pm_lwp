import os
import math
import numpy as np
from grand.io.root_trees import *
from misc_functions import are_trees_unique
from PM_functions.LNA import LNA_get
from PM_functions.filters import filter_get


def preevent_func(pipeline, in_root_file, output_dir, **kwargs):
    """Called before each event loop
    Assuming this won't change much, thus I put it in this file"""
    ret_dict = {}

    filename_only = os.path.split(in_root_file)[-1]

    # If given, create the output directory and use it
    if output_dir != "":
        os.makedirs(output_dir, exist_ok=True)

    # Open the ROOT tree with simulation shower data
    tshower = ShowerEventSimdataTree(in_root_file)
    ret_dict["tshower"] = tshower
    # Open the ROOT tree with traces
    tefield = EfieldEventTree(in_root_file)
    ret_dict["tefield"] = tefield
    # Open the ROOT tree with simulation run info
    trunefieldsimdata = EfieldRunSimdataTree(in_root_file)
    ret_dict["trunefieldsimdata"] = trunefieldsimdata

    # Generate the list of files/trees to be written
    output_trees = []
    for (key, part) in pipeline.items():
        if part["type"] == "store":
            # Get the tree class object from its name
            class_object = globals()[part["tree_type"]]

            prefix = ""
            if output_dir != "":
                prefix = output_dir + "/"

            # Create the tree instance in a filename with suffix
            if "filename_suffix" in part:
                output_filename = prefix + filename_only.split(".root")[0] + part["filename_suffix"] + ".root"
            # Create the tree instance in the input file (dangerous, can corrupt data in a rare, bad case)
            else:
                output_filename = prefix + filename_only

            part["output_tree"] = class_object(output_filename)

            # Set the tree name if specified
            if "tree_name" in part:
                part["output_tree"].tree_name(part["tree_name"])

            output_trees.append(part["output_tree"])

    # Check if there are no trees with the same name in the same files
    if not are_trees_unique(output_trees):
        raise RuntimeError("Can not have trees with identical names in the same files. Please check the pipeline output trees' file names and tree names!")

    ret_dict["output_trees"] = output_trees

    # Get the time bin size
    trunefieldsimdata.get_entry(0)
    sampling_time = trunefieldsimdata.t_bin_size

    # Adjust traces function in pipeline extends the trace length to 2000
    if "adjust_traces" in pipeline:
        # fs = 1 / sampling_time * 1000  # sampling frequency, MHZ
        # N = math.ceil(fs)
        trace_length = 2000
    else:
        # Get the traces length (assuming it is constant for all events)
        tefield.get_entry(0)
        trace_length = np.array(tefield.trace_x).shape[-1]

    N = trace_length

    freqs_tot = np.fft.rfftfreq(N, d=sampling_time * 1e-9) / 1e6  # MHz

    ret_dict["sampling_time"] = sampling_time
    ret_dict["freqs"] = freqs_tot

    # Read LNA coefficient
    ret_dict.update(LNA_get(freqs_tot, 1))

    # Read cable and filter coefficient
    ret_dict.update(filter_get(freqs_tot, 1))

    return ret_dict