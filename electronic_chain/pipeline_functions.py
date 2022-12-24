"""Functions used in the pipline execution"""

import os
import importlib
import numpy as np
import math
import electronic_chain.ec_config
from electronic_chain.trace_functions import add_traces, multiply_traces, add_traces_randomized
import grand.io.root_trees

def store_traces(traces_t, tree, copy_tree=None):
    """Stores provided traces in the provided tree"""
    # Copy contents of another tree if requested
    if copy_tree:
        tree.copy_contents(copy_tree)

    # Different traces fields for ADC tree
    if "ADC" in tree.type.upper():
        tree.trace_0 = traces_t[:, 0, :].astype(np.int16)
        tree.trace_1 = traces_t[:, 1, :].astype(np.int16)
        tree.trace_2 = traces_t[:, 2, :].astype(np.int16)
    else:
        tree.trace_x = traces_t[:, 0, :].astype(np.float32)
        tree.trace_y = traces_t[:, 1, :].astype(np.float32)
        tree.trace_z = traces_t[:, 2, :].astype(np.float32)

    tree.fill()

# ToDo: there should also be a kind of tree names list. Probably will be solved with usage of DataFile
def execute_pipeline(pipeline, filelist, output_dir=""):
    """Execute the pipeline"""

    # Changes the returns of functions to dictionaries - needed for the pipeline
    electronic_chain.ec_config.in_pipeline = True

    # Execute the prep_func before looping through files
    # Only one prep func allowed in pipeline
    if "prefileloop_call" in pipeline:
        prep_func = importlib.__import__(pipeline["prefileloop_call"]["module"], fromlist=("prefileloop_call")).prep_func
        var_dict = prep_func(**pipeline["prefileloop_call"]["kwargs"])

    # Loop through files
    for in_root_file in filelist:
        print("************** Analysing file", in_root_file)

        # Removing the trees from the previous file from memory
        # ToDo: This should not be up to a user, at least not in this ugly way
        grand.io.root_trees.grand_tree_list = []

        # Execute all the preeventloop_calls in the pipeline
        for (key,part) in pipeline.items():
            if part["type"]=="preeventloop_call":
                # Find the specified preeventloop function in the module
                func = getattr(importlib.__import__(part["module"]), key)
                # Call the function
                var_dict.update(func(pipeline=pipeline, in_root_file=in_root_file, output_dir=output_dir, **var_dict))

        output_trees = var_dict["output_trees"]
        tshower, tefield = var_dict["tshower"], var_dict["tefield"]

        # Loop through events
        for i in range(tshower.get_entries()):
            tshower.get_entry(i)
            # ToDo: A bug in root_trees? Get entry should not be necessary after get entry on tshower (friends!)
            tefield.get_entry(i)
            # tvoltage.copy_contents(tefield)

            # Loop through the elements of the pipeline
            for (key,part) in pipeline.items():
                print("Applying ", key)
                # Skip the prefileloop_call and preeventloop_call
                if key=="prefileloop_call" or key=="preeventloop_call": continue

                # Take action depending on the type of the pipeline element
                # ToDo: when we upgrade to Python >=3.10, "match" conditional should be used
                if part["type"]=="call":
                    # ToDo: slightly more optimal to do the import before all the looping (but not much)
                    func = getattr(importlib.__import__(part["module"]), key)
                    # Merge the function arguments with the output dictionary
                    # ToDo: It should be dict union "|" for python >=3.9
                    if "kwargs" in part:
                        input_dict = {**part["kwargs"], **var_dict}
                    else:
                        input_dict = var_dict

                    res = func(**input_dict)
                    # Update the results dictionary
                    var_dict.update(res)

                elif part["type"]=="add":
                    res = add_traces(addend=var_dict[key], **var_dict)
                    var_dict.update(res)

                elif part["type"]=="add_randomized":
                    res = add_traces_randomized(addend=var_dict[key], **var_dict)
                    var_dict.update(res)

                elif part["type"]=="multiply":
                    res = multiply_traces(multiplier=var_dict[key], **var_dict)
                    var_dict.update(res)

                # Store the results in a tree if requested
                elif part["type"]=="store":
                    if "copy_tefield" in part and part["copy_tefield"] == True:
                        store_traces(var_dict["traces_t"], part["output_tree"], tefield)
                    else:
                        store_traces(var_dict["traces_t"], part["output_tree"])

        # Write all the trees that are to be written
        for tree in output_trees:
            print("Writing", tree.tree_name)
            tree.write()
