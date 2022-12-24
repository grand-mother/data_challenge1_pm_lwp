#!/usr/bin/python

import sys

from electronic_chain import execute_pipeline

import yaml
import argparse

def main():

    clparser = argparse.ArgumentParser()
    clparser.add_argument("filename", nargs="+")
    clparser.add_argument("-od", "--output_dir", default="")
    clargs = clparser.parse_args()

    # XDU_config.XDU_files_path = "electronic_chain/XDU_electronic_chain/XDU_files"

    pipeline = \
        {
            # Preparation function - run before the files loop
            "prefileloop_call": {"type": "prefileloop_call", "kwargs": {"PM_files_path": "PM_functions/PM_files"}, "module": "PM_functions"},
            # Preevent loop functions - called once per file
            "preevent_func": {"type": "preeventloop_call", "module": "PM_functions.preevent_func"},
            # Stuff below is repeated for every file and event
            # Read shower angles from the shower tree
            "read_angles": {"type": "call", "module": "electronic_chain", "kwargs": {"tree": "tshower"}},
            # Read sampling time (time bin) from the simulation run tree
            "read_sampling_time": {"type": "call", "module": "electronic_chain", "kwargs": {"tree": "trunefieldsimdata"}},
            # Read event traces
            "read_traces": {"type": "call", "module": "electronic_chain", "kwargs": {"tree": "tefield"}},
            # Adjust traces sizes - XDU chain needs this stuff
            "adjust_traces": {"type": "call", "module": "electronic_chain"},
            # Calculate effective length of antenna (needs angles, so for every event)
            # "CEL": {"type": "call", "module": "XDU_electronic_chain.antenna_effective_length"},
            # Convert efield to voltage
            "efield2voltage_pm": {"type": "call", "module": "PM_functions"},
            # Readout/prepare galactic noise
            "generate_galacticnoise": {"type": "call", "module": "PM_functions.noisemodel"},
            # Add galactic noise
            "galactic_noise": {"type": "add_randomized"},
            # Readout/prepare LNA coefficient
            # "LNA_get": {"type": "call", "module": "XDU_electronic_chain.LNA"},
            # Add LNA noise
            "LNA_coefficient": {"type": "multiply"},
            # Readout/prepare cable and filter coefficients
            # "filter_get": {"type": "call", "module": "XDU_electronic_chain.filters"},
            # Add cable noise
            "cable_coefficient": {"type": "multiply"},
            # Add filter noise
            "filter_coefficient": {"type": "multiply"},
            # Turn the strange XDU traces lengths into the original lengths
            "restore_traces_length": {"type": "call", "module": "electronic_chain"},
            # Store the voltage at this stage of the pipeline
            "store_voltage": {"type": "store", "tree_type": "VoltageEventTree", "filename_suffix": "_voltage", "copy_tefield": True},
            # Compute ADC counts from voltage at this stage of the pipeline
            "voltage2adc": {"type": "call", "module": "XDU_electronic_chain"},
            # Store the ADC counts
            "store_adc": {"type": "store", "tree_type": "ADCEventTree", "filename_suffix": "_adc", "copy_tefield": True}
        }

    execute_pipeline(pipeline, clargs.filename, clargs.output_dir)

    # yaml.safe_dump(pipeline, sort_keys=False, stream=open("PM_pipeline.yaml", "w"))


if __name__ == '__main__':
    main()
