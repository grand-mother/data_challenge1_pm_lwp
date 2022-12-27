#!/usr/bin/python
# Created by Lech Wiktor Piotrowski

from electronic_chain import execute_pipeline

import yaml
import argparse


def main():

    # Parse the command line arguments
    clparser = argparse.ArgumentParser()
    clparser.add_argument("filename", nargs="+")
    clparser.add_argument("-od", "--output_dir", default="")
    clargs = clparser.parse_args()

    # Define the pipeline dictionary
    pipeline = \
        {
            # Preparation function - run before the files loop
            "prefileloop_call": {"type": "prefileloop_call", "kwargs": {"XDU_files_path": "XDU_electronic_chain/XDU_files"}, "module": "XDU_electronic_chain"},
            # Preevent loop functions - called once per file
            "preevent_func": {"type": "preeventloop_call", "module": "XDU_electronic_chain.preevent_func"},
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
            "CEL": {"type": "call", "module": "XDU_electronic_chain.antenna_effective_length"},
            # Convert efield to voltage
            "efield2voltage": {"type": "call", "module": "XDU_electronic_chain.efield2voltage"},
            # Readout/prepare galactic noise
            "gala": {"type": "call", "module": "XDU_electronic_chain.galactic_noise"},
            # Add galactic noise
            "galactic_noise": {"type": "add_randomized", "module": "XDU_electronic_chain.galactic_noise"},
            # Readout/prepare LNA coefficient
            "LNA_get": {"type": "call", "module": "XDU_electronic_chain.LNA"},
            # Add LNA noise
            "LNA_coefficient": {"type": "multiply"},
            # Readout/prepare cable and filter coefficients
            "filter_get": {"type": "call", "module": "XDU_electronic_chain.filters"},
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

    # Execute the pipeline dictionary
    execute_pipeline(pipeline, clargs.filename, clargs.output_dir)

    # Store the pipeline dictionary as a YAML file
    # yaml.safe_dump(pipeline, sort_keys=False, stream=open("XDU_pipeline.yaml", "w"))


if __name__ == '__main__':
    main()
