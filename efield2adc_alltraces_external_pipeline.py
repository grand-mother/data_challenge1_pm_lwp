#!/usr/bin/python
# Created by Lech Wiktor Piotrowski

import sys

from electronic_chain import execute_pipeline

import yaml
import argparse
import numpy as np

def main():

    # Parse the command line arguments
    clparser = argparse.ArgumentParser()
    clparser.add_argument("filename", nargs="+")
    clparser.add_argument("-p", "--pipeline_file_name", required=True)
    clparser.add_argument("-od", "--output_dir", default="")
    clparser.add_argument("-se", "--seed", default=None)
    clargs = clparser.parse_args()

    # Load the pipeline from the specified external YAML file into a Python dictionary
    pipeline = yaml.safe_load(open(clargs.pipeline_file_name, "r"))

    # Set the random seed if requested on command line (useful to avoid differences between calls due to changing galactic noise)
    if clargs.seed is not None:
        np.random.seed(int(clargs.seed))

    # Execute the loaded pipeline dictionary
    execute_pipeline(pipeline, clargs.filename, clargs.output_dir)

if __name__ == '__main__':
    main()
