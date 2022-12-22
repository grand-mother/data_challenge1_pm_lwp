#!/usr/bin/python

import sys

from electronic_chain import execute_pipeline

import yaml
import argparse

def main():

    clparser = argparse.ArgumentParser()
    clparser.add_argument("filename", nargs="+")
    clparser.add_argument("-p", "--pipeline_file_name", required=True)
    clparser.add_argument("-od", "--output_dir", default="")
    clargs = clparser.parse_args()

    # XDU_config.XDU_files_path = "electronic_chain/XDU_electronic_chain/XDU_files"

    pipeline = yaml.safe_load(open(clargs.pipeline_file_name, "r"))

    execute_pipeline(pipeline, clargs.filename, clargs.output_dir)

if __name__ == '__main__':
    main()
