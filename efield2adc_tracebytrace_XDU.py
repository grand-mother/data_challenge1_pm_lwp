#!/usr/bin/python

import grand.io.root_trees
from grand.io.root_trees import *
from electronic_chain import *
from misc_functions import *
import argparse

from XDU_electronic_chain import manual_pipeline
import XDU_electronic_chain.config as config


def main():

    config.XDU_files_path = "XDU_electronic_chain/XDU_files"

    clparser = argparse.ArgumentParser()
    clparser.add_argument("filename", nargs="+")
    clparser.add_argument("-os", "--output_file_suffix", default="_voltage_adc")
    clparser.add_argument("-od", "--output_dir", default="")
    clargs = clparser.parse_args()

    # Loop through files
    for in_root_file in clargs.filename:

        # Removing the trees from the previous file from memory
        # ToDo: This should not be up to a user, at least not in this ugly way
        grand.io.root_trees.grand_tree_list = []

        filename_only = os.path.split(in_root_file)[-1]

        # If given, create the output directory and use it
        if clargs.output_dir != "":
            os.makedirs(clargs.output_dir, exist_ok=True)
            out_root_file = clargs.output_dir+"/"+filename_only.split(".root")[0]+clargs.output_file_suffix+".root"
        else:
            out_root_file = filename_only.split(".root")[0] + clargs.output_file_suffix + ".root"

        # Open the ROOT tree with simulation shower data
        tshower = ShowerEventSimdataTree(in_root_file)
        # Open the ROOT tree with traces
        tefield = EfieldEventTree(in_root_file)
        # Open the ROOT tree with simulation run info
        trunefieldsimdata = EfieldRunSimdataTree(in_root_file)
        trunefieldsimdata.get_entry(0)
        # In this file antennas are misordered. This is the dictionary that gives ROOT file trace number for given trace file number
        # trace_file_to_root = {95: 62, 94: 95, 93: 42, 92: 19, 91: 21, 90: 63, 89: 14, 88: 41, 87: 0, 86: 45, 85: 54, 84: 22, 83: 20, 82: 7, 81: 69, 80: 77, 79: 68, 78: 81, 77: 30, 76: 76, 75: 17, 74: 3, 73: 18, 72: 32, 71: 53, 70: 55, 69: 4, 68: 90, 67: 82, 66: 40, 65: 12, 64: 93, 63: 13, 62: 38, 61: 44, 60: 84, 59: 15, 58: 70, 57: 58, 56: 75, 55: 78, 54: 23, 53: 8, 52: 34, 51: 29, 50: 9, 49: 37, 48: 59, 47: 79, 46: 16, 45: 83, 44: 86, 43: 47, 42: 36, 41: 92, 40: 51, 39: 66, 38: 49, 37: 31, 36: 57, 35: 6, 34: 50, 33: 61, 32: 56, 31: 89, 30: 48, 29: 65, 28: 67, 27: 64, 26: 71, 25: 28, 24: 52, 23: 27, 22: 25, 21: 91, 20: 35, 19: 85, 18: 33, 17: 87, 16: 24, 15: 88, 14: 46, 13: 39, 12: 73, 11: 72, 10: 5, 9: 94, 8: 26, 7: 43, 6: 11, 5: 60, 4: 80, 3: 74, 2: 10, 1: 1, 0: 2}

        # tefield1 = EfieldEventTree(out_root_file)
        tvoltage = VoltageEventTree(out_root_file)
        tadc = ADCEventTree(out_root_file)

        for i in range(tshower.get_entries()):
            tefield.get_entry(i)
            # A bug in root_trees? Get entry should not be necessary
            tshower.get_entry(i)
            tvoltage.copy_contents(tefield)

            # e_theta = float(list1[5])
            # e_phi = float(list1[6])
            e_theta = 180-tshower.shower_zenith
            e_phi = tshower.shower_azimuth-180
            if e_theta<0: e_theta=360-e_theta
            if e_phi<0: e_phi=360+e_phi

            print("theta is:", e_theta, "degree", tshower.shower_azimuth)
            print("phi is:", e_phi, "degree", tshower.shower_zenith)

            # ===================================Arbitrary input file first generates input parameters needed by subroutine================================================
            Ts = trunefieldsimdata.t_bin_size

            tvoltage.trace_x.clear()
            tvoltage.trace_y.clear()
            tvoltage.trace_z.clear()

            time_passed(True)

            for trace_num in range(len(tefield.trace_x)):
                traces_t = np.stack([tefield.trace_x[trace_num], tefield.trace_y[trace_num], tefield.trace_z[trace_num]], axis=-2)

                original_trace_length = traces_t.shape[-1]

                traces_t = adjust_traces(traces_t, Ts)

                all_traces_t, all_traces_f = manual_pipeline(traces_t, e_phi, e_theta, du_count=tefield.du_count, original_trace_length=original_trace_length, sampling_time=0.5, return_all_traces=True)

                tvoltage.trace_x.append(all_traces_t[-2][0,:].astype(np.float32).tolist())
                tvoltage.trace_y.append(all_traces_t[-2][1,:].astype(np.float32).tolist())
                tvoltage.trace_z.append(all_traces_t[-2][2,:].astype(np.float32).tolist())

                tadc.trace_0.append(all_traces_t[-1][0,:].astype(np.int16).tolist())
                tadc.trace_1.append(all_traces_t[-1][1,:].astype(np.int16).tolist())
                tadc.trace_2.append(all_traces_t[-1][2,:].astype(np.int16).tolist())


            print(time_passed())

            tvoltage.fill()
            tadc.fill()
            tadc.write()
            print("Wrote trees")


if __name__ == '__main__':
    main()
